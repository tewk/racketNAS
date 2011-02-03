#lang racket/base
(require racket/list)
(require racket/match)
(require racket/place)
(require (for-syntax racket/base))
(require (only-in racket/unsafe/ops [unsafe-fx+ fx+]
                                    [unsafe-fx- fx-]
                                    [unsafe-fx* fx*]
                                    [unsafe-fx= fx=]
                                    [unsafe-fx< fx<]
                                    [unsafe-fxabs fxabs]))

(provide CGspawn CG-n0-only CG-B CGfor CGfor/fold CGid CGnp CGSingle CGSerial CGpipeline CG-Parallel-Only CGfor/stride)

(define-syntax-rule (define-syntax-case (N a ...) b ...)
  (define-syntax (N stx)
    (syntax-case stx ()
      [(_ a ...) b ...])))

(define-syntax-rule (define-syntax-cases (n [(a ...) b ...] ...))
     (define-syntax(n stx)
       (syntax-case stx ()
         [(a ...) b ...] ...)))

;[genname (datum->syntax #f (gensym))]
(define-syntax-case (place/anon (ch) body ...)
 (with-syntax ([interal-def-name
                (syntax-local-lift-expression #'(lambda (ch) body ...))]
               [funcname #'OBSCURE_FUNC_NAME_%#%])
  (syntax-local-lift-provide #'(rename interal-def-name funcname))
  #'(let ([module-path (resolved-module-path-name
          (variable-reference->resolved-module-path
           (#%variable-reference)))])
   (place module-path (quote funcname)))))

(define-syntax-rule (!or= x b ...)
  (not (ormap (lambda (y) (fx= x y)) (list b ...))))

(define (stripe cg i st en step np)
  (when  (!or= step 1 -1)
    (error "step must be 1 or -1, not ~a" step))
  (define ec (fxabs (fx- en st)))
  (define step-sign (if (negative? step) fx- fx+))
  (define-values (chunk-size rem) (quotient/remainder ec np))
  (define-values (soff eoff) (if (i . fx< . rem) (values i (fx+ i 1)) (values rem rem)))
  (define start (step-sign st (fx+ (fx* i chunk-size) soff)))
  (define end (step-sign st (fx+ (fx* (fx+ i 1) chunk-size) eoff)))
;  (printf "SR ~a ~a ~a ~a ~a ~a ~a\n" st en step i np start end)
;  (flush-output)
;  (in-range start end step))
  (values start end step))

(define-struct CG (id np pls))

(define (CGid cg) (CG-id cg))
(define (CGnp cg) (CG-np cg))

(define (CG-0-send cg)
  (match cg
    [(CG 0  np (list-rest _ pls)) (for ([ch pls]) (place-channel-send ch 0))]
    [(CG id np (list-rest ch  _)) (place-channel-recv ch)]))

(define (CG-0-recv cg)
  (match cg
    [(CG 0  np (list-rest _ pls)) (for ([ch pls]) (place-channel-recv ch))]
    [(CG id np (list-rest ch  _)) (place-channel-send ch 1)]))

(define (CG-B cg)
  (match cg
    [(CG _ 0 _) (void)]
    [else
      (CG-0-recv cg)
      (CG-0-send cg)]))

(define (CGSingle) (make-CG 0 0 #f))


(define-syntax-rule (CGspawn NPE func args ...)
  (match NPE
    [0 (func (CGSingle) args ...)]
    [np 
      (define pls (for/list ([i (in-range 1 np)])
        (place/anon (ch)
          (match (place-channel-recv ch)
            [(list id np pls args ...)
              (func (make-CG id np (cons ch pls)) args ...)]))))

      (for ([i (in-range 1 np)]
            [ch pls])
        (place-channel-send ch (list i np pls args ...))) 

      (func (make-CG 0 np (cons #f pls)) args ...)]))

(define-syntax-rule (CG-n0-only cg body ...)
  (match cg
    [(CG _ 0 _) body ...]
    [(CG id _ _) 
      (CG-0-recv cg)
      (when (= id 0)
        body ...)
      (CG-0-send cg)]))

(define-syntax-cases (CGfor
  [(_ cg ([V (in-range st en step)]) body ...) 
    #'(match cg [(CG _ 0 _) (for ([V (in-range st en step)]) body ...)]
                [(CG id np _) 
                 (let-values ([(ST EN STEP) (stripe cg id st en step np)])
                  (for ([V (in-range ST EN STEP)]) body ...))])]
  [(_ cg ([V (in-range st en)])      body ...) #'(CGfor cg ([V (in-range st en 1)]) body ...)]
  [(_ cg ([V (in-range en)])         body ...) #'(CGfor cg ([V (in-range  0 en 1)]) body ...)]))

(define-syntax-cases (CGfor/fold
  [(_ cg (VARS ...) ([V (in-range st en step)]) body ...) 
    #'(match cg [(CG _ 0 _) (for/fold (VARS ...) ([V (in-range st en step)]) body ...)]
                [(CG id np _) 
                 (let-values ([(ST EN STEP) (stripe cg id st en step np)])
                  (for/fold (VARS ...) ([V (in-range ST EN STEP)]) body ...))])]
  [(_ cg (VARS ...) ([V (in-range st en)])      body ...) #'(CGfor/fold cg (VARS ...) ([V (in-range st en 1)]) body ...)]
  [(_ cg (VARS ...) ([V (in-range en)])         body ...) #'(CGfor/fold cg (VARS ...) ([V (in-range  0 en 1)]) body ...)]))

(define-syntax-cases (CGfor/stride
  [(_ cg ([V (in-range st en step)]) body ...) 
    #'(match cg [(CG _ 0 _) (for ([V (in-range st en step)]) body ...)]
                [(CG id np _) (for ([V (in-range (fx+ st id) en (fx* step np))]) body ...)])]
  [(_ cg ([V (in-range st en)])      body ...) #'(CGfor cg ([V (in-range st en 1)]) body ...)]
  [(_ cg ([V (in-range en)])         body ...) #'(CGfor cg ([V (in-range  0 en 1)]) body ...)]))

(define-syntax-rule (ASSERT A B)
  (when (not (= A B))
    (eprintf "PIPELINE ERROR\n")))

(define-syntax-rule (CGpipeline cg k body ...)
  (match cg
    [(CG _ 0 _) body ...]
    [(CG id np pls) 
;      (unless (and (= k 10) (= id 0))
;        (define idx (if (= id 0) (sub1 np) 0))
;        (match (place-channel-recv (list-ref pls idx))
      (unless (= id 0)
        (match (place-channel-recv (car pls))
          [(list sid sk srid) (void)])) ;(printf "ME: ~a S: ~a K: ~a SK:~a ~a\n" id sid k sk srid)]))
      body ...
;      (unless (and (= k 1) (= id (sub1 np)))
;        (define idx (if (= id (sub1 np)) 0 (add1 id))) 
;        (place-channel-send (list-ref pls idx) (list id k idx)))]))
      (unless (= id (sub1 np)) (place-channel-send (list-ref pls (add1 id)) (list id k (add1 id))))]))

(define-syntax-rule (CGSerial cg body ...)
  (match cg
    [(CG _ 0 _) body ...]
    [(CG (and 0 id) np (list-rest _ pls))
     body ... 
     (for ([ch pls]) (place-channel-send ch 2)
                     (place-channel-recv ch))]
    [(CG id np (list-rest ch _)) 
      (place-channel-recv ch)
      body ...
       (place-channel-send ch 3)]))

(define-syntax-rule (CG-Parallel-Only cg body ...)
  (match cg
    [(CG _ 0 _) (void)]
    [else body ...])) 
