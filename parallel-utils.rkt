#lang racket/base
(require racket/match)
(require racket/place)
(require racket/place-utils)
(require (for-syntax racket/base))

(provide CGspawn CG-n0-only CG-B CGfor CGid CGSingle CGforS CGSerial p-range CGforSS)

(define-syntax-rule (!or= x b ...)
  (not (ormap (lambda (y) (= x y)) (list b ...))))

(define (stripe-range cg i st en step np)
  (when  (!or= step 1 -1)
    (error "step must be 1 or -1, not ~a" step))
  (define ec (abs (- en st)))
  (define step-sign (if (negative? step) - +))
  (define-values (chunk-size rem) (quotient/remainder ec np))
  (define-values (soff eoff) (if (i . < . rem) (values i (+ i 1)) (values rem rem)))
  (define start (step-sign st (+ (* i chunk-size) soff)))
  (define end (step-sign st (+ (* (+ i 1) chunk-size) eoff)))
;;;  (printf "SR ~a ~a ~a ~a ~a ~a ~a\n" st en step i np start end)
;;;  (flush-output)
  (in-range start end step))

(define-struct CG (id np ch))

(define (CGid cg) (match cg [(CG id _ _) id]))

(define (CG-0-send cg)
  (match cg
    [(CG 0  np pls) (for ([ch pls]) (place-channel-send ch 0))]
    [(CG id np ch ) (place-channel-recv ch)])-)

(define (CG-0-recv cg)
  (match cg
    [(CG 0  np pls) (for ([ch pls]) (place-channel-recv ch))]
    [(CG id np ch ) (place-channel-send ch 1)]))

(define-syntax-rule (when0 np body ...)
   (unless (and (pair? np) (not (= (car np) 0)))
    body ...))

(define (CG-B cg)
  (match cg
    [(CG _ 0 _) (void)]
    [(CG id _ _) 
      (CG-0-recv cg)
      (CG-0-send cg)]))

(define (CGSingle) (make-CG 0 0 #f))


(define-syntax-rule (CGspawn NPE func args ...)
  (match NPE
    [0 (func (CGSingle) args ...)]
    [np 
      (define pls (for/list ([i (in-range 1 np)])
        (place/main (pwkr ch)
          (match (place-channel-recv ch)
            [(list id np args ...)
              (func (make-CG id np ch) args ...)]))))

      (for ([i (in-range 1 np)]
            [ch pls])
        (place-channel-send ch (list i np args ...))) 

      (func (make-CG  0 np pls) args ...)]))

(define-syntax-rule (CG-n0-only cg body ...)
  (match cg
    [(CG _ 0 _) body ...]
    [(CG id _ _) 
      (CG-0-recv cg)
      (when (= id 0)
        body ...)
      (CG-0-send cg)]))

(define-syntax (p-range stx)
  (syntax-case stx (in-range)
    [(_ cg (in-range st en step)) #'(match cg [(CG _ 0 _) (in-range st en step)]
                                              [(CG id np _) (stripe-range cg id st en step np)])]
    [(_ cg (in-range en)) #'(p-range cg (in-range 0 en 1))]
    [(_ cg (in-range st en)) #'(p-range cg (in-range st en 1))]))

(define-syntax-rule (CGfor cg ([V (in-range ST EN)]) body ...)
  (match cg
    [(CG _ 0 _) (for ([V (in-range ST EN)]) body ...)]
    [(CG id np _) 
      (for ([V (stripe-range cg id ST EN 1 np)])
        body ...)]))

(define-syntax-rule (CGSerial cg body ...)
  (match cg
    [(CG _ 0 _) body ...]
    [(CG (and 0 id) np pls)
     body ... 
     (for ([ch pls]) (place-channel-send ch 0)
                     (place-channel-recv ch))]
    [(CG id np ch) 
      (place-channel-recv ch)
      body ...
      (place-channel-send ch 1)]))

(define-syntax-rule (CGforS cg ([V (in-range ST EN)]) body ...)
  (match cg
    [(CG _ 0 _) (for ([V (in-range ST EN)]) body ...)]
    [(CG (and 0 id) np pls) 
      (for ([V (stripe-range cg id ST EN 1 np)]) body ...)
      (for ([ch pls]) (place-channel-send ch 0)
                      (place-channel-recv ch))]
    [(CG id np ch) 
      (place-channel-recv ch)
      (for ([V (stripe-range cg id ST EN 1 np)]) body ...)
      (place-channel-send ch 1)]))

(define-syntax-rule (CGforSS cg ([V R]) body ...)
  (match cg
    [(CG _ 0 _) (for ([V R]) body ...)]
    [(CG (and 0 id) np pls) 
      (for ([V R]) body ...)
      (for ([ch pls]) (place-channel-send ch 0)
                      (place-channel-recv ch))]
    [(CG id np ch) 
      (place-channel-recv ch)
      (for ([V R]) body ...)
      (place-channel-send ch 1)]))
