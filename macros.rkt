#lang racket/base
(require (for-syntax racket/base))
(require (only-in racket/flonum flvector-length [flvector-ref fr] [flvector-set! f!]))
(require (rename-in scheme/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-flvector-ref flvr] 
                    [unsafe-flvector-set! flvs!]
                    [unsafe-fl+ fl+op]
                    [unsafe-fl- fl-op]
                    [unsafe-fl* fl*op]
                    [unsafe-fl/ fl/]
                    [unsafe-fx+ fx+op]
                    [unsafe-fx- fx-op]
                    [unsafe-fx* fx*op]
                    [unsafe-fx= fx=]
                    [unsafe-flmax flmax]
))
(provide define-syntax-case
         PICK2M
         PICK3M
         with-syntax-values
         REPDET*
         PICK2
         PICK3
         <epsilon-vmap

        fx++ fx-- fx+ fx- fx* fl+ fl- fl* fl/ f!+ f!- f!* f!/ flmax flmax* fx=)

(define-syntax-rule (define-syntax-case (N a ...) b ...)
  (define-syntax (N stx)
    (syntax-case stx ()
      [(_ a ...) b ...])))

(define-syntax-rule (with-syntax-values ([(a ...) b] ...) body ...)
(syntax-case (list b ...) ()
  [( (a ...) ...) (let () body ...)]))

(define-syntax-case (REPDET* R2 V1 V2 V3 VV1 VV2 VV3)
#'(case (syntax->datum R2)
  [(1) #'(VV1 V2 V3)]
  [(2) #'(V1 VV2 V3)]
  [(3) #'(V1 V2 VV3)]
  [else (raise (format "invalid REPDEP value ~a" (syntax->datum #'R2)))]))

(define-syntax-case (PICK2 R2 V1 V2)
#'(begin
  (case (syntax->datum R2) [(1) #'V1] [(2) #'V2]
  [else (raise (format "~a" (syntax->datum R2)))])))

(define-syntax-case (PICK3 R2 V1 V2 V3)
#'(begin
  (case (syntax->datum R2) [(1) #'V1] [(2) #'V2] [(3) #'V3]
  [else (raise (format "~a" (syntax->datum R2)))])))

(define-syntax-case (PICK2M R2 V1 V2)
  (case (syntax->datum #'R2) [(1) #'V1] [(2) #'V2]
  [else (raise (format "~a" (syntax->datum #'R2)))]))

(define-syntax-case (PICK3M R2 V1 V2 V3)
  (case (syntax->datum #'R2) [(1) #'V1] [(2) #'V2] [(3) #'V3]
  [else (raise (format "~a" (syntax->datum #'R2)))]))

(define-syntax-rule (<epsilon-vmap v epsilon)
  (let ([vmax (flvector-length v)])
    (let loop ([t #t]
               [i 0])
      (if (and t (i . < . vmax))
        (loop (and t ((fr v i) . < . epsilon)) (add1 i))
        t))))

(define-syntax-rule (fx++ a) (fx+ a 1))
(define-syntax-rule (fx-- a) (fx- a 1))
(define-syntax (fx+ stx)
  (syntax-case stx ()
    [(_ a) #'a]
    [(_ a b) #'(fx+op a b)]
    [(_ a b c ...) #'(fx+op (fx+op a b) (fx+ c ...))]))

(define-syntax (fx- stx)
  (syntax-case stx ()
    [(_ a) #'(fx-op 0 a)]
    [(_ a b) #'(fx-op a b)]
    [(_ a b c ...) #'(fx- (fx-op a b) c ...)]))

(define-syntax (fx* stx)
  (syntax-case stx ()
    [(_ a) #'a]
    [(_ a b) #'(fx*op a b)]
    [(_ a b c ...) #'(fx*op (fx*op a b) (fx* c ...))]))


(define-syntax (fl+ stx)
  (syntax-case stx ()
    [(_ a) #'a]
    [(_ a b) #'(fl+op a b)]
    [(_ a b c ...) #'(fl+op (fl+op a b) (fl+ c ...))]))

(define-syntax (fl- stx)
  (syntax-case stx ()
    [(_ a) #'(fl-op 0.0 a)]
    [(_ a b) #'(fl-op a b)]
    [(_ a b c ...) #'(fl- (fl-op a b) c ...)]))

(define-syntax (fl* stx)
  (syntax-case stx ()
    [(_ a) #'a]
    [(_ a b) #'(fl*op a b)]
    [(_ a b c ...) #'(fl*op (fl*op a b) (fl* c ...))]))

(define-syntax-rule (f!+ v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (fl+ (fr v idx) val ...))))
(define-syntax-rule (f!- v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (fl- (fr v idx) val ...))))
(define-syntax-rule (f!* v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (fl* (fr v idx) val ...))))
(define-syntax-rule (f!/ v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (fl/ (fr v idx) val ...))))
(define-syntax (flmax* stx)
  (syntax-case stx ()
    [(_ a b) #'(flmax a b)]
    [(_ a b ...) #'(flmax a (flmax* b ...))]))
