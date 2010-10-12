#lang racket/base
(require (for-syntax racket/base))
(provide define-syntax-case
         PICK2M
         PICK3M
         with-syntax-values
         REPDET*
         PICK2
         PICK3)

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

