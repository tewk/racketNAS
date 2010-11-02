#lang racket/base
(require racket/flonum)
(provide pflv pflve pflvi)

(define vr vector-ref)
(define flvr flvector-ref)

(define (pflvi v d)
  (for ([i (in-range (vector-length v))])
    (printf "~a ~a ~a\n" d i (vr v i)))
  (flush-output))

(define (pflve v d)
  (for ([i (in-range (flvector-length v))])
    (printf "~a ~a ~a\n" d i (flvr v i)))
  (flush-output)
  (exit 0))

(define (pflv v d)
  (for ([i (in-range (flvector-length v))])
    (printf "~a ~a ~a\n" d i (flvr v i)))
  (flush-output))
 
