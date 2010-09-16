#lang racket/base

(provide stripe)

;(stripe id element-count partitions)
(define (stripe i s1 s2)
  (define-values (chunk-size rem) (quotient/remainder s1 s2))
;  (define chunk-size (quotient s1 s2))
;  (define rem (remainder s1 s2))
  (define offset (if (i . < . rem) i rem))
  (define start (+ (* i chunk-size) offset))
  (define end (+ (- (* (+ i 1) chunk-size) 1) offset))
  (values start end))

