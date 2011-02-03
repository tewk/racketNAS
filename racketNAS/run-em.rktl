#lang racket/base

(require racket/system)

(define benchmarks (list "IS" "FT" "CG" "MG" "SP" "BT" "LU"))
(define classes (list "S" "W" "A" "B" "C"))
(define class "S")

(for ([b benchmarks])
  (define cmd (format "rktl -tm ~a.rkt CLASS=~a ~a" (string-downcase b) class "SERIAL"))
  (printf "~a\n" cmd)
  (system cmd))
