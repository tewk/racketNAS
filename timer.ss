#lang scheme
(provide read-timer
         timer-start 
         timer-stop) 
(define timers-last-start (make-vector 30 0.0))
(define running-total (make-vector 30 0.0))

(define (read-timer arg)
  (vector-ref running-total arg))

(define (timer-start index)
  (vector-set! timers-last-start index (current-inexact-milliseconds)))

(define (timer-stop index) 
  (let* ([stop (current-inexact-milliseconds)]
         [diff (- stop (vector-ref timers-last-start index))])
    (vector-set! running-total index
                 (+ (vector-ref running-total index) diff))))