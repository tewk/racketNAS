#lang racket
(require tests/eli-tester
         "places-mpi.rkt")

;(for ([n (in-range 9)])
;  (displayln (test-build-hypercube-channels n)))

(test
  (for/list ([n (in-range 9)])
    (test-build-hypercube-channels n))
  => '(
#() 
#(()) 
#((1) (0)) 
#((1 2) (0) (0)) 
#((1 2) (0) (3 0) (2))
#((1 2 4) (0) (3 0) (2) (0))
#((1 2 4) (0) (3 0) (2) (5 0) (4))
#((1 2 4) (0) (3 0) (2) (5 6 0) (4) (4))
#((1 2 4) (0) (3 0) (2) (5 6 0) (4) (7 4) (6))))

(define (splat txt fn)
  (call-with-output-file fn #:exists 'replace
      (lambda (out)
        (fprintf out "~a" txt))))
(define (places-wait places) (for ([pl places]) (place-wait pl)))

(define (scatter-test)
  (splat
  #<<END
  (module pct1 racket
    (require "places-mpi.rkt")
    (provide place-main)

    (define (place-main ch)
      (define comgrp (receive-hypercube-comgrp ch))
      (let ([v (scatter comgrp)])
        (place-channel-send ch (vector-append v (vector (comgrp-id comgrp)))))
      (let ([v (broadcast comgrp)])
        (place-channel-send ch (comgrp-append-id comgrp v))))
  )
END
  "pct1.ss")

  (define places (for/list ([i (in-range 7)]) (place "pct1.ss" 'place-main)))
  (define comgrp (build-hypercube-comgrp 8 places))

  (test (scatter comgrp '#(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)) => '#(0 1))
  (for ([i (in-range 1 8)]
        [pl places])
    (let ([v (place-channel-recv pl)])
      (test v => (vector (* i 2) (+ ( * i 2) 1) i))))

  (test (broadcast comgrp "foobar") => "foobar")
  (for ([i (in-range 1 8)]
        [pl places])
    (let ([v (place-channel-recv pl)])
      (test v => (string-append "foobar" (number->string i)))))

  (places-wait places)
)

(define (reduce-test)
  (splat
  #<<END
  (module pct1 racket
    (require "places-mpi.rkt")
    (provide place-main)

    (define (place-main ch)
      (define comgrp (receive-hypercube-comgrp ch))
      (let ([v (reduce/vector comgrp + 0 (vector (comgrp-id comgrp) 1 2 3 (comgrp-id comgrp)))])
        (place-channel-send ch v)))
  )
END
  "pct1.ss")

  (define places (for/list ([i (in-range 7)]) (place "pct1.ss" 'place-main)))
  (define comgrp (build-hypercube-comgrp 8 places))

  (test (reduce/vector comgrp + 0 (vector (comgrp-id comgrp) 1 2 3 (comgrp-id comgrp))) => #(28 8 16 24 28))
  
  (for ([i (list #(1 1 2 3 1) #(5 2 4 6 5) #(3 1 2 3 3) #(22 4 8 12 22) #(5 1 2 3 5) #(13 2 4 6 13) #(7 1 2 3 7))]
        [pl places])
    (let ([v (place-channel-recv pl)])
      (test i => v)))


  (places-wait places)
)

(scatter-test)
(reduce-test)
