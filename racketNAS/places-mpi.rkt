#lang racket 

(provide (struct-out comgrp)
         build-hypercube-comgrp
         build-hypercube-channels
         test-build-hypercube-channels
         receive-hypercube-comgrp
         scatter
         broadcast
         reduce/vector
         comgrp-append-id comgrp)

(define-struct comgrp (id cnt peers))

(define (log-2-ceiling n)
  (let loop ([p 0]
             [m 1])
    (if (m . > . n)
        p
        (loop (add1 p) (arithmetic-shift m 1)))))

;;id : nonzero-positive-integer
(define (subtree-size id)
  (when (id . <= . 0) (raise "id must be > 0"))
  (let loop ([p 1]
             [m 1])
    (if (bitwise-and m id)
        p
        (loop (add1 p) (arithmetic-shift m 1)))))

(define-syntax-rule (vsa! v idx val)
  (vector-set! v idx (append (vector-ref v idx) (list val))))

(define-syntax-rule (buildx chs j peer-id)
  (let-values ([(ch1 ch2) (place-channel)])
    (vsa! chs j (list peer-id ch1))
    (vsa! chs peer-id (list j ch2))))

(define-syntax-rule (testx chs j peer-id)
  (begin
    (vsa! chs j peer-id)
    (vsa! chs peer-id j)))

(define-syntax-rule (build-hypercube-channels n)
  (build-hypercube-channels-worker buildx n))

(define-syntax-rule (test-build-hypercube-channels n)
  (build-hypercube-channels-worker testx n))

(define-syntax-rule (build-hypercube-channels-worker builder n)
  (let ([chs (make-vector n '())])
    (let loop ([two_2_i 1])
      (when (two_2_i . < . n)
        (for ([j (in-range 0 n 2)])
          (define peer-id (bitwise-xor two_2_i j))
          (when (and (peer-id . < . n)
                     (zero? (bitwise-and peer-id (sub1 two_2_i)))
                     (not (zero? (bitwise-and peer-id two_2_i))))
              (builder chs j peer-id)))
        (loop (arithmetic-shift two_2_i 1))))
  chs))

(define (build-hypercube-comgrp cnt places)
  (define chs (build-hypercube-channels cnt))
  (for ([i (in-range 1 cnt)]
        [pl places])
    (place-channel-send pl (list i cnt (vector-ref chs i))))
  (make-comgrp 0 cnt (vector-ref chs 0)))

(define (receive-hypercube-comgrp ch)
  (apply make-comgrp (place-channel-recv ch)))

(define (comgrp-append-id comgrp v)
  (string-append v (number->string (comgrp-id comgrp))))
(define (mkerror comgrp)
  (string-append "OOPS_ERROR" (number->string (comgrp-id comgrp))))

;(eprintf "~a: ~a~n" id (reverse (comgrp-peers comgrp)))

(define broadcast
  (case-lambda 
    [(comgrp val)
      (let ([id (comgrp-id comgrp)])
        (for/fold ([v val]) ([p (reverse (comgrp-peers comgrp))])
          (match p 
            [(list peerid pch) 
              (if (peerid . > . id)
                (begin (place-channel-send pch v) v)
                (place-channel-recv pch))])))]
    [(comgrp) (broadcast comgrp (mkerror comgrp))]))

(define scatter
  (case-lambda 
    [(comgrp val)
      (let ([id (comgrp-id comgrp)]
            [cnt (comgrp-cnt comgrp)])
        (for/fold ([v val]) ([p (reverse (comgrp-peers comgrp))])
          (match p 
            [(list peerid pch) 
              (if (peerid . > . id)
                (begin 
                  (let* ([left (- cnt peerid)]
                         [possible (subtree-size peerid)]
                         [split-point (/ (min left possible) (+ (min left possible) possible))])
                    (let-values ([(v1 v2) (vector-split-at v (ceiling (* split-point (vector-length v))))]) 
                      (place-channel-send pch v2) v1)))
                (place-channel-recv pch))])))]
    [(comgrp) (scatter comgrp (void))]))

(define (reduce/vector comgrp op right-identity val)
  (let ([id (comgrp-id comgrp)])
    (for/fold ([v val]) ([p (comgrp-peers comgrp)])
      (match p 
        [(list peerid pch) 
          (if (peerid . < . id)
            (begin (place-channel-send pch v) v)
            (vector-map op (place-channel-recv pch) v))]))))

(define-syntax-rule (values->list body ...)
  (call-with-values (lambda () body ...) list))
  
