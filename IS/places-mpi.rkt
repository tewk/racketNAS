#lang racket 
;
(define-struct comgrp (id cnt peers))

(define (log-2-ceiling n)
  (let loop ([p 0]
             [m 1])
    (if (m . > . n)
        p
        (loop (add1 p) (arithmetic-shift m 1)))))

(define-syntax-rule (vsa! v idx val)
  (vector-set! v idx (append (vector-ref v idx) (list val))))

(define (build-hypercube-channels n)
  (define chs (make-vector n '()))
  (let ([p2n (p2 (- n 1))])
    (for ([i (in-range p2n)])
      (for ([j (in-range n)])
        (define peer-id (bitwise-xor (expt 2 i) j))
        (when (and (zero? (bitwise-and peer-id (- (expt 2 i) 1)))
                   (not (zero? (bitwise-and peer-id (expt 2 i)))))
            (let-values ([(ch1 ch2) (place-channel)])
          (vsa! chs j (list peer-id ch1)
          (vsa! chs peer-id (list j ch2))))))))
  chs)

(define (build-hypercube-comgrp cnt places)
  (define chs (build-hypercube-ch cnt))
  (for ([i (in-range 1 cnt)]
        [pl places])
    (place-channel-send pl (list i (vector-ref chs i))))
  (make-comgrp 0 (vector-ref chs 0)))

(define broadcast
  (case-lambda 
    [(comgrp v)
      (let ([id (comgrp-id comgrp)])
        (for/fold ([val v]) ([p (reverse (comgrp-peers comgrp))])
          (match p 
            [(list peerid pch) 
              (if (peerid . > . id)
                (begin (place-channel-send pch v) v)
                (place-channel-recv pch))])))]
    [(comgrp) (broadcast comgrp (void))]))

(define scatter
  (case-lambda 
    [(comgrp v)
      (let ([id (comgrp-id comgrp)]
            [cnt (comgrp-vnt comgrp)])
        (for/fold ([val v]) ([p (reverse (comgrp-peers comgrp))])
          (match p 
            [(list peerid pch) 
              (if (peerid . > . id)
                (begin 
                  (let* ([left (- cnt peerid)]
                         [possible (p2n peerid)]
                         [split-point (/ (min left possible) (+ (min left possible) possible))])
                    (let-values ([(v1 v2) (vector-split-at v (ceiling (* split-point (vector-length v))))]) 
                      (place-channel-send pch v2) v1)))
                (place-channel-recv pch))])))]
    [(comgrp) (split comgrp (void))]))

(define (reduce/vector comgrp op right-identity v)
  (let ([id (comgrp-id comgrp)])
    (for/fold ([val v]) ([p (comgrp-peers comgrp)])
      (match p 
        [(list peerid pch) 
          (if (peerid . < . id)
            (begin (place-channel-send pch v) v)
            (vector-zip + (place-channel-recv pch) v))]))))

(define-syntax-rule (values->list body ...)
  (call-with-values (lambda () body ...) list))
  
(define places (for/list ([i (in-range 1 wrk-cnt)]) (place "integer-sort.ss" 'sort-start-worker)))

(define (places-sort places wrk-cnt key-array max-key)
  (define comgrp (build-hypercube-comgrp wrk-cnt places))
  (sort-worker comgrp (scatter comgrp key-array) (broadcast comgrp max-key)))

(define (sort-start-worker ch)
  (define comgrp (receive-hypercube-comgrp ch))
  (sort-worker comgrp (scatter comgrp) (broadcast comgrp)))

(define (sort-worker comgrp key-array max-key)
    (define local-hist (make-vector max-key 0))
    
    (for ([i (in-range (vector-length key-array))])
      (v++! local-hist (vr key-array i)))

    (for/fold ([init (vr local-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
      (define sum (fx+ init (vr local-hist i)))
      (vs! local-hist i sum)
      sum)

    (reduce/vector comgrp + 0 local-hist))
