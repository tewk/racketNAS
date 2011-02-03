#lang racket/base

(require racket/future)
(require racket/vector)
(require "bm-args.rkt")
(require "bm-results.rkt")
(require "rand-generator.rkt")
(require "timer.rkt")
(require "parallel-utils.rkt")

(provide main)

;; safe primitives
;;#;
(require (rename-in racket [vector-ref vr] [vector-set! vs!])
         racket/fixnum)

;; unsafe ones
#;
(require (rename-in racket/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-fxvector-ref fxr] 
                    [unsafe-fxvector-set! fx!]
                    [unsafe-fx+ fx+]
                    [unsafe-fx- fx-]
                    [unsafe-fx= fx=]
                    [unsafe-fx<= fx<=]))

(define fxr fxvector-ref)
(define fx! fxvector-set!)
(define (initial-conditions cls)
  (define S-test-index-array #(48427 17148 23627 62548 4431))
  (define S-test-rank-array #(0 18 346 64917 65463))
  (define W-test-index-array #(357773 934767 875723 898999 404505))
  (define W-test-rank-array #(1249 11698 1039987 1043896 1048018))
  (define A-test-index-array #(2112377 662041 5336171 3642833 4250760))
  (define A-test-rank-array #(104 17523 123928 8288932 8388264))
  (define B-test-index-array #(41869 812306 5102857 18232239 26860214))
  (define B-test-rank-array #(33422937 10244 59149 33135281 99))
  (define C-test-index-array #(44172927 72999161 74326391 129606274 21736814))
  (define C-test-rank-array #(61147 882988 266290 133997595 133525895))
  (case cls 
    [(#\S) (values S-test-index-array S-test-rank-array 16 11 9)]       
    [(#\W) (values W-test-index-array W-test-rank-array 20 16 10)]
    [(#\A) (values A-test-index-array A-test-rank-array 23 19 10)]
    [(#\B) (values B-test-index-array B-test-rank-array 25 21 10)]
    [(#\C) (values C-test-index-array C-test-rank-array 27 23 10)]
    [else (error (format "Unknown class ~a" cls))]))

(define (main . argv) 
  (let* ([args (parse-cmd-line-args argv "Integer Sort")]
         [class (BMArgs-class args)]
         [serial (BMArgs-serial args)]
         [num-threads (BMArgs-num-threads args)]
         [bmname "IS"]) 
    (define-values (test-index-array test-rank-array total-keys-log-2 max-key-log-2 num-buckets-log-2)
      (initial-conditions class))
    (let* ([num-keys (arithmetic-shift 1 total-keys-log-2)]
           [max-key (arithmetic-shift 1 max-key-log-2)]
           [num-buckets (arithmetic-shift 1 num-buckets-log-2)]
           [test-array-size 5]
           [max-iterations 10]
           [key-array (make-shared-fxvector num-keys 0)]
           [master-hist (make-shared-fxvector max-key 0)]
           [partial-verify-vals (make-shared-fxvector test-array-size 0)])

  (define amult 1220703125.0)
  (define (init-keys a) 
    (define k (/ max-key 4)) 
    (for ([i (in-range 0 num-keys)])
      (define x (+ (randlc a) (randlc a) (randlc a) (randlc a)))
      (fx! key-array i (inexact->exact (floor (* x k))))))




  (print-banner "Integer Sort" args) 
  (printf "Size: ~a Iterations: ~a~n" num-keys max-iterations) 

  ;Generate random number sequence and subsequent keys on all procs 
  (init-keys amult) ;Random number gen seed 

  (when (not (eq? #\S class)) (printf "~n    Iteration #~n"))

  (define local-hists (for/vector ([i (in-range num-threads)]) (make-shared-fxvector max-key 0)))

  (if serial
  (let ()
  (define-syntax-rule (sitit it)
    (begin
      (define (prep-iteration iteration)
        (fx! key-array iteration iteration) 
        (fx! key-array (fx+ iteration max-iterations) (fx- max-key iteration)) 
        (for ([i (in-range test-array-size)])
          (fx! partial-verify-vals i (fxr key-array (vr test-index-array i)))))
      (define (serial-sort)
          (for ([i (in-range max-key)])
            (fx! master-hist i 0))

          (for ([i (in-range num-keys)])
            (fx++! master-hist (fxr key-array i)))

          (for/fold ([init (fxr master-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
            (define sum (fx+ init (fxr master-hist i)))
            (fx! master-hist i sum)
            sum))
      (printf "    Iteration ~a\n" it)
      (prep-iteration it)
      (serial-sort)
      (partial-verify class it master-hist partial-verify-vals test-rank-array num-keys)))

  (sitit 1)
  (timer-start 1)
  (for ([it (in-range 1 (+ max-iterations 1))]) (sitit it))
  (timer-stop 1))

  (CGspawn (if serial 0 num-threads) is-body key-array test-rank-array test-array-size 
    test-index-array partial-verify-vals master-hist num-keys max-key max-iterations class num-threads local-hists))

  ;This tests that keys are in sequence: sorting of last 
  ;ranked key seq occurs here, but is an untimed operation 
  (let ([verified (full-verify master-hist key-array max-key num-keys)])
    (print-verification-status class verified bmname)
    (let* ([tm-sec (/ (read-timer 1) 1000)]
           [results (new-BMResults "IS" 
                                    class
                                    num-keys 
                                    0 
                                    0 
                                    max-iterations 
                                    (exact->inexact tm-sec)
                                    (get-mops tm-sec max-iterations num-keys) 
                                    "keys ranked" 
                                    verified
                                    serial
                                    num-threads
                                    0)])                             
      (print-results results))))))

(define (is-body cg key-array test-rank-array test-array-size test-index-array partial-verify-vals master-hist num-keys max-key max-iterations class num-threads local-hists)
  (define (prep-iteration iteration)
    (fx! key-array iteration iteration) 
    (fx! key-array (fx+ iteration max-iterations) (fx- max-key iteration)) 
    (for ([i (in-range test-array-size)])
      (fx! partial-verify-vals i (fxr key-array (vr test-index-array i)))))
  (define local-hist (if (= (CGnp cg) 0) master-hist (vr local-hists (CGid cg))))
  (define-syntax-rule (itit it)
    (begin
    (CG-n0-only cg
      (printf "    Iteration ~a\n" it)
      (prep-iteration it))

    (for ([i (in-range max-key)])
        (fx! local-hist i 0))

    ;In this section, the keys themselves are used as their 
    ;own indices to determine how many of each there are: their 
    ;individual population 
    (CGfor cg ([i (in-range num-keys)])
        (fx++! local-hist (fxr key-array i)))

    ;Now they have individual key population 
    ;Do density to distribution conversion 
    (for/fold ([init (fxr local-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
      (define sum (fx+ init (fxr local-hist i)))
      (fx! local-hist i sum)
      sum)

    (CG-B cg)

    (CG-Parallel-Only cg 
      ;Parallel calculation of the master's histogram
      (CGfor cg ([i (in-range max-key)])
        (fx! master-hist i 
             (for/fold ([accum 0]) ([j (in-range num-threads)])
               (fx+ accum (fxr (vr local-hists j) i))))))
    (CG-n0-only cg
      (partial-verify class it master-hist partial-verify-vals test-rank-array num-keys))))

  ;Do one iteration "free" (i.e. untimed) to guarantee 
  ;initialization of all data and code pages and respective tables 
  (itit 1)
  ;(when (not (eq? #\S class)) (printf "~n    iteration#~n"))

  (timer-start 1)

  ;This is the main iteration 
  (for ([it (in-range 1 (+ max-iterations 1))]) 
    (itit it))
  (timer-stop 1))


(define (get-mops total-time niter num-keys)  
  (if (> total-time 0) 
      (* (+ niter num-keys) (/ niter (* total-time 1000000.0)))
      0.0))

;; partial-verify : int -> void
(define (partial-verify class iteration master-hist partial-verify-vals test-rank-array num-keys) 
  (for ([i (in-range (fxvector-length partial-verify-vals))]) 
    (let ([k (fxr partial-verify-vals i)])
      (if (and (<= 0 k) (<= k (- num-keys 1))) 
          (let ([offset
            (case class
              [(#\S) 
               (if (<= i 2) 
                   iteration
                   (- 0 iteration))]
              [(#\W) 
               (if (< i 2) 
                   (- iteration 2)
                   (- 0 iteration))]
              [(#\A) 
               (if (<= i 2) 
                   (- iteration 1)
                   (+ (- 0 iteration) 1))]
              [(#\B) 
               (if (or (= i 1) (= i 2) (= i 4)) 
                   iteration
                   (- 0 iteration))]
              [(#\C) 
               (if (<= i 2) 
                   iteration
                   (- 0 iteration))]
              [else (error "Unknown class")])])
            (let ([mh (fxr master-hist (- k 1))]
                  [th (+ (vr test-rank-array i) offset)])
              (when (not (= mh th))
                (error (format "Failed partial verification: iteration ~a, test key ~a ~a ~a~n" iteration i mh th)))))
          #t))))

(define (full-verify master-hist key-array max-key num-keys) 
  (let ([key 0] 
        [idx 0]) 
    (for ([i (in-range 0 num-keys)]) 
      (let/ec k
        (let loop ()
          (when (= idx (fxr master-hist key))
            (set! key (+ key 1))
            (when (or (>= key max-key)
                      (>= idx num-keys))
              (k (void))))))
      (fx! key-array idx key)
      (set! idx (+ idx 1))))
    
    (let-values ([(prev count)
        (for/fold ([prev (fxr key-array 0)]
                 [count 0])
                 ([i (in-range 1 num-keys)])
          (let ([curr (fxr key-array i)])
            (values curr 
                    (if (prev . > . curr)
                        (fx+ count 1)
                        count))))])
      (if (= count 0)
          #t
          (begin 
            (printf "Full_verify: number of keys out of sort: ~s\n" count)
            #f))))

;; full-verify : -> void 
#;(define (full-verify) 
    ;To save copy and memory sorting can be done directly 
    (let ([key 0]
          [idx 0]) 
      (for ([i (in-range 0 num-keys)]) 
        (for/and ([j (in-naturals)]) 
          (if (= idx (fxr master-hist key)) 
              (begin 
                (set! key (+ key 1)) 
                (if (or (>= key max-key) (>= idx num-keys)) 
                    #f 
                    #t))
              #f))
        (vs! key-array idx key) 
        (set! idx (+ idx 1)))))


(define-syntax-rule (v++! v idx)
  (vs! v idx (fx+ (vr v idx) 1)))

(define-syntax-rule (fx++! v idx)
  (fx! v idx (fx+ (fxr v idx) 1)))


(define-syntax-rule (future-for func v)
  (let ([fs (for/list ([x (in-range v)])
              (future (λ () (func x))))])
    (for ([x fs])
      (touch x))))
  
;; PLACES MPI SORT
#|
    [#f
      (let ()
        (define places (for/list ([i (in-range 1 num-threads)]) (place "integer-sort.rkt" 'sort-start-worker)))
        (define comgrp (build-hypercube-comgrp num-threads places))
        (broadcast comgrp (list max-key (add1 max-iterations)))

        (lambda (key-array num-keys master-hist max-key)
          (set! master-hist (places-sort comgrp key-array max-key))
          master-hist))]


(define (places-sort comgrp key-array max-key)
  (sort-worker comgrp (scatter comgrp key-array) max-key))

(define (sort-start-worker ch)
  (define comgrp (receive-hypercube-comgrp ch))
  (match (broadcast comgrp)
    [(list max-key its)
      (for ([i (in-range its)])
        (sort-worker comgrp (scatter comgrp) max-key))]))

(define (sort-worker comgrp key-array max-key)
  (define local-hist (make-vector max-key 0))
  
  (for ([i (in-range (vector-length key-array))])
    (v++! local-hist (vr key-array i)))

  (for/fold ([init (vr local-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
    (define sum (fx+ init (vr local-hist i)))
    (vs! local-hist i sum)
    sum)

  (reduce/vector comgrp + 0 local-hist))
|#
;;;        (let ([local-hists (for/vector ([i (in-range num-threads)]) (make-vector max-key 0))])
;;;          (lambda (key-array num-keys master-hist max-key)
;;;            (future-for (lambda (x) (step1 x key-array master-hist max-key (vr local-hists x) num-keys num-threads)) num-threads)
;;;            (future-for (lambda (x) (step2 x master-hist max-key local-hists num-keys num-threads)) num-threads)
;;;            master-hist)))]

