#lang racket 

;; racket -tm integer-sort.rkt CLASS=S SERIAL
;; racket -tm integer-sort.rkt CLASS=S NP=2
;; racket -tm integer-sort.rkt CLASS=S NP=2

(require racket/future)
(require racket/vector)
(require "../bm-args.rkt")
(require "../bm-results.rkt")
(require "../rand-generator.rkt")
(require "../timer.rkt")
(require "places-mpi.rkt")

(provide main
         sort-start-worker)

;; safe primitives
;;#;
(require (rename-in racket [vector-ref vr] [vector-set! vs!])
         racket/fixnum)

;; unsafe ones
#;
(require (rename-in racket/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-fx+ fx+]
                    [unsafe-fx- fx-]
                    [unsafe-fx= fx=]
                    [unsafe-fx<= fx<=]))

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

(define-struct IS-Params (serial num-threads max-key num-buckets num-keys))

(define-struct RankThread (id 
                           local-hist 
                           start 
                           end 
                           rstart 
                           rend) 
  #:mutable 
  #:transparent)


(define rankthreads #())
(define rankthreads-local-hist #())
(define local-hist (make-vector 0))

(define (setup-threads num-threads max-key num-keys) 
  (define (new-RankThread id s1 e1 s2 e2) 
    (make-RankThread id 
                     (make-vector max-key 0) 
                     s1 
                     e1 
                     s2 
                     e2))     

  (set! rankthreads (make-vector num-threads 0))
  (let ([chunk-size (quotient num-keys num-threads)]
        [rsize (quotient max-key num-threads)])

      (for/fold ([remainder (remainder num-keys num-threads)]
                 [offset 0] 
                 [rremainder (remainder max-key num-threads)]
                 [roffset 0]) 
                ([i (in-range 0 num-threads)])
        (define (segmenter chunk-size offset remainder)
          (define start (+ (* i chunk-size) offset)) 
          (define end (+ (- (+ (* i chunk-size) chunk-size) 1) offset))
          (if (> remainder 0)
            (values start (- remainder 1) (+ offset 1) (+ end 1))
            (values start remainder offset end)))
        (define-values (start nr no end) (segmenter chunk-size offset remainder))
        (define-values (rstart nrr nro rend) (segmenter rsize roffset rremainder))
        (define-values (s1 e2) (strip-indicies i num-keys num-threads))

        (vs! rankthreads i (new-RankThread i start end rstart rend))
        (values nr no nrr nro)))
    (set! rankthreads-local-hist (vector-map RankThread-local-hist rankthreads)))

(define (strip-indicies i s1 s2)
  ;(define-values (chunk-size remainder) (quotient/remainder s1 s2))
  (define chunk-size (quotient s1 s2))
  (define rem (remainder s1 s2))
  (define offset (if (i . < . rem) i rem))
  (define start (+ (* i chunk-size) offset))
  (define end (+ (- (* (+ i 1) chunk-size) 1) offset))
  (values start end))


;IS code   
(define bid -1)

; The benchmark entry point.
;;main : list-of-string -> void
(define (main . argv) 
  (let* ([args (parse-cmd-line-args argv "Integer Sort")]
         [class (BMArgs-class args)]
         [serial (BMArgs-serial args)]
         [num-threads (BMArgs-num-threads args)]) 
    (define-values (test-index-array test-rank-array total-keys-log-2 max-key-log-2 num-buckets-log-2)
      (initial-conditions class))
    (let* ([num-keys (arithmetic-shift 1 total-keys-log-2)]
           [max-key (arithmetic-shift 1 max-key-log-2)]
           [num-buckets (arithmetic-shift 1 num-buckets-log-2)]
           [test-array-size 5]
           [max-iterations 10]
           [key-array (make-vector num-keys 0)]
           [master-hist (make-vector max-key 0)]
           [partial-verify-vals (make-vector test-array-size 0)])

  (define amult 1220703125.0)
  (define (init-keys a) 
    (define k (/ max-key 4)) 
    (for ([i (in-range 0 num-keys)])
      (define x (+ (randlc a) (randlc a) (randlc a) (randlc a)))
      (vs! key-array i (inexact->exact (floor (* x k))))))


  (define (prep-iteration iteration)
    (vs! key-array iteration iteration) 
    (vs! key-array (fx+ iteration max-iterations) (fx- max-key iteration)) 
    (for ([i (in-range test-array-size)])
      (vs! partial-verify-vals i (vr key-array (vr test-index-array i)))))


  (print-banner "Integer Sort" args) 
  (printf "Size: ~a Iterations: ~a~n" num-keys max-iterations) 
  ;Generate random number sequence and subsequent keys on all procs 


  (init-keys amult) ;Random number gen seed 
  ;Do one iteration "free" (i.e. untimed) to guarantee 
  ;initialization of all data and code pages and respective tables 
  ;(print (vector-take key-array 1000))

  (define sorter (match serial
    [#t (lambda (key-array num-keys master-hist max-key)
      (serial-sort key-array num-keys master-hist max-key)
      master-hist)]
    [#f
      (let ()
        (define places (for/list ([i (in-range 1 num-threads)]) (place "integer-sort.rkt" 'sort-start-worker)))
        (define comgrp (build-hypercube-comgrp num-threads places))
        (broadcast comgrp (list max-key (add1 max-iterations)))

        (lambda (key-array num-keys master-hist max-key)
          (set! master-hist (places-sort comgrp key-array max-key))
          master-hist))]
    ["FUTURES"
      (let ()
        (setup-threads num-threads max-key)
        (let ([local-hists (for/vector ([i (in-range num-threads)]) (make-vector max-key 0))])
          (lambda (key-array num-keys master-hist max-key)
            (future-for (lambda (x) (step1 x key-array master-hist max-key (vr local-hists x) num-keys num-threads)) num-threads)
            (future-for (lambda (x) (step2 x master-hist max-key local-hists num-keys num-threads)) num-threads)
            master-hist)))]
    [else
      (raise "OOPS")]))

  (prep-iteration 1)
  (sorter key-array num-keys master-hist max-key)
  (partial-verify class 1 master-hist partial-verify-vals test-rank-array num-keys) 

  (when (not (eq? #\S class)) (printf "~n    iteration#~n"))

  (timer-start 1)

  ;This is the main iteration 
  (let ([master-hist (for/fold ([master-hist master-hist]) ([it (in-range 1 (+ max-iterations 1))]) 
      (printf "    Iteration ~a\n" it)
      (prep-iteration it)
      (let ([master-hist (sorter key-array num-keys master-hist max-key)])
         (partial-verify class it master-hist partial-verify-vals test-rank-array num-keys)
         master-hist))])

    (timer-stop 1)

    ;This tests that keys are in sequence: sorting of last 
    ;ranked key seq occurs here, but is an untimed operation 
    (let ([verified (full-verify master-hist key-array max-key num-keys)])
      (print-verification-status class verified "Integer Sort")
      (let* ([tm-sec (/ (read-timer 1) 1000)]
             [results (make-BMResults "IS" 
                                      "Machine Name?" 
                                      "PLT Racket" 
                                      class
                                      num-keys 
                                      0 
                                      0 
                                      max-iterations 
                                      (exact->inexact tm-sec)
                                      0 
                                      0
                                      (get-mops tm-sec max-iterations num-keys) 
                                      0 
                                      0 
                                      0 
                                      "keys ranked" 
                                      num-threads
                                      serial
                                      0 
                                      verified)])                             
        (print-results results)))))))


(define (get-mops total-time niter num-keys)  
  (if (> total-time 0) 
      (* (+ niter num-keys) (/ niter (* total-time 1000000.0)))
      0.0))

;NEEDS max-key master-hist key-array num-keys
(define (serial-sort key-array num-keys master-hist max-key) 
    ;Clear the work array 
    (for ([i (in-range max-key)])
        (vs! master-hist i 0))
    ;In this section, the keys themselves are used as their 
    ;own indices to determine how many of each there are: their 
    ;individual population 
    (for ([i (in-range num-keys)])
        (v++! master-hist (vr key-array i)))
            
    ;Now they have individual key population 
    ;Do density to distribution conversion 
    (for/fold ([init (vr master-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
      (define sum (fx+ init (vr master-hist i)))
      (vs! master-hist i sum)
      sum))

;; partial-verify : int -> void
;; deps num-keys 
(define (partial-verify class iteration master-hist partial-verify-vals test-rank-array num-keys) 
  (for ([i (in-range (vector-length partial-verify-vals))]) 
    (let ([k (vr partial-verify-vals i)])
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
            (let ([mh (vr master-hist (- k 1))]
                  [th (+ (vr test-rank-array i) offset)])
              (when (not (= mh th))
                (error (format(printf "Failed partial verification: iteration ~a, test key ~a ~a ~a~n" iteration i mh th))))))
          #t))))

(define (full-verify master-hist key-array max-key num-keys) 
  (let ([key 0] 
        [idx 0]) 
    (for ([i (in-range 0 num-keys)]) 
      (let/ec k
        (let loop ()
          (when (= idx (vr master-hist key))
            (set! key (+ key 1))
            (when (or (>= key max-key)
                      (>= idx num-keys))
              (k (void))))))
      (vs! key-array idx key)
      (set! idx (+ idx 1))))
    
    (let-values ([(prev count)
        (for/fold ([prev (vr key-array 0)]
                 [count 0])
                 ([i (in-range 1 num-keys)])
          (let ([curr (vr key-array i)])
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
          (if (= idx (vr master-hist key)) 
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

;; step1 : void
;; deps max-key local-hist start end key-array 
(define (step1 id key-array master-hist max-key local-hist num-keys num-threads)
  (define-values (start end) (strip-indicies id num-keys num-threads))
  (for ([i (in-range max-key)])
      (vs! local-hist i 0))
  ;In this section, the keys themselves are used as their 
  ;own indices to determine how many of each there are: their 
  ;individual population 
  (for ([i (in-range start (fx+ end 1))])
      (v++! local-hist (vr key-array i)))
          
  ;Now they have individual key population 
  ;Do density to distribution conversion 
  (for/fold ([init (vr local-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
    (define sum (fx+ init (vr local-hist i)))
    (vs! local-hist i sum)
    sum))

;Parallel calculation of the master's histogram
;; deps num-threads local-hist rstart rend
(define (step2 id master-hist max-key local-hists num-keys num-threads)
  (define-values (rstart rend) (strip-indicies id max-key num-threads))
    (for ([i (in-range rstart (fx+ rend 1))])
      (vs! master-hist 
           i 
           (for/fold ([accum 0]) ([j (in-range num-threads)])
             (fx+ accum (vr (vr local-hists j) i))))))

(define-syntax-rule (future-pmap/vector func v)
  (let ([fs (vector-map (λ (x) (future (λ () (func x)))) v)])
    (for ((i (in-range (vector-length fs))))
      (touch (vr fs i)))))
  
(define-syntax-rule (future-for func v)
  (let ([fs (for/list ([x (in-range v)])
              (future (λ () (func x))))])
    (for ([x fs])
      (touch x))))
  
(define (future-sort2 a key-array num-threads master-hist max-key) 
   ;; deps max-key local-hist start end key-array
  (define (step1 rt)
    (let ([local-hist (RankThread-local-hist rt)]
          [start (RankThread-start rt)] 
          [end (RankThread-end rt)]
          [key-array key-array])   ;; make this reference faster
      (for ([i (in-range max-key)])
          (vs! local-hist i 0))
      ;In this section, the keys themselves are used as their 
      ;own indices to determine how many of each there are: their 
      ;individual population 
      (for ([i (in-range start (fx+ end 1))])
          (v++! local-hist (vr key-array i)))
              
      ;Now they have individual key population 
      ;Do density to distribution conversion 
      (for/fold ([init (vr local-hist 0)]) ([i (in-range 1 (fx- max-key 1))])
        (define sum (fx+ init (vr local-hist i)))
        (vs! local-hist i sum)
        sum)))

   ;Parallel calculation of the master's histogram
   ;; deps num-threads local-hist rstart rend
  (define (step2 rt) 
    (let ([local-hist (RankThread-local-hist rt)] 
          [rstart (RankThread-rstart rt)] 
          [rend (RankThread-rend rt)])
      
       (for ([i (in-range rstart (fx+ rend 1))])
         (vs! master-hist
              i
              (for/fold ([accum 0]) ([j (in-range num-threads)])
               (fx+ accum (vr (vr rankthreads-local-hist j) i)))))))

  (future-pmap/vector step1 rankthreads)
  (future-pmap/vector step2 rankthreads))


;; PLACES MPI SORT

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

;; PLACES SHARED MEMORY SORT

(define (sort-start-worker2 ch)
  (define comgrp (receive-hypercube-comgrp ch))
  (match (broadcast comgrp)
    [(list master-hist max-key key-array num-keys local-hists num-threads its)
      (for ([i (in-range its)])
        (sort-worker2 comgrp master-hist max-key key-array num-keys local-hists num-threads))]))

(define (sort-worker2 comgrp master-hist max-key key-array num-keys local-hists num-threads)
    (let* ([id (comgrp-id comgrp)]
           [local-hist (vr local-hists id)])
;          (barrier compgrp)
        (step1 id key-array master-hist max-key local-hists num-keys num-threads)
;          (barrier compgrp)
        (step2 id master-hist max-key local-hists num-keys num-threads)))
