#lang racket 

(provide main)
(require racket/future)
(require racket/vector)
(require "../bm-args.ss")
(require "../bm-results.ss")
(require "../rand-generator.ss")


;; safe primitives

(require (rename-in scheme [vector-ref vr] [vector-set! vs!])
         scheme/fixnum)

;; unsafe ones
#;
(require (rename-in scheme/unsafe/ops
                    [unsafe-vector-ref vr] [unsafe-vector-set! vs!]
                    [unsafe-fx+ fx+]
                    [unsafe-fx- fx-]
                    [unsafe-fx= fx=]
                    [unsafe-fx<= fx<=]))

(define max-iterations 10)
(define test-array-size 5)

(define total-keys 0)
(define max-key 0)
(define num-buckets 0)
(define num-keys 0)

(define test-index-array #())
(define test-rank-array #())

(define passed-verification 0)

;These are the three main arrays
(define master-hist #())
(define key-array #())
(define partial-verify-vals #())

(define num-threads 0)

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
    [(#\C) (values C-test-index-array C-test-rank-array 27 23 10)]))

;; init-is: character boolean -> void
(define (init-is cls np ser)
  (set! serial ser)
  (set! num-threads np)
  (define-values (tia tra total-keys-log-2 max-key-log-2 num-buckets-log-2)
    (initial-conditions cls))
  ;Common variables 
  (set! test-index-array tia)
  (set! test-rank-array tra)
  (set! total-keys (arithmetic-shift 1 total-keys-log-2))
  (set! max-key (arithmetic-shift 1 max-key-log-2))
  (set! num-buckets (arithmetic-shift 1 num-buckets-log-2))
  (set! num-keys total-keys)
  (set! key-array (make-vector num-keys 0))
  (set! master-hist (make-vector max-key 0))
  (set! partial-verify-vals (make-vector test-array-size 0)))

(define-struct RankThread (id 
                           s1 
                           e1 
                           s2 
                           e2 
                           local-hist 
                           start 
                           end 
                           rstart 
                           rend) 
  #:mutable 
  #:transparent)

(define (new-RankThread id s1 e1 s2 e2) 
  (make-RankThread id 
                   s1 
                   e1 
                   s2 
                   e2 
                   (make-vector max-key 0) 
                   s1 
                   e1 
                   s2 
                   e2))     


(define rankthreads #())
(define rankthreads-local-hist #())
(define local-hist (make-vector 0))

(define (setup-threads) 
  (set! rankthreads (make-vector num-threads 0))
  (let ([chunk-size (quotient total-keys num-threads)]
        [rsize (quotient max-key num-threads)])

      (for/fold ([remainder (remainder total-keys num-threads)]
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

        (vs! rankthreads i (new-RankThread i start end rstart rend))
        (values nr no nrr nro)))
    (set! rankthreads-local-hist (vector-map RankThread-local-hist rankthreads)))

;IS code   
(define serial #f)
(define bid -1)
(define amult 1220703125.0)


; The benchmark entry point.
;;main : list-of-string -> void
(define (main argv) 
  (let ([args (parse-cmd-line-args argv "Integer Sort")])
    (init-is (BMArgs-class args) 
             (BMArgs-num-threads args) 
             (BMArgs-serial args)) 
    (run-benchmark (BMArgs-class args) args)))

(define (run-benchmark class args) 
  (print-banner "Integer Sort" args) 
  (printf "Size: ~a Iterations: ~a~n" total-keys max-iterations) 
  ;Generate random number sequence and subsequent keys on all procs 
  (init-keys amult) ;Random number gen seed 
  ;Do one iteration "free" (i.e. untimed) to guarantee 
  ;initialization of all data and code pages and respective tables 
  ;(print (vector-take key-array 1000))
  (if serial 
      (rank class 1) 
      (begin 
        (setup-threads) 
        (do-sort 1) 
        (partial-verify class 1))) 
  ;Start verification counter 
  (set! passed-verification 0) 
  (if (not (eq? #\S class))
      (printf "~n    iteration#~n")
      void) 
  ;Start timing 
  ;This is the main iteration 
  (let-values ([(r cpu-ms real-ms gc-ms) 
                (time-apply 
                 (λ (x)
                   (for ([it (in-range 1 (+ max-iterations 1))]) 
                     (printf "    Iteration ~a\n" it)
                     (if serial 
                         (rank class it) 
                         (begin 
                           (do-sort it)
                           (partial-verify class it))))) 
                 '(1))])
    ;Stop timing 
    ;This tests that keys are in sequence: sorting of last 
    ;ranked key seq occurs here, but is an untimed operation 
    (full-verify)
    (let ([verified 0])
      (when (= passed-verification (+ (* 5 max-iterations) 1)) 
        (set! verified 1))
      (print-verification-status class verified "Integer Sort")
      (let* ([tm-sec (/ real-ms 1000)]
             [results (make-BMResults "Integer Sort" 
                                      "Machine Name?" 
                                      "PLT Scheme" 
                                      class
                                      total-keys 
                                      0 
                                      0 
                                      max-iterations 
                                      (exact->inexact tm-sec)
                                      0 
                                      0
                                      (get-mops tm-sec 
                                                max-iterations 
                                                total-keys) 
                                      0 
                                      0 
                                      0 
                                      "keys ranked" 
                                      (BMArgs-num-threads args) 
                                      (BMArgs-serial args) 
                                      0 
                                      verified)])                             
        (print-results results)))))   


(define (get-mops total-time niter num-keys)  
  (if (> total-time 0) 
      (* (+ niter num-keys) (/ niter (* total-time 1000000.0)))
      0.0))

(define (rank class iteration) 
  (vs! key-array iteration iteration) 
  (vs! key-array (fx+ iteration max-iterations) (fx- max-key iteration)) 
  (let loop ([i 0])
    (unless (fx= i test-array-size)
      (vs! partial-verify-vals 
           i 
           (vr key-array (vr test-index-array i)))
      (loop (fx+ i 1))))
  ;Clear the work array 
  (let loop ([i 0])
    (unless (fx= i max-key)
      (vs! master-hist i 0)
      (loop (fx+ i 1))))
  ;In this section, the keys themselves are used as their 
  ;own indices to determine how many of each there are: their 
  ;individual population 
  (let ([master-hist master-hist]
        [key-array key-array]
        [num-keys num-keys])
    (let loop ([i 0])
      (unless (fx= i num-keys)
        (let ([ki (vr key-array i)])
          (vs! master-hist 
               ki
               (fx+ (vr master-hist ki)  
                    1))
          (loop (fx+ i 1))))))
  ;Now they have individual key population 
  ;Do density to distribution conversion 
  (let ([end (- max-key 1)])
    (let loop ([i 0])
      (unless (fx= i end)
        (vs! master-hist 
             (fx+ i 1) 
             (fx+ (vr master-hist (fx+ i 1)) 
                  (vr master-hist i)))
        (loop (fx+ i 1)))))
  (partial-verify class iteration))

;; partial-verify : int -> void
(define (partial-verify class iteration) 
  (for ([i (in-range 0 test-array-size)]) 
    (let ([k (vr partial-verify-vals i)] 
          [offset iteration])
      (if (and (<= 0 k) (<= k (- num-keys 1))) 
          (begin 
            (case class
              [(#\S) 
               (if (<= i 2) 
                   (set! offset iteration) 
                   (set! offset (- 0 iteration)))]
              [(#\W) 
               (if (< i 2) 
                   (set! offset (- iteration 2)) 
                   (set! offset (- 0 iteration)))]
              [(#\A) 
               (if (<= i 2) 
                   (set! offset (- iteration 1)) 
                   (set! offset (+ (- 0 iteration) 1)))]
              [(#\B) 
               (if (or (= i 1) (= i 2) (= i 4)) 
                   (set! offset iteration) 
                   (set! offset (- 0 iteration)))]
              [(#\C) 
               (if (<= i 2) 
                   (set! offset iteration) 
                   (set! offset (- 0 iteration)))])
            (if (not (= (vr master-hist (- k 1)) 
                        (+ (vr test-rank-array i) offset)))
                (begin 
                  (printf "Failed partial verification: iteration~a, test key ~a~n" 
                          iteration 
                          i) 
                  #;(printf "Expected ~a at master-hist[~a], but got ~a instead~n" 
                            (+ (vr test-rank-array i) offset) 
                            (- k 1) 
                            (vr master-hist (- k 1)))
                  #;(if (> (vector-count (λ (x) (= x (- k 1))) key-array) 0) 
                        (printf "~a is in the original key array.~n" (- k 1)) 
                        (printf "~a is NOT in the original key array.~n" (- k 1)))
                  #;(printf "Surround: ~a~n" 
                            (for/list ([n (in-range (- k 6) (+ k 5))]) 
                              (vr master-hist n))))
                (set! passed-verification (+ passed-verification 1))))
          void))))

(define (full-verify) 
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
      (set! idx (+ idx 1)))
    
    (let ([count 0])
      (for ([i (in-range 1 num-keys)])
        (when (> (vr key-array (- i 1))
                 (vr key-array i))
          (set! count (+ count 1))))
      (if (= count 0)
          (set! passed-verification (+ passed-verification 1))
          (printf "Full_verify: number of keys out of sort: ~s\n" count)))))

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



(define (init-keys a) 
  (define k (/ max-key 4)) 
  (for ([i (in-range 0 num-keys)])
    (define x (+ (randlc a) (randlc a) (randlc a) (randlc a)))
    (vs! key-array i (inexact->exact (floor (* x k))))))

;; step1 : void
(define (step1 rt)
  (let ([local-hist (RankThread-local-hist rt)]
        [start (RankThread-start rt)] 
        [end (RankThread-end rt)]
        [key-array key-array])   ;; make this reference faster
    (let loop ([i 0])
      (unless (fx= i max-key)
        (vs! local-hist i 0)
        (loop (fx+ i 1)))) 
    (let loop ([i start])
      (when (fx<= i end)
        (let ([ki (vr key-array i)])
          (vs! local-hist 
               ki
               (fx+ (vr local-hist ki) 1))
          (loop (fx+ i 1)))))
    (let ([end (- max-key 1)])
      (let loop ([i 0])
        (unless (fx= i end)
          (vs! local-hist 
               (fx+ i 1) 
               (fx+ (vr local-hist (fx+ i 1)) 
                    (vr local-hist i)))
          (loop (fx+ i 1)))))))


;Parallel calculation of the master's histogram
(define (step2 rt) 
  (let ([local-hist (RankThread-local-hist rt)] 
        [rstart (RankThread-rstart rt)] 
        [rend (RankThread-rend rt)])
    
    (let loop ([i rstart])
      (when (fx<= i rend)
        (let loop ([j 0])
          (unless (fx= j num-threads)
            (vs! master-hist 
                 i 
                 (fx+
                  (vr master-hist i) 
                  (vr (vr rankthreads-local-hist j) i)))
            (loop (fx+ j 1))))
        (loop (fx+ i 1))))))

(define (do-sort iteration) 
  (vs! key-array iteration iteration)
  (vs! key-array (fx+ iteration max-iterations) (fx- max-key iteration))
  (for ([i (in-range 0 test-array-size)])
    (vs! partial-verify-vals i (vr key-array (vr test-index-array i))))
  (let ([step1-futures (vector-map (λ (rthread) 
                                     (future (λ () 
                                               (step1 rthread))))
                                   rankthreads)])
    (for ((i (in-range 0 (vector-length step1-futures))))
      (touch (vr step1-futures i))))
  (let loop ([i 0])
    (unless (fx= i max-key)
      (vs! master-hist i 0)
      (loop (fx+ i 1))))
  (let ([step2-futures (vector-map (λ (rthread) 
                                     (future (λ () 
                                               (step2 rthread))))
                                   rankthreads)]) 
    (for ((i (in-range 0 (vector-length step2-futures))))
      (touch (vr step2-futures i)))))

(main (current-command-line-arguments))
