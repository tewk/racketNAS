#lang racket/base

(require racket/place/distributed
         racket/place/distributed/rmpi
         racket/vector
         racket/match
         racket/flonum
         racket/format
         "bm-args.rkt"
         "bm-results.rkt"
         "rand-generator.rkt"
         "timer.rkt")

;; safe primitives
(define-syntax define/provide
  (syntax-rules ()
    [(_ (name x ...) body ...)
     (begin (provide name)
            (define (name x ...) body ...))]
    [(_ name val)
     (begin (provide name)
            (define name val))]))

(require (only-in racket [vector-ref vr] [vector-set! vs!])
                 
         (rename-in racket/fixnum [fxvector-ref fxr] [fxvector-set! fx!])
         (rename-in racket/flonum [flvector-ref flr] [flvector-set! fl!])
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

(provide main)

(define-values (T_TOTAL T_RANK T_RCOMM T_VERIFY T_LAST) (values 0 1 2 3 4))

(define (initial-conditions cls)
  (define S-test-index-array (fxvector 48427 17148 23627 62548 4431))
  (define S-test-rank-array (fxvector 0 18 346 64917 65463))
  (define W-test-index-array (fxvector 357773 934767 875723 898999 404505))
  (define W-test-rank-array (fxvector 1249 11698 1039987 1043896 1048018))
  (define A-test-index-array (fxvector 2112377 662041 5336171 3642833 4250760))
  (define A-test-rank-array (fxvector 104 17523 123928 8288932 8388264))
  (define B-test-index-array (fxvector 41869 812306 5102857 18232239 26860214))
  (define B-test-rank-array (fxvector 33422937 10244 59149 33135281 99))
  (define C-test-index-array (fxvector 44172927 72999161 74326391 129606274 21736814))
  (define C-test-rank-array (fxvector 61147 882988 266290 133997595 133525895))
  (define D-test-index-array (fxvector 1317351170 995930646 1157283250 1503301535 1453734525))
  (define D-test-rank-array (fxvector 1 36538729 1978098519 2145192618 2147425337))
  (case cls 
    [(#\S) (values S-test-index-array S-test-rank-array 16 11 9 1)]       
    [(#\W) (values W-test-index-array W-test-rank-array 20 16 10 1)]
    [(#\A) (values A-test-index-array A-test-rank-array 23 19 10 1)]
    [(#\B) (values B-test-index-array B-test-rank-array 25 21 10 1)]
    [(#\C) (values C-test-index-array C-test-rank-array 27 23 10 1)]
    [(#\D) (values D-test-index-array D-test-rank-array 29 27 10 4)]
    [else (error (format "Unknown class ~a" cls))]))

(define (find-my-seed id cnt nn seed a)
  (define (log-2-floor x)
    (let loop ([p2 0]
               [x x])
      (if (> x 1)
          (loop (add1 p2) (quotient x 2))
          p2)))

  (define-values (an _)
    (for/fold ([s a]
               [a2  a])
              ([i (in-range 1 (add1 (log-2-floor (quotient nn cnt))))])
      (randlc/2 s s)))

  (let loop ([t1 seed]
             [t2 an]
             [kk id]
             [i 1])
    (cond
      [(= i 101) t1]
      [else
        (define ik (quotient kk 2))
        (define (finish t1 t2)
           (cond
             [(= ik 0) t1]
             [else
               (let-values ([(t2 _) (randlc/2 t2 t2)])
                 (loop t1 t2 ik (add1 i)))]))
        (cond 
          [(not (= (* 2 ik) kk))
           (let-values ([(t1 _) (randlc/2 t1 t2)])
             (finish t1 t2))]
          [else
            (finish t1 t2)])])))

(define (create-seq seed a max-key num-keys key-array)
  (define k (/ max-key 4))
  (for/fold ([seed seed])
            ([i num-keys])
    (let-values ([(seed x) (for/fold ([seed seed]
                                      [x 0])
                                     ([i 4])
                             (let-values ([(seed cx) (randlc/2 seed a)])
                               (values seed (+ x cx))))])
      (fx! key-array i (inexact->exact (floor (* k x))))
      seed)))

(define (rank iteration comm id cnt class max-iterations num-keys max-key num-buckets 
              key-array key-buff1 key-buff2 
              bucket-size bucket-ptrs
              process-bucket-distrib-ptr1 process-bucket-distrib-ptr2
              send-count recv-count send-displ recv-displ
              test-index-array test-rank-array
              test-array-size max-key-log-2 num-buckets-log-2)

  (define-syntax-rule (when-id-0 body ...)
    (when (= id 0)
      body ...))

  (timer-start T_RANK)

  (define lshift (- max-key-log-2 num-buckets-log-2))
  (define rshift (- (- max-key-log-2 num-buckets-log-2)))

  (when-id-0
    (fx! key-array iteration iteration) 
    (fx! key-array (fx+ iteration max-iterations) (fx- max-key iteration)))

  (for ([i (+ num-buckets test-array-size)])
    (fx! bucket-size i 0)
    (fx! process-bucket-distrib-ptr1 i 0)
    (fx! process-bucket-distrib-ptr2 i 0))

  (for ([i test-array-size])
    (when (= (quotient (fxr test-index-array i) num-keys) id)
      (fx! bucket-size (fx+ num-buckets i) 
           (fxr key-array (modulo
                            (fxr test-index-array i)
                            num-keys)))))

  (for ([i num-keys])
    (fx++! bucket-size (arithmetic-shift (fxr key-array i) rshift)))

  (fx! bucket-ptrs 0 0)
  (for/fold ([prev 0]) ([i (in-range 1 num-buckets)])
    (define nprev (+ prev (fxr bucket-size (fx- i 1))))
    (fx! bucket-ptrs i nprev)
    nprev)

  (for ([i num-keys])
    (define key (fxr key-array i))
    (define bkptr-index (arithmetic-shift key rshift))
    (define bkptr-val (fxr bucket-ptrs bkptr-index))
    (fx! bucket-ptrs bkptr-index (fx+ bkptr-val 1))
    (fx! key-buff1 bkptr-val key))

  (timer-stop T_RANK)
  (timer-start T_RCOMM)

  (define bucket-size-totals (rmpi-allreduce comm fx+ bucket-size))

  (timer-stop T_RCOMM)
  (timer-start T_RANK)

  (let-values ([(j bucket-sum-accumulator local-bucket-sum-accumulator)
    (for/fold ([j 0]
               [bucket-sum-accumulator 0]
               [local-bucket-sum-accumulator 0]) ([i num-buckets])
      (let ([bucket-sum-accumulator (fx+ bucket-sum-accumulator
                                         (fxr bucket-size-totals i))]
            [local-bucket-sum-accumulator (fx+ local-bucket-sum-accumulator
                                               (fxr bucket-size i))])
        (cond
          [(bucket-sum-accumulator . >= . (* num-keys (fx+ j 1)))
           (fx! send-count j local-bucket-sum-accumulator)
           (cond
             [(not (zero? j))
              (fx! send-displ j (+ (fxr send-displ (fx- j 1))
                                  (fxr send-count (fx- j 1))))
              (fx! process-bucket-distrib-ptr1 j (fx+ (fxr process-bucket-distrib-ptr2 (fx- j 1)) 1))])
           (fx! process-bucket-distrib-ptr2 j i)
           (values (fx+ j 1) bucket-sum-accumulator 0)]
          [else
            (values j bucket-sum-accumulator local-bucket-sum-accumulator)])))])

    (for ([jj (in-range j cnt)])
      (fx! send-count j 0)
      (fx! process-bucket-distrib-ptr1 j 1)))

  (timer-stop T_RANK)
  (timer-start T_RCOMM)

  (define recv-count (rmpi-alltoall comm send-count))

  (fx! recv-displ 0 0)
  (for/fold ([prev-displ 0]) ([i (in-range 1 cnt)])
    (define nprev (+ prev-displ (fxr recv-count (fx- i 1))))
    (fx! recv-displ i nprev)
    nprev)

  (rmpi-alltoallv comm key-buff1 send-count send-displ
                  key-buff2 recv-count recv-displ)

  (timer-stop T_RCOMM)
  (timer-start T_RANK)

  (define min-key-val (arithmetic-shift (fxr process-bucket-distrib-ptr1 id) lshift))
  (define max-key-val (fx- (arithmetic-shift (fx+ (fxr process-bucket-distrib-ptr2 id) 1) lshift) 1))

  (for ([i (in-range (fx+ (fx- max-key-val min-key-val) 1))])
    (fx! key-buff1 i 0))

  (let ([m (for/fold ([m 0]) ([k (in-range id)])
             (for/fold ([m m])
                       ([i (in-range (fxr process-bucket-distrib-ptr1 k)
                                     (fx+ (fxr process-bucket-distrib-ptr2 k) 1))])
             (fx+ m (fxr bucket-size-totals i))))]
        [j (for/fold ([j 0]) ([i (in-range (fxr process-bucket-distrib-ptr1 id)
                                           (fx+ (fxr process-bucket-distrib-ptr2 id) 1))])
             (fx+ j (fxr bucket-size-totals i)))])

    (for ([i (in-range j)])
      (fx++! key-buff1 (fx- (fxr key-buff2 i) min-key-val)))

    (for/fold ([prev (fxr key-buff1 min-key-val)])
              ([i (in-range (fx+ min-key-val 1) (fx+ max-key-val 1))])
      (define nprev (fx+ (fxr key-buff1 i) prev))
      (fx! key-buff1 i nprev)
      nprev)
    
    ;parial verify test
    (define (partial-verify class iteration key-buff1 bucket-size-totals test-rank-array min-key-val max-key-val)
      (for/and ([i test-array-size]) ;num-buckets (fx+ num-buckets test-array-size))]) 
        (define k (fxr bucket-size-totals (fx+ i num-buckets)))
        (if (and (<= min-key-val k) (<= k max-key-val)) 
            (let ([offset
              (case class
                [(#\S) (if (<= i 2) iteration (- 0 iteration))]
                [(#\W) (if (< i 2) (- iteration 2) (- 0 iteration))]
                [(#\A) (if (<= i 2) (- iteration 1) (+ (- 0 iteration) 1))]
                [(#\B) (if (or (= i 1) (= i 2) (= i 4)) iteration (- 0 iteration))]
                [(#\C) (if (<= i 2) iteration (- 0 iteration))]
                [(#\D) (if (< i 2) iteration (- 0 iteration))]
                [else (error "Unknown class")])])
              (let ([mh (fx+ (fxr key-buff1 (fx- (fx- k 1) min-key-val)) m)]
                    [th (fx+ (fxr test-rank-array i) offset)])
                (cond
                  [(not (= mh th))
                   (printf "Failed partial verification: iteration ~a, processor ~a test key ~a ~a ~a~n" iteration id i mh th)
                    #f]
                  [else #t])))
            #t)))

    (let ([x (partial-verify class iteration key-buff1 bucket-size-totals test-rank-array min-key-val max-key-val)])
      (timer-stop T_RANK)
      (values x min-key-val j))))

(define (full-verify comm id cnt min-key-val total-local-keys key-array key-buff1 key-buff2)
  (define total-lesser-keys 0)

  (timer-start T_VERIFY)

  (for ([i (in-range total-local-keys)])
    (define key-buff2-val (fxr key-buff2 i))
    (define key-buff-ptr-val (fx--! key-buff1 (fx- key-buff2-val min-key-val)))
    (fx! key-array (fx- key-buff-ptr-val total-lesser-keys) key-buff2-val))

  (define last-local-key (if (total-local-keys . < . 1) 0 (fx- total-local-keys 1)))

  (let*-values ([(k) (cond
                       [(= id 0)
                        (rmpi-send comm (fx+ id 1) (fxr key-array last-local-key))
                        #f]
                       [(= id (fx- cnt 1))
                        (rmpi-recv comm (fx- id 1))]
                       [else
                         (begin0
                           (rmpi-recv comm (fx- id 1))
                           (rmpi-send comm (fx+ id 1) (fxr key-array last-local-key)))])]
                [(j) (if (and (id . > . 0)
                              (total-local-keys . > . 0)
                              (k . > . (fxr key-array 0)))
                         1
                         0)]
                [(j prev)
                 (for/fold ([j j]
                            [prev (fxr key-array 0)])
                   ([i (in-range 1 total-local-keys)])
                   (define nprev (fxr key-array i))
                   (values 
                     (if (prev . > . nprev)
                         (fx+ j 1) 
                         j)
                     prev))
                 ])

    (timer-stop T_VERIFY)

    (cond
      [(= j 0) #t]
      [else
        (printf "Processor ~a  Full verify: number of keys out of sort ~a\n" id j)
        #f])))


(define/provide (is-place ch)
  (define-values (comm args tc) (rmpi-init ch))
  (match-define (list class bmname) args)
  (define id (rmpi-id comm))
  (define cnt (rmpi-cnt comm))
  (define-values (test-index-array test-rank-array total-keys-log-2 max-key-log-2 num-buckets-log-2 MIN-PROCS)
    (initial-conditions class))
  (define num-keys (arithmetic-shift 1 total-keys-log-2))
  (define max-key (arithmetic-shift 1 max-key-log-2))
  (define num-buckets (arithmetic-shift 1 num-buckets-log-2))
  (define test-array-size 5)
  (define max-iterations 10)

  (define SIZE-OF-BUFFERS
    (cond 
      [(cnt . <  . 256) (* 3/2 num-keys)]
      [(cnt . <  . 512) (* 5/2 num-keys)]
      [(cnt . <  . 1024) (* 4 num-keys)]
      [else (* 13/2 num-keys)]))

  (define key-array (make-fxvector SIZE-OF-BUFFERS))
  (define key-buff1 (make-fxvector SIZE-OF-BUFFERS))
  (define key-buff2 (make-fxvector SIZE-OF-BUFFERS))
  (define bucket-size (make-fxvector (+ num-buckets test-array-size)))
  (define bucket-ptrs (make-fxvector num-buckets))
  (define process-bucket-distrib-ptr1 (make-fxvector (+ num-buckets test-array-size)))
  (define process-bucket-distrib-ptr2 (make-fxvector (+ num-buckets test-array-size)))
  (define send-count (make-fxvector cnt))
  (define send-displ (make-fxvector cnt))
  (define recv-count (make-fxvector cnt))
  (define recv-displ (make-fxvector cnt))


  (define TOTAL-KEYS (arithmetic-shift 1 total-keys-log-2))
  ;(define key-array (make-fxvector num-keys 0))

  (create-seq (find-my-seed id cnt (* 4 TOTAL-KEYS MIN-PROCS) 314159265.00 1220703125.00)
              1220703125.00
              max-key
              num-keys
              key-array)

  (define (do-rank iteration)
    (rank iteration comm id cnt class max-iterations num-keys max-key num-buckets 
              ;fxvecotrs
              key-array key-buff1 key-buff2 
              bucket-size bucket-ptrs
              process-bucket-distrib-ptr1 process-bucket-distrib-ptr2
              send-count recv-count send-displ recv-displ

              test-index-array test-rank-array
              test-array-size max-key-log-2 num-buckets-log-2))

  (do-rank 1)

  (printf "Size: ~a Iterations: ~a~n" num-keys max-iterations) 
  (when (and (= id 0) (not (eq? class #\S))) (printf "\n   iteration\n"))

  (timer-start 0)

  (define-values (partial-verified min-key-val j)
    (for/fold ([partial-verified #t]
               [min-key-val 0]
               [j 0]) 
              ([iteration (in-range 1 (add1 max-iterations))])
      (when (and (= id 0) (not (eq? class #\S))) (printf "        ~a\n" iteration))
      (define-values (lpv lmkv lj) (do-rank iteration))
      (values (and partial-verified lpv) lmkv j)))

  (timer-stop 0)

  (define max-time (rmpi-reduce comm 0 max (timer-read 0)))

  (define (boolean-and a b) (and a b))

  (define verified (and
                     (rmpi-reduce comm 0 boolean-and (full-verify comm id cnt min-key-val j key-array key-buff2 key-buff2))
                     partial-verified))

  (when (= id 0)
      (print-verification-status class verified bmname)
      (let* ([tm-sec (/ (timer-read 1) 1000)]
             [results (new-BMResults "IS" 
                                      class
                                      TOTAL-KEYS
                                      0 
                                      0 
                                      max-iterations 
                                      (exact->inexact tm-sec)
                                      (get-mops tm-sec max-iterations num-keys) 
                                      "keys ranked" 
                                      verified
                                      #f 
                                      cnt ;num-threads
                                      0)])                             
        (print-results results)))
    
  (begin
    (define t1 
      (for/flvector ([i T_LAST])
        (timer-read i)))

    (define tmin (rmpi-reduce comm 0 min t1))
    (define tsum (rmpi-reduce comm 0 + t1))
    (define tmax (rmpi-reduce comm 0 max t1))

    (when (= id 0)
      (printf " nprocs = ~a          minimum     maximum     average\n" cnt)
      (for ([i T_LAST]
            [d '(total rcomp rcomm verify)])
        (printf " timer ~a (~a): ~a ~a ~a\n"
                (~r (fx+ i 1) #:min-width 2)
                (~a d #:width 8)
                (~r (flr tmin i) #:precision '(= 4) #:min-width 10)
                (~r (flr tmax i) #:precision '(= 4) #:min-width 10)
                (~r (/ (flr tsum i) cnt) #:precision '(= 4) #:min-width 10)
                ))
      (printf "\n")))


  
  (rmpi-finish comm tc))


(define (main . argv) 
  (let* ([args (parse-cmd-line-args argv "Integer Sort")]
         [class (BMArgs-class args)]
         [serial (BMArgs-serial args)]
         [num-threads (BMArgs-num-threads args)]
         [bmname "IS"])

    (print-banner "Integer Sort" args) 

    (define place-args (list class bmname))

    (rmpi-launch
      (rmpi-build-default-config
        #:mpi-module (quote-module-path)
        #:mpi-func 'is-place
        #:mpi-args place-args)
      (rmpi-make-localhost-config 4 6341 'is))))

(define (get-mops total-time niter num-keys)  
  (if (> total-time 0) 
      (* (+ niter num-keys) (/ niter (* total-time 1000000.0)))
      0.0))

(define-syntax-rule (v++! v idx)
  (vs! v idx (fx+ (vr v idx) 1)))

(define-syntax-rule (fx++! v idx)
  (fx! v idx (fx+ (fxr v idx) 1)))

(define-syntax-rule (fx--! v idx)
  (fx! v idx (fx- (fxr v idx) 1)))

(module+ test
         (find-my-seed 6 16 100 314159265.00 1220703125.00)
         (find-my-seed 7 16 100 314159265.00 1220703125.00)
         (find-my-seed 8 16 100 314159265.00 1220703125.00))
