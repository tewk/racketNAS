#lang racket/base

(require racket/place/distributed
         racket/place/distributed/rmpi
         racket/vector
         racket/match
         racket/flonum
         racket/format
         "bm-args.rkt"
         "bm-results.rkt"
         (except-in "rand-generator.rkt" ipow46)
         "timer.rkt"
         "print-results.rkt")

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
                 
         (rename-in racket/fixnum [fxvector-ref fxr] [fxvector-set! fx!] [fxquotient fx/]
                    [fxvector-ref fxvr] [fxvector-set! fxvs!])
         (rename-in racket/flonum [flvector-ref flr] [flvector-set! fl!]
                    [flvector-ref flvr] [flvector-set! flvs!])
         racket/fixnum)

;; unsafe ones
#;
(require (rename-in racket/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-flvector-ref flr] 
                    [unsafe-flvector-set! fl!]
                    [unsafe-fxvector-ref fxr] 
                    [unsafe-fxvector-set! fx!]
                    [unsafe-fx+ fx+]
                    [unsafe-fx- fx-]
                    [unsafe-fx= fx=]
                    [unsafe-fx<= fx<=])

         (only-in racket/fixnum fxvector in-fxvector make-fxvector))

(provide main)

(define-values (T_TOTAL T_CONJG T_RCOMM T_NCOMM T_COMP T_COMM T_LAST) (values 0 1 2 3 4 5 6))

(define pi 3.141592653589793238)

(define-syntax-rule (when0 x body ...)
  (when (= 0 x)
    body ...))

(define-syntax-rule (with-timer T body ...)
  (begin
    (timer-start T)
    (begin0
      (let ()
        body ...)
      (timer-stop T))))

(define-syntax-rule (with-timer! T body ...)
  (begin
    (timer-start T)
    (begin0
      (let ()
        body ...)
      (timer-stop T))))

(define-syntax-rule (flv+= v o val)
  (fl! v o (fl+ (flr v o) val)))

(define-syntax-rule (flv-= v o val)
  (fl! v o (fl- (flr v o) val)))

(define-syntax-rule (fxvs!++ v o)
  (fx! v o (fx+1 (fxr v o))))

(define (rmpi-exit-error comm)
  (exit -1))

(define amult 1220703125.0)
(define-syntax-rule (fx++ a) (fx+ a 1))

(define (get-class-size class)
  (case class 
    [(#\S) (values    1400  7    10  15 .1 8.5971775078648)]
    [(#\W) (values    7000  8    12  15 .1 10.362595087124)] 
    [(#\A) (values   14000  11   20  15 .1 17.130235054029)]
    [(#\B) (values   75000  13   60  75 .1 22.712745482631)]
    [(#\C) (values  150000  15  110  75 .1 28.973605592845)]
    [(#\D) (values 1500000  21  500 100 .1 52.514532105794)]
    [(#\E) (values 9000000  26 1500 100 .1 77.522164599383)]
    [else (error "Unknown class")]))

(define (fx+1 x) (fx+ 1 x))

(define (ilog2 x)
  (cond
    [(<= x 0) #f]
    [else
      (let loop ([x x]
                 [i 0])
        (cond 
          [(= x 1) i]
          [(not (zero? (remainder x 2))) #f]
          [else (loop (quotient x 2) (+ 1 i))]))]))

(define (ipow2 i)
  (cond
    [(< i 0) #f]
    [(= i 0) 1]
    [else (arithmetic-shift 1 i)]))


(define (setup-proc-info comm id nprocs num-procs num-proc-rows num-proc-cols)
  (when (not (= nprocs num-procs))
    (when0 id 
      (printf "Error: num of procs allocated ~a is not equal to compiled number of procs ~a\n" nprocs num-procs))
    (rmpi-exit-error comm))

  (when (not (ilog2 num-proc-cols))
    (when0 id 
      (printf "Error: num-proc-cols is ~a which is not a power of two\n" num-proc-cols))
    (rmpi-exit-error comm))

  (when (not (ilog2 num-proc-rows))
    (when0 id 
      (printf "Error: num-proc-rows is ~a which is not a power of two\n" num-proc-rows))
    (rmpi-exit-error comm))

  (define log2nprocs (ilog2 nprocs))
  (when (not log2nprocs)
    (when0 id 
      (printf "Error: num-procs is ~a which is not a power of two\n" num-procs))
    (rmpi-exit-error comm))

  (values log2nprocs num-proc-cols num-proc-rows))

(define (setup-submatrix-info id naa nprows npcols reduce-exch-proc 
                              reduce-send-starts reduce-send-lengths reduce-recv-starts reduce-recv-lengths)

  (define col-size 0)
  (define row-size 0)

  ;/partit-size/
  ;ya set (define naa 0)
  ;ya set (define nzz 0)
  ;ya set (define npcols 0)
  ;ya set (define nprows 0)
  (define proc-row (fx/ id npcols))
  (define proc-col (fx- id (fx* proc-row npcols)))

  (define firstrow 0)
  (define lastrow 0)
  (define firstcol 0)
  (define lastcol 0)

  (define exch-proc 0)
  (define exch-recv-length 0)
  (define send-start 0)
  (define send-len 0)


  (cond 
    ; npcols = 1
    [(= (fx* (fx/ naa npcols) npcols) naa)
     (set! col-size (fx/ naa npcols))
     (set! firstcol (fx+ (fx* proc-col col-size) 1))
     (set! lastcol  (fx+ (fx- firstcol 1) col-size))
     (set! row-size (fx/ naa nprows))
     (set! firstrow (fx+ (fx* proc-row row-size) 1))
     (set! lastrow (fx+ (fx- firstrow 1) row-size))
     ;(printf/f "A1 ~a ~a ~a ~a ~a ~a\n" id npcols proc-col col-size firstcol lastcol)
    ]
    [else
      (define X (fx- naa (fx* (fx/ naa nprows) nprows)))
      (cond
        [(< proc-row X)
         (printf/f "~a (< proc-row X)\n" id)
         (set! row-size (fx+ (fx/ naa nprows) 1))
         (set! firstrow (fx+ (fx* proc-row row-size) 1))
         (set! lastrow  (fx+ (fx- firstrow 1) row-size))]
        [else
         (printf/f "~a NOT (< proc-row X)\n" id)
         (set! row-size (fx/ naa nprows))
         (set! firstrow (fx+ (fx+ (fx* X (fx+ row-size 1))
                                  (fx* (fx- proc-row X) row-size))
                             1))
         (set! lastrow  (fx+ (fx- firstrow 1) row-size))])

      (cond
        [(= npcols nprows)
         (printf/f "~a (= npcols nprows)\n" id)
         (define X (fx- naa (fx* (fx/ naa npcols) npcols)))
         (cond
           [(< proc-col X)
            (set! col-size (fx+ (fx/ naa npcols) 1))
            (set! firstcol (fx+ (fx* proc-col col-size) 1))
            (set! lastcol  (fx+ (fx- firstcol 1) col-size))]
           [else
            (set! col-size (fx/ naa npcols))
            (set! firstcol (fx+ (fx+ (fx* X (fx+ col-size 1))
                                     (fx* (fx- proc-col X) col-size))
                                1))
            (set! lastcol  (fx+ (fx- firstcol 1) col-size))])]
         [else
          (printf/f "~a NOT (= npcols nprows)\n" id)
           (define npcols/2 (fx/ npcols 2))
           (define proc-col/2 (fx/ proc-col 2))
           (define X (fx- naa (fx* (fx/ naa npcols/2) npcols/2)))
           (cond
             [(< (fx/ proc-col 2) X)
              (printf/f "~a (< (/ proc-col 2) X)\n" id)
              (set! col-size (fx+ (fx/ naa npcols/2) 1))
              (set! firstcol (fx+ (fx* proc-col/2 col-size) 1))
              (set! lastcol (fx+ (fx- firstcol 1) col-size))
              ]
             [else
              (printf/f "~a NOT (< (/ proc-col 2) X)\n" id)
              (set! col-size (fx/ naa npcols/2))
              (set! firstcol (fx+ (fx+ (fx* X (fx+ col-size 1))
                                       (fx* (fx- proc-col/2 X) col-size))
                                  1))
              (set! lastcol (fx+ (fx- firstcol 1) col-size))
              (printf/f "~a ~a ~a ~a ~a ~a\n" id npcols proc-col col-size firstcol lastcol)
               ])
           (cond
             [(= (modulo id 2) 0)
              (set! lastcol (fx+ (fx+ (fx- firstcol 1) (fx/ (fx- col-size 1) 2)) 1))]
             [else
              (set! firstcol (fx+ (fx+ firstcol (fx/ (fx- col-size 1) 2)) 1))
              (set! lastcol (fx+ (fx- firstcol 1) (fx/ col-size 2)))
               ])
           ])])

  (cond
    [(= npcols nprows)
     (set! send-start 1)
     (set! send-len (fx+ (fx- lastrow firstrow) 1))]
    [else
     (cond
       [(= (modulo id 2) 0)
        (set! send-start 1)
        (set! send-len (fx/ (fx+ 1 (fx+ (fx- lastrow firstrow) 1)) 2))]
       [else
        (set! send-start (fx+ (fx/ (fx+ 1 (fx+ (fx- lastrow firstrow) 1)) 2) 1))
        (set! send-len (fx/ (fx+ (fx- lastrow firstrow) 1) 2))
         ])])
  
  (cond
    [(= npcols nprows)
     (set! exch-proc (fx+ (fx* (modulo id nprows) nprows) (fx/ id nprows)))]
    [else
      (define id/2 (fx/ id 2))
     (set! exch-proc (fx+ (fx* (fx+ (fx* (modulo id/2 nprows) nprows) 
                                    (fx/ id/2 nprows))
                               2)
                          (modulo id 2)))])

  (define l2npcols (ilog2 npcols))

  (for/fold ([div-factor npcols]) ([i (in-range 1 (fx+1 l2npcols))])
    (define j (fx+ (modulo (fx+ proc-col (fx/ div-factor 2)) div-factor)
                   (fx* (fx/ proc-col div-factor) div-factor)))
    (fx! reduce-exch-proc i (fx+ (fx* proc-row npcols) j))
    (fx/ div-factor 2))

  (for ([i (in-range l2npcols 0 -1)])
    (cond
      [(= nprows npcols)
       (fx! reduce-send-starts i send-start)
       (fx! reduce-send-lengths i send-len)
       (fx! reduce-recv-lengths i (fx+ (fx- lastrow firstrow) 1))]
      [else
       (fx! reduce-recv-lengths i send-len)
       (cond
         [(= i l2npcols)
          (fx! reduce-send-lengths i (fx- (fx+ (fx- lastrow firstrow) 1) send-len))
          (cond
            [(= (fx* (fx/ id 2) 2) id)
              (fx! reduce-send-starts i (fx+ send-start send-len))]
            [else
              (fx! reduce-send-starts i 1)])]
         [else
           (fx! reduce-send-lengths i send-len)
           (fx! reduce-send-starts i send-start)])])
    (fx! reduce-recv-starts i send-start))

  (set! exch-recv-length (fx+ (fx- lastcol firstcol) 1))

  (values l2npcols firstrow lastrow firstcol lastcol exch-proc exch-recv-length send-start send-len)
)

(define (rmpi-send/flvo comm peer v o l)
  (rmpi-send comm peer (for/flvector ([i l]) (flr v (fx+ o i)))))

(define (rmpi-recv/flvo comm peer v o l)
  (define tv (rmpi-recv comm peer))
  ;(printf "RECV ~a ~a ~a\n" (flvector-length tv) l o)
  (for ([i l])
    (fl! v (fx+ o i) (flr tv i))))

(define-syntax-rule (two-way id peerid send recv)
  (cond
    [(< id peerid) send recv]
    [else recv send]))

(define (printfl->file fn v)
  (with-output-to-file fn #:exists 'replace (lambda ()
    (for ([i (in-naturals)]
          [x (in-flvector v)])
      (printf "~a ~a\n" 
              (~r i #:min-width 12) 
              (~r x #:precision 16 #:notation 'exponential)))))
  )

(define (conj-grad id comm naa nprows npcols firstrow lastrow firstcol lastcol l2npcols
                   colidx rowstr x z a p q r w 
                   reduce-exch-proc reduce-send-starts reduce-send-lengths
                   reduce-recv-starts reduce-recv-lengths
                   exch-proc exch-recv-length send-start send-len)
  (define cgitmax 25)
  (with-timer T_CONJG
    (for ([j (in-range 1 (fx+ 1 (fx/ naa nprows)))])
      (define xj (flr x j))
      (fl! q j 0.0)
      (fl! z j 0.0)
      (fl! r j xj)
      (fl! p j xj)
      (fl! w j 0.0))

    (define rho
      (with-timer T_RCOMM
        (for/fold ([sum (for/fold ([sum 0.0]) ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
                          (define rj (flr r j))
                          (fl+ sum (fl* rj rj)))]) 
                  ([i (in-range 1 (fx+1 l2npcols))])
          ;(printf "R ~a ~a ~a\n" id i sum)
          (define peerid (fxr reduce-exch-proc i))
          (rmpi-send comm peerid sum)
          (fl+ sum (rmpi-recv comm peerid)))))
    ;(printf "RHO ~a ~a\n" id rho)
    #;(when (= id 0)
      (printfl->file "A.dat" a))

    (for/fold ([rho rho]) ([cgit (in-range cgitmax)])
      #;(when (= id 0)
      (for ([i (in-naturals)]
            [x (in-flvector a)])
        (printf "~a ~a\n" i x)))
      ;q = a * p
      (for ([j (in-range 1 (fx+ (fx- lastrow firstrow) 2))])
        (flvs! w j (for/fold ([sum 0.0]) ([k (in-range (fxvr rowstr j) (fxvr rowstr (fx++ j)))])
                    (fl+ sum (fl* (flvr a k) (flvr p (fxvr colidx k))))))
        ;(printf "WJ ~a ~a ~a\n" id j (flvr w j))
        )
        
      (for ([i (in-range l2npcols 0 -1)])
        (with-timer T_RCOMM
          (define peerid (fxr reduce-exch-proc i))
          ;(printf "T1 ~a ~a ~a ~a ~a ~a\n" id peerid (fxr reduce-send-starts i) (fxr reduce-send-lengths i) (fxr reduce-recv-starts i) (fxr reduce-recv-lengths i))
          (two-way id peerid
            (rmpi-send/flvo comm peerid w (fxr reduce-send-starts i) (fxr reduce-send-lengths i))
            (rmpi-recv/flvo comm peerid q (fxr reduce-recv-starts i) (fxr reduce-recv-lengths i))))
        (for ([j (in-range send-start (fx+ send-start (fxr reduce-recv-lengths i)))])
          (flv+= w j (flr q j))))

      ;Exchange piece of q with transpose processor:
      (if (not (fx= l2npcols 0))
          (with-timer T_RCOMM
            (cond 
              [(= exch-proc id)
               (for ([x send-len])
                 (fl! q (fx+1 x) (flr w (fx+ send-start x))))]
              [else
                (two-way id exch-proc
                  (rmpi-send/flvo comm exch-proc w send-start send-len)
                  (rmpi-recv/flvo comm exch-proc q 1 exch-recv-length))]
                ))
          (for ([j (in-range 1 (fx+1 exch-recv-length))])
            (fl! q j (flr w j))))

      ;clear w for reuse
      (for ([j (in-range 1 (max (fx+ (fx- lastrow firstrow) 2)
                                (fx+ (fx- lastcol firstcol) 2)))])
        (fl! w j 0.0))

      ;obtain p * q
      (define d
        (for/fold ([sum (for/fold ([sum 0.0]) ([j (in-range 1  (fx+ (fx- lastcol firstcol) 2))])
                          (fl+ sum (fl* (flr p j) (flr q j))))])
                  ([i (in-range 1 (fx+1 l2npcols))])

          (fl+ 
            sum
            (with-timer T_RCOMM
              (rmpi-send comm (fxr reduce-exch-proc i) sum)
              (rmpi-recv comm (fxr reduce-exch-proc i))))))

      (define alpha (fl/ rho d))
      ;(printf "ALPHA ~a ~a ~a ~a\n" id alpha d rho)

      (define rho0 rho)

      ;z = z + (alpha * p)
      ;r = r - (alpha * q)
      (for ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
        (flv+= z j (fl* alpha (flr p j)))
        (flv-= r j (fl* alpha (flr q j))))

      ;rho = r.r
      (define rhon
        (for/fold ([sum (for/fold ([sum 0.0]) ([j (fx+ (fx- lastcol firstcol) 2)])
                          (define rj (flr r j))
                          (fl+ sum (fl* rj rj)))])
                  ([i (in-range 1 (fx+1 l2npcols))])
          (fl+
            sum
            (with-timer T_RCOMM
              (rmpi-send comm (fxr reduce-exch-proc i) sum)
              (rmpi-recv comm (fxr reduce-exch-proc i))))))

      (define beta (fl/ rhon rho0))

      (for ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
        (fl! p j (fl+ (flr r j) (fl* beta (flr p j)))))
      
      rhon)

    (for ([j (in-range 1 (fx+ (fx- lastrow firstrow) 2))])
      (fl! w j (for/fold ([sum 0.0]) ([k (in-range (fxr rowstr j) (fxr rowstr (fx+1 j)))])
                         (fl+ sum (fl* (flr a k) (flr z (fxr colidx k)))))))

    (for ([i (in-range l2npcols 0 -1)])
      (with-timer T_RCOMM
        (define peerid (fxr reduce-exch-proc i))
        (two-way id peerid
          (rmpi-send/flvo comm peerid w (fxr reduce-send-starts i) (fxr reduce-send-lengths i))
          (rmpi-recv/flvo comm peerid r (fxr reduce-recv-starts i) (fxr reduce-recv-lengths i))))
      (for ([j (in-range send-start (fx+ send-start (fxr reduce-recv-lengths i)))])
        (flv+= w j (flr r j))))

    (cond
      [(not (= l2npcols 0))
       (with-timer T_RCOMM
            (cond 
              [(= exch-proc id)
               (for ([x send-len])
                 (fl! r (fx+1 x) (flr w (fx+ send-start x))))]
              [else
         (two-way id exch-proc
           (rmpi-send/flvo comm exch-proc w send-start send-len)
           (rmpi-recv/flvo comm exch-proc r 1 exch-recv-length))]))]
      [else
        (for ([j (in-range 1 (fx+1 exch-recv-length))])
          (fl! r j (flr w j)))])

    ;obtain d with sum-reduce
    ;rnorm
    (sqrt 
      (for/fold ([sum (for/fold ([sum 0.0]) ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
                        (define d (fl- (flr x j) (flr r j)))
                        (fl+ sum (fl* d d)))])
                ([i (in-range 1 (fx+1 l2npcols))])
        (fl+
          sum
          (with-timer T_RCOMM
            (rmpi-send comm (fxr reduce-exch-proc i) sum)
            (rmpi-recv comm (fxr reduce-exch-proc i))))))
  ))

(define (makea n nz a colidx rowstr nonzer 
               firstrow lastrow firstcol lastcol 
               rcond arow acol aelt v iv shift lrandlc)
  (define (warn nnza nz iouter)
    (printf "Space for matrix elements exceeded in makea~n")
    (printf "nnza, nzmax = ~a, ~a~n" nnza nz)
    (printf " iouter = ~a~n" iouter)
    (exit 0))
  (define-syntax-rule (++ v) (set! v (add1 v)))
  (define ratio (expt rcond (/ 1.0 (exact->inexact n))))
  (define nn1 (let loop ([x 1])
                (if (< x n)
                  (loop (* 2 x))
                  x)))

  (for ([i (in-range 1 (add1 n))]) (fxvs! iv (+ n i) 0))

  (define nnza2 0)
  (let-values ([(size nnza)
    (for/fold ([size 1.0]
               [nnza 0]) ([iouter  (in-range 1 (add1 n))])
      (sprnvc n nn1 nonzer v colidx iv 1 iv (add1 n) lrandlc)
      (let* ([nzv (vecset n v colidx nonzer iouter 0.5) ]
             [nnza
        (begin 
        (for/fold ([nnza nnza]) ([ivelt (in-range 1 (add1 nzv))])
          (define jcol (fxvr colidx ivelt))
          (let ([nnza (if (and (jcol . >= . firstcol) (jcol . <= . lastcol))
            (for/fold ([nnza nnza]) ([ivelt1 (in-range 1 (add1 nzv))])
              (define irow (fxvr colidx ivelt1))
              (if (and (irow . >= . firstrow) (irow . <= . lastrow))
                (let ([nnza (add1 nnza)])
                  (++ nnza2)
                  (when (nnza2 . > . nz) (warn nnza2 nz iouter))
                  (fxvs! acol nnza2 jcol)
                  (fxvs! arow nnza2 irow)
                  (flvs! aelt nnza2 (* size (flvr v ivelt) (flvr v ivelt1)))
                  nnza)
                nnza))
            nnza)])
          nnza)))])
        (values (* size ratio) nnza)))])
    
    (let ([nnza (for/fold ([nnza nnza]) ([i (in-range firstrow (add1 lastrow))])
      (if (and (i . >= . firstcol) (i . <= . lastcol))
        (let ([nnza (add1 nnza)])
          (++ nnza2)
          (define iouter (+ n i))
          (when (nnza2 . > . nz) (warn nnza2 nz iouter))
          (fxvs! acol nnza2 i)
          (fxvs! arow nnza2 i)
          (flvs! aelt nnza2 (- rcond shift))
          nnza)
        nnza))])
      (sparse firstrow lastrow a colidx rowstr n arow acol aelt v iv 1 iv (add1 n) nnza2))))

(define (sparse firstrow lastrow a colidx rowstr n arow acol aelt x mark mark-offset nzloc nzloc-offset nnza)
  (define nrows (add1 (- lastrow firstrow)))
  (for ([j (in-range 1 (fx+1 n))])
    (fxvs! rowstr j 0)
    (fxvs! mark (fx+ j mark-offset) 0))
  (fxvs! rowstr (fx+1 n) 0)

  (for ([nza (in-range 1 (add1 nnza))])
    (let ([j (+ (- (fxvr arow nza) firstrow) 2)])
      (fxvs!++ rowstr j)))

  (fxvs! rowstr 1 1)
  (for/fold ([oldval 1])([j (in-range 2 (+ nrows 2))])
    (let ([newval (+ (fxvr rowstr j) oldval)])
      (fxvs! rowstr j newval)
      newval))
  
  (for ([nza (in-range 1 (add1 nnza))])
    (let* ([j (+ (- (fxvr arow nza) firstrow) 1)]
           [k (fxvr rowstr j)])
      (flvs! a k (flvr aelt nza))
      (fxvs! colidx k (fxvr acol nza))
      (fxvs!++ rowstr j) ))
    
  (for ([j (in-range nrows 0 -1)])
      (fxvs! rowstr (add1 j) (fxvr rowstr j)))
  (fxvs! rowstr 1 1)

  (for ([i (in-range 1 (add1 n))])
    (flvs! x i 0.0)
    (fxvs! mark (+ i mark-offset) 0))


  (for/fold ([nza 0]
             [jajp1 (fxvr rowstr 1)])
            ([j (in-range 1 (add1 nrows))])
    (let ([nzrow
      (for/fold ([nzrow 0])
                ([k (in-range jajp1 (fxvr rowstr (add1 j)))])
        (let* ([i (fxvr colidx k)]
               [nx (+ (flvr x i) (flvr a k))])
          (flvs! x i nx)
          (if (and (zero? (fxvr mark (+ i mark-offset))) (not (zero? nx)))
            (let ([nzrow (add1 nzrow)])
              (fxvs! mark (+ i mark-offset) 1)
              (fxvs! nzloc (+ nzrow nzloc-offset) i)
              nzrow)
            nzrow)))])

      (let ([nza
        (for/fold ([nza nza])
                  ([k (in-range 1 (add1 nzrow))])
          (let ([i (fxvr nzloc (+ k nzloc-offset))])
            (fxvs! mark (+ i mark-offset) 0)
            (let ([xi (flvr x i)])
              (flvs! x i 0.0)
              (if (not (zero? xi))
                (let ([nza (add1 nza)])
                  (flvs! a nza xi)
                  (fxvs! colidx nza i)
                  nza)
                nza))))])

        (let ([jajp1 (fxvr rowstr (add1 j))])
          (fxvs! rowstr (add1 j) (+ nza (fxvr rowstr 1)))
          (values nza jajp1))))))

(define (sprnvc n nn1 nz v iv nzloc nzloc-offset mark mark-offset lrandlc)
  (let loop ([nzv 0]
             [nzrow 0])
    (if ( nzv . >= . nz )
      (for ([ii (in-range 1 (add1 nzrow))])
        (fxvs! mark (+ (fxvr nzloc (+ ii nzloc-offset)) mark-offset) 0))
      (let* ([vecelt (lrandlc)]
             [vecloc (lrandlc)]
             [idx (inexact->exact (floor (fl+ 1.0 (* vecloc nn1))))])
        (if (and (idx . <= . n) (zero? (fxvr mark (+ idx mark-offset))))
          (let ([nzrow (add1 nzrow)]
                [nzv (add1 nzv)])
            (fxvs! mark (+ idx mark-offset) 1)
            (fxvs! nzloc (+ nzrow nzloc-offset) idx)
            (flvs! v nzv vecelt)
            (fxvs! iv nzv idx)
            (loop nzv nzrow))
          (loop nzv nzrow))))))

(define (vecset n v iv nzv ival val)
  (if (not (for/fold ([set #f]) ([k (in-range 1 (fx+1 nzv))])
              (if (= (fxvr iv k) ival)
                (begin 
                  (flvs! v k val)
                  #t)
                set)))
    (let ([nzv (add1 nzv)])
      (flvs! v nzv val)
      (fxvs! iv nzv ival)
      nzv)
    nzv))


(define/provide (cg-place ch)
  (define-values (comm args tc) (rmpi-init ch))
  (match-define (list class bmname nprocs2) args)
  (define id (rmpi-id comm))
  (define nprocs (rmpi-cnt comm))
  (define lognprocs (ilog2 nprocs))
  (when (not lognprocs)
    (when0 id
      (printf "num-threads ~a is not a power of 2\n" nprocs)
      (rmpi-exit-error comm)))
  (define-values (num-proc-cols num-proc-rows)
    (let* ([c (fx/ lognprocs 2)]
           [r (fx/ lognprocs 2)]
           [c (if (not (= (+ r c) lognprocs)) (+ 1 c) c)])
      (values (ipow2 c) (ipow2 r))))

  (define num-procs (* num-proc-cols num-proc-rows))
  ;(printf/f "HERE ~a ~a ~a ~a\n" id num-proc-cols num-proc-rows num-procs)

  (define-values (na nonzer shift niter rcond zeta-verify-value) 
                   (get-class-size class))


  (define naa na)
  (define nz (+ (* na (fx/ (fx+1 nonzer) num-procs) (fx+1 nonzer))
                nonzer
                (fx/ (* na (+ nonzer 2 (fx/ num-procs 256))) num-proc-cols)))
  (define nzz nz)
  (define colidx (make-fxvector (fx+1 nz)))
  (define rowstr (make-fxvector (fx+ na 2)))
  (define iv (make-fxvector (fx+ (fx* na 2) 2)))
  (define arow (make-fxvector (fx+1 nz)))
  (define acol (make-fxvector (fx+1 nz)))
  (define v (make-flvector (fx+ na 2)))
  (define aelt (make-flvector (fx+1 nz)))
  (define a (make-flvector (fx+1 nz)))
  (define _idx (fx+ (fx/ na num-proc-rows) 3))
  (define x (make-flvector _idx))
  (define z (make-flvector _idx))
  (define p (make-flvector _idx))
  (define q (make-flvector _idx))
  (define r (make-flvector _idx))
  (define w (make-flvector _idx))

  (define reduce-exch-proc (make-fxvector num-proc-cols))
  (define reduce-send-starts (make-fxvector num-proc-cols))
  (define reduce-send-lengths (make-fxvector num-proc-cols))
  (define reduce-recv-starts (make-fxvector num-proc-cols))
  (define reduce-recv-lengths (make-fxvector num-proc-cols))




  (when (= id 0)
   (printf "NAS Parallel Benchmarks 3.3 -- CG Benchmark\n")
   (printf "Class: ~a\n" class)
   (printf "Size: ~a\n" na)
   (printf "Iterations: ~a\n" niter)
   (printf "Number of active processes: ~a\n" nprocs)
   (printf "Number of nonzeroes per row: ~a\n" nonzer)
   (printf "Eigenvalue shift: ~a\n" shift)
   )

  ;(printf " naa ~a\n na ~a\n nz ~a\n nzz ~a\n" naa na nz nzz)

  (define-values (log2nprocs npcols nprows) (setup-proc-info comm id nprocs num-procs num-proc-rows num-proc-cols))
  ;(printf "log2nprocs ~a\n npcols ~a\n nprows ~a\n" log2nprocs npcols nprows)

  (define-values (l2npcols firstrow lastrow firstcol lastcol exch-proc exch-recv-length send-start send-len)
    (setup-submatrix-info id naa nprows npcols reduce-exch-proc
                          reduce-send-starts
                          reduce-send-lengths
                          reduce-recv-starts
                          reduce-recv-lengths))
  ;(printf " l2npcols ~a\n firstrow ~a\n lastrow ~a\n firstcol ~a\n lastcol ~a\n exch-proc ~a\n exch-recv-length ~a\n send-start ~a\n send-len ~a\n" l2npcols firstrow lastrow firstcol lastcol exch-proc exch-recv-length send-start send-len)

  (define ncols (add1 (- lastcol firstcol)))
  (define nrows (add1 (- lastrow firstrow)))

  (for ([i T_LAST]) (timer-clear i))

  (define tran 314159265.0)
  (define amult 1220703125.0)
  (define lrandlc (mk-randlc tran amult))
  (define zeta (lrandlc))

  ;Set up partition's sparse random matrix for given class size 
  (makea naa nzz a colidx rowstr nonzer 
         firstrow lastrow firstcol lastcol 
         rcond arow acol aelt v iv shift lrandlc)


  (for* ([j (in-range 1 (fx+ (fx- lastrow firstrow) 2))]
         [k (in-range (fxvr rowstr j) (fxvr rowstr (fx+1 j)))])
    (fxvs! colidx k (fx+1 (fx- (fxvr colidx k) firstcol))))

  (define (cg-body its silent)
    (for ([i (in-range 1 (fx+ (fx/ na num-proc-rows) 2))]) (fl! x i 1.0))
    (for/fold ([zeta 1.0])
              ([it its])
      (define rnorm (conj-grad id comm naa nprows npcols firstrow lastrow firstcol lastcol l2npcols
                 colidx rowstr x z a p q r w 
                 reduce-exch-proc reduce-send-starts reduce-send-lengths
                 reduce-recv-starts reduce-recv-lengths
                 exch-proc exch-recv-length send-start send-len))


      (let-values ([(norm-temp10 norm-temp11)
        (for/fold ([norm-temp10 0.0]
                   [norm-temp11 0.0])
                  ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
          (define zj (flr z j))
          (values (fl+ norm-temp10 (fl* (flr x j) zj))
                  (fl+ norm-temp11 (fl* zj zj))))])

        (let-values ([(norm-temp10 norm-temp11)
          (for/fold ([norm-temp10 norm-temp10]
                     [norm-temp11 norm-temp11])
                    ([i (in-range 1 (fx+1 l2npcols))])
            (define _ans
              (with-timer T_NCOMM
                (rmpi-send comm (fxr reduce-exch-proc i) (flvector norm-temp10 norm-temp11))
                (rmpi-recv comm (fxr reduce-exch-proc i))))
            (define norm-temp20 (flvector-ref _ans 0))
            (define norm-temp21 (flvector-ref _ans 1))
            (values (fl+ norm-temp10 norm-temp20)
                    (fl+ norm-temp11 norm-temp21)))])
           
           (define norm_temp11 (fl/ 1.00 (flsqrt norm-temp11)))

           (begin0
             (cond 
               [(and (not silent) (= id 0))
                 (define zeta (+ shift (fl/ 1.0 norm-temp10)))
                 (when (= it 0)
                   (printf "   iteration           ||r||                 zeta\n"))
                 (printf " ~a ~a ~a\n" it rnorm zeta)
                 zeta]
               [else zeta]
               )

             (for ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
               (fl! x j (fl* norm_temp11 (flr z j)))))))))

  (cg-body 1 #t)
  (for ([i T_LAST]) (timer-clear i))
  (rmpi-barrier comm)
  (let ([zeta
            (with-timer T_TOTAL 
                        (cg-body niter #f))])

    (define tmax (rmpi-reduce comm 0 max (timer-read 1)))

    (when (= id 0)
      (printf " Benchmark Complete\n")
      (define err (flabs (fl/ (fl- zeta zeta-verify-value) zeta-verify-value)))
      (define verified (fl< err 0.0000000001))
      (cond
        [(fl< err 0.0000000001)
          (printf " VERIFICATION SUCCESSFUL\n")
          (printf " Zeta                ~a\n" zeta)
          (printf " Error is            ~a\n" err)
          ]
        [else
          (printf " VERIFICATION FAILED\n")
          (printf " Zeta                ~a\n" zeta)
          (printf " The correct zeta is ~a\n" zeta-verify-value)])

      (define mflops
        (if (not (= tmax 0.0))
            (/ (/ (* (*  2 niter na)
                    (+ 3.0 (* nonzer (+ nonzer 1))
                       (* 25.0 (+ 5 (* nonzer (+ nonzer 1))))
                       3.0))
                  (/ tmax 1000))
               1000000.0)
            0.0))

      (print-results-fortran "CG" class na 0 0 niter nprocs nprocs (/ tmax 1000) mflops
                     "          floating point" verified 0.1)))


  (define t1 (for/flvector #:length T_LAST ([i T_LAST]) (timer-read i)))
  (fl! t1 T_CONJG (fl- (flr t1 T_CONJG) (flr t1 T_RCOMM)))
  (fl! t1 T_COMM  (fl+ (flr t1 T_RCOMM) (flr t1 T_NCOMM)))
  (fl! t1 T_COMP  (fl- (flr t1 T_TOTAL) (flr t1 T_COMM)))

  (define tsum (rmpi-reduce comm 0 + t1))
  (define tming (rmpi-reduce comm 0 min t1))
  (define tmaxg (rmpi-reduce comm 0 max t1))

  (when (= id 0)
    (printf " nprocs = ~a          minimum     maximum     average\n" nprocs)
    (for ([i T_LAST]
          [d '(total conjg rcomm ncomm totcomp totcomm)])
      (printf " timer ~a (~a): ~a ~a ~a\n"
              (~r (fx+ i 1) #:min-width 2)
              (~a d #:width 8)
              (~r (flr tming i) #:precision '(= 4) #:min-width 10)
              (~r (flr tmaxg i) #:precision '(= 4) #:min-width 10)
              (~r (/ (flr tsum i) nprocs) #:precision '(= 4) #:min-width 10)
              ))
    (printf "\n"))
          
  (rmpi-finish comm tc))

(define (main . argv) 
  (let* ([args (parse-cmd-line-args argv "Conjugate Gradient")]
         [class (BMArgs-class args)]
         [serial (BMArgs-serial args)]
         [num-threads (BMArgs-num-threads args)]
         [bmname "CG"])

    (print-banner "Conjugate Gradient" args) 

    (define place-args (list class bmname num-threads))

    (rmpi-launch
      (rmpi-build-default-config
        #:mpi-module (quote-module-path)
        #:mpi-func 'cg-place
        #:mpi-args place-args)
      (rmpi-make-localhost-config num-threads 6341 'cg))))

(define-syntax-rule (v++! v idx)
  (vs! v idx (fx+ (vr v idx) 1)))

(define-syntax-rule (fx++! v idx)
  (fx! v idx (fx+ (fxr v idx) 1)))

(define-syntax-rule (fx--! v idx)
  (fx! v idx (fx- (fxr v idx) 1)))
