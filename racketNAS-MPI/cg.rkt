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
                 
         (rename-in racket/fixnum [fxvector-ref fxr] [fxvector-set! fx!] [fxquotient fx/])
         (rename-in racket/flonum [flvector-ref flr] [flvector-set! fl!])
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

(define pi 3.141592653589793238)

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
    [(= x 1) 0]
    [else]
      (let loop ([x x]
                 [i 0])
        (cond 
          [(zero? x) i]
          [(not (zero? (remainder x 2))) #f]
          [else (loop (quotient x 2) (+ 1 i))]))))

(define (ispow2 i)
  (cond
    [(< i 0) #f]
    [(= i 0) 1]
    (let loop ([i i]
               [x 1])
      (cond 
        [(zero? i) x]
        [else (loop (- i 1) (arithmetic-shift x 1))]))))


(define (setup-proc-info num-procs num-proc-rows num-proc-cols)
  (when (not (= nprocs num-procs))
    (when0 id 
      (printf "Error: num of procs allocated ~a is not equal to compiled number of procs ~a\n" nprocs num-procs))
    (rmpi-exit-error comm))

  (when (not (ispow2 num-proc-cols))
    (when0 id 
      (printf "Error: num-proc-cols is ~a which is not a power of two\n" num-proc-cols))
    (rmpi-exit-error comm))

  (define log2nprocs (ispow2 nprocs))
  (when (not log2nprocs)
    (when0 id 
      (printf "Error: num-procs is ~a which is not a power of two\n" num-procs))
    (rmpi-exit-error comm))

  (values log2nprocs num-proc-cols num-proc-rows))

(define (setup-submatric-info l2npcols reduce-exch-proc reduce-send-starts redcue-send-lengths reduce-recv-starts reduce-recv-lengths)

  (define col-size 0)
  (define row-size 0)

  ;/partit-size/
  ;ya set (define naa 0)
  ;ya set (define nzz 0)
  ;ya set (define npcols 0)
  ;ya set (define nprows 0)
  (define proc-col (fx/ id npcols))
  (define proc-row (fx- id (fx* proc-row npcols)))

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
    [(= (fx/ naa (fx* npcols npcols)) naa)
     (set! col-size (fx/ naa npcols))
     (set! firstcol (fx+ (fx* proc-col col-size) 1))
     (set! lastcol  (fx+ (fx- firstcol 1) col-size))
     (set! firstrow (fx+ (fx* proc-row row-sze)) 1)
     (set! lastrow (fx+ (fx- firstrow 1) row-size))
    ]
    [else
      (define X (fx- naa (fx/ naa (fx* nprows nprows))))
      (cond
        [(< proc-row X)
         (set! row-size (fx+ (fx/ naa nprows) 1))
         (set! firstrow (fx+ (fx* proc-row row-size) 1))
         (set! lastrow  (fx+ (fx- firstrow 1) row-size))]
        [else
         (set! row-size (fx/ naa nprows))
         (set! firstrow (fx+ (fx* X (fx+ row-size 1))
                             (fx* (fx- proc-row X) row-size)
                             1))
         (set! lastrow  (fx+ (fx- firstrow 1) row-size))])

      (cond
        [(= npcols nprows)
         (cond
           [(< proc-col X)
            (set! col-size (fx+ (fx/ naa npcols) 1))
            (set! firstcol (fx+ (fx* proc-col col-size) 1))
            (set! lastcol  (fx+ (fx- firstcol 1) col-size))]
           [else
            (set! col-size (fx/ naa npcols))
            (set! firstcol (fx+ (fx* X (fx+ col-size 1))
                                (fx* (fx- proc-col X) col-size)
                                1))
            (set! lastcol  (fx+ (fx- firstcol 1) col-size))])
         [else
           (error "NOT IMPLEMENTED (= npcols nprows) must be true")]
           ])])

  (cond
    [(= npcols nprows)
     (set! send-start 1)
     (set! send-len (fx+ (fx- lastrow firstrow) 1))]
    [else
      (error "NOT IMPLEMENTED (= npcols nprows) must be true")])

  
  (cond
    [(= npcols nprows)
     (set! exch-proc ((fx* (modulo d nprows) nprows) (fx/ id nprows)))]
    [else
      (error "NOT IMPLEMENTED (= npcols nprows) must be true")])

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
        (error "NOT IMPLEMENTED (= npcols nprows) must be true")])
    (fx! reduce-recv-length (fx+ (fx- lastcol firstcol) 1)))

  (set! exch-recv-length (fx+ (fx- lastcol firstcol) 1))
)

(define (rmpi-send/flvo comm peer v o l)
  (rmpi-send comm peer (for/flvector ([i l]) (flr v (fx+ o i)))))

(define (rmpi-recv/flvo comm peer v o l)
  (define tv (rmpi-recv comm peer))
  (for ([i l])
    (fl! v (fx+ o i) (flr tv i))))

(define (conj-grad colidx rowstr x z a p q r w rnorm lp2ncols
                   reduce-exch-proc reduce-send-starts reduce-send-lengths
                   reduce-recv-starts reduce-recv-lengths)
  (with-timer T_CONJG
    (for ([j (in-range 1 (fx+ 1 (fx/ naa nprows)))])
      (define xj (flr x j))
      (fl! q j 0.0)
      (fl! z j 0.0)
      (fl! r j xj)
      (fl! p j xj)
      (fl! w j 0.0))

    (define sum (for/fold ([sum 0.0]) ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
      (define rj (flr r j))
      (fl+ sum (fl* rj rj))))

    (define rho
      (with-timer T_RCOMM
        (for/fold ([sum sum]) ([i (in-range 1 (fx+1 l2npcols))])
          (define peerid (fxr reducde-exch-proc i))
          (rmpi-send comm peerid sum)
          (fl+ sum (rmpi-recv comm peerid)))))

    (for/fold ([rho rho]) ([cgit (in-range cgitmax)])
      ;q = a * p
      (for ([j (in-range 1 (fx+ (fx- lastrow firstrow) 2))])
        (flvs! w j (for/fold ([sum 0.0]) ([k (in-range (fxvr rowstr j) (fxvr rowstr (fx++ j)))])
                    (fl+ sum (fl* (flvr a k) (flvr p (fxvr colidx k)))))))
      (for ([i (in-range 1 (fx+1 l2npcols))])
        (with-timer T_RCOMM
          (define peerid (fxr reducde-exch-proc i))
          (rmpi-send comm peerid (flr w (fxr reduce-send-starts i)))
          (fl! q (fxr reduce-recv-starts i) (rmpi-recv comm peerid)))
        (for ([j (in-range send-start (fx+ send-start (fxr reduce-recv-lengths i)))])
          (flv+ w j (flr q j))))

      ;Exchange piece of q with transpose processor:
      (if (not (fx= l2npcols 0))
          (with-timer T_RCOMM
            (rmpi-send comm exch-proc (fx-subvec w send-start send-length))
            (rmpi-recv/flvo comm exch-proc q 0 exch-recv-length ))
          (for ([j 1 (fx+1 exch-recv-length)])
            (fl! q j (flr w j))))

      ;clear w for reuse
      (for ([j (in-range 1 (max (fx+ (fx- lastrow firstrow) 2)
                                (fx+ (fx- lastcol firstcol) 2)))])
        (fl! w j 0.0))

      ;obtain p * q
      (define alpha (fl/ rho
        (for/fold ([sum (for/fold ([sum 0.0]) ([j (in-range 1  (fx+ (fx- lastcol firstcol) 2))])
                          (fl+ sum (fl* (flrp j) (flr q j))))])
                  ([i (in-range 1 (fx+1 l2npcls))])

          (fl+ 
            sum
            (with-timer T_RCOMM
              (rmpi-send comm (fxr reduce-exch-proc i) sum)
              (rmpi-recv comm (fxr reduce-exch-proc i)))))))

      (define rho0 rho)

      ;z = z + (alpha * p)
      ;r = r - (alpha * q)
      (for ([j (in-range 1 (fx+ (fx- lastcol firstcl) 2))])
        (fl+= z j (fl* alpha (flv p j)))
        (fl-= r j (fl* alpha (flv q j))))

      ;rho = r.r
      (define rho
        (for/fold ([sum (for/fold ([sum 0.0]) ([j (fx+ (fx- lastcol firstcol) 2)])
                          (define rj (flr r j))
                          (fl+ (fl* rj rj)))])
          (fl+
            sum
            (with-timer T_RCOMM
              (rmpi-send comm (fxr reduce-exch-proc i) sum)
              (rmpi-recv comm (fxr reduce-exch-proc i))))))

      (define beta (fl/ rho rho0))

      (for ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
        (fl! p j (fl+ (flr r j) (fl* beta (flr p j))))))

    (for ([j (in-range 1 (fx+ (fx- lastrow firstrow) 2))])
      (fl! w j (for/fold ([sum 0.0]) ([k (in-range (fxr rowstr j) (fxr rowstr (fx+1 j)))])
                         (fl+ sum (fl* (flr a k) (flr z (fxr colidx k)))))))

    (for ([i (in-range l2npcols 0 -1)])
      (with-timer T_RCOMM
        (rmpi-send/flvo comm (fxr reduce-exch-proc i) w (fxr reduce-send-starts i) (fxr reduce-send-lengths i))
        (rmpi-recv/flvo comm (fxr reduce-exch-proc i) r (fxr reduce-recv-starts i) (fxr reduce-recv-lengths i)))
      (for ([j (in-range send-start (fx+ send-start (fxr reduce-recv-lengths i)))])
        (fl+= w j (flr r j))))

    (cond
      [(not (= l2npcols 0))
       (with-timer T_RCOMM
         (rmpi-send/flvo comm exch-proc w send-start send-len)
         (rmpi-recv/flvo comm exhc-proc r 0 exch-recv-length))]
      [else
        (for ([j (in-range 1 (fx+1 exch-recv-length))])
          (fl! r j (flr w j)))])

    (sqrt 
      (for/fold ([sum (for/fold ([sum 0.0]) ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
                        (define d (fl- (flr x j) (flr r j)))
                        (fl+ sum (fl* d d)))])
        (fl+
          sum
          (with-timer T_RCOMM
            (rmpi-send comm (fxr reduce-exch-proc i) sum)
            (rmpi-recv-comm (fxr reduce-exch-proc i))))))
  ))

(define (conj-grad nrows ncols naa shift colidx rowstr x z a p q r)
  
  (define rho (for/fold ([rho 0.0]) ([j (in-range 1 (fx+ 1 naa))])
    (flvs! q j 0.0)
    (flvs! z j 0.0)
    (let ([xj (flvr x j)])
      (flvs! r j xj)
      (flvs! p j xj)
      (fl+ rho (fl* xj xj)))))

  (for/fold ([rho rho])  ([cgit (in-range cgitmax)])
    (for ([j (in-range 1 (fx++ nrows))])
      (flvs! q j (for/fold ([sum 0.0]) ([k (in-range (fxvr rowstr j) (fxvr rowstr (fx++ j)))])
                  (fl+ sum (fl* (flvr a k) (flvr p (fxvr colidx k)))))))

    (define d (for/fold ([d 0.0]) ([j (in-range 1 (fx++ ncols))])
                  (fl+ d (fl* (flvr p j) (flvr q j)))))

    (define alpha (fl/ rho d))
                    
    (for ([j (in-range 1 (fx++ ncols))])
      (flvs!+ z j (fl* alpha (flvr p j)))
      (flvs!- r j (fl* alpha (flvr q j))))

    (define rhon (for/fold ([rho 0.0]) ([j (in-range 1 (fx++ ncols))])
                  (let ([rj (flvr r j)])
                    (fl+ rho (fl* rj rj)))))

    (define beta (fl/ rhon rho))

    (for ([j (in-range 1 (fx++ ncols))])
      (flvs! p j (fl+ (flvr r j) (fl* beta (flvr p j)))))
    rhon)

  (for ([j (in-range 1 (fx++ nrows))])
    (flvs! r j (for/fold ([sum 0.0]) ([k (in-range (fxvr rowstr j) (fxvr rowstr (fx++ j)))])
                (fl+ sum (fl* (flvr a k) (flvr z (fxvr colidx k)))))))

  (define rnorm (flsqrt (for/fold ([sum 0.0]) ([j (in-range 1 (fx++ ncols))])
    (let ([xj-rj (fl- (flvr x j) (flvr r j))])
      (fl+ sum (fl* xj-rj xj-rj))))))

  (let-values ([(tnorm1 tnorm2)
    (let-values ([(tnorm1 tnorm2)
        (for/fold ([tnorm1 0.0]
                   [tnorm2 0.0])
                ([j (in-range 1 (fx++ ncols))])
          (let ([zj (flvr z j)])
            (values (fl+ tnorm1 (fl* (flvr x j) zj))
                    (fl+ tnorm2 (fl* zj zj)))))])
        (values tnorm1 (fl/ 1.0 (flsqrt tnorm2))))])

    (for ([j (in-range 1 (fx++ ncols))])
      (flvs! x j (fl* tnorm2 (flvr z j))))
    (values rnorm (fl+ (exact->inexact shift) (fl/ 1.0 tnorm1)))))


;I think I switched the iv and colidx vectors by accident
(define (makea firstrow lastrow firstcol lastcol n nz a colidx rowstr nonzer rcond arow acol aelt v iv shift)
  (define (warn nnza nz iouter)
    (printf "Space for matrix elements exceeded in makea~n")
    (printf "nnza, nzmax = ~a, ~a~n" nnza nz)
    (printf " iouter = ~a~n" iouter)
    (exit 0))
  (define-syntax-rule (++ v) (set! v (add1 v)))
  (define ratio (expt rcond (/ 1.0 (exact->inexact n))))
  (define nn1 (expt 2 (ilog2 n)))

  (for ([i (in-range 1 (add1 n))]) (fxvs! colidx (+ n i) 0))

  (define nnza2 0)
  (let-values ([(size nnza)
    (for/fold ([size 1.0]
               [nnza 0]) ([iouter  (in-range 1 (add1 n))])
      (sprnvc n nn1 nonzer v iv colidx 0 colidx n)
      (let* ([nzv (vecset n v iv nonzer iouter 0.5) ]
             [nnza
        (begin 
        (for/fold ([nnza nnza]) ([ivelt (in-range 1 (add1 nzv))])
          (define jcol (fxvr iv ivelt))
          (let ([nnza (if (and (jcol . >= . firstcol) (jcol . <= . lastcol))
            (for/fold ([nnza nnza]) ([ivelt1 (in-range 1 (add1 nzv))])
              (define irow (fxvr iv ivelt1))
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
      (sparse firstrow lastrow a colidx rowstr n arow acol aelt v iv 0 iv n nnza2))))

;(define (sparse a collidx rowstr n arow acol aelt firstrow lastrow x mark mark-offset nzloc nzloc-offset nnza)
)
(define (sparse firstrow lastrow a collidx rowstr n arow acol aelt x mark mark-offset nzloc nzloc-offset nnza)
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
    
  (for ([j (in-range (sub1 nrows) -1 -1)])
    ;(let ([curval (fxvr rowstr j)])
      (fxvs! rowstr (add1 j) (fxvr rowstr j)))
    ;)
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

      (sprnvc n nn1 nonzer v iv colidx 0 colidx n)
(define (sprnvc n nn1 nz v iv nzloc nzloc-offset mark mark-offset)
  (let loop ([nzv 0]
             [nzrow 0])
    (if ( nzv . >= . nz )
      (for ([ii (in-range 1 (add1 nzrow))])
        (fxvs! mark (+ (fxvr nzloc (+ ii nzloc-offset)) mark-offset) 0))
      (let* ([vecelt (randlc amult)]
             [vecloc (randlc amult)]
             [idx (inexact->exact (floor (fx+1 (* vecloc nn1))))])
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
  (match-define (list class bmname nprocs) args)
  (define id (rmpi-id comm))
  (define nprocs (rmpi-cnt comm))
  (define lognprocs (ilog2 nprocs))
  (when (not lognprocs)
    (when0 id
      (printf "num-threads ~a is not a power of 2\n" nprocs)
      (rmpi-exit-error comm))
  (define num-proc-cols (ipow2 (- lognprocs 1)))
  (define num-proc-rows (ipow2 (- lognprocs 1)))
  (define num-procs (* num-proc-cols num-proc-rows))

  (define-values (na nonzer shift niter rcond zeta-verify-value) 
                   (get-class-size CLASS))


  (define naa na)
  (define nz (+ (* na (fx/ (fx+1 nonzer) num-procs) (fx+1 nonzer))
                nonzer
                (fx/ (* na (+ nonzer 2 (fx/ num-procs 256)) num-proc-cols))))
  (define nzz nz)


  (define ncols (add1 (- lastcol firstcol)))
  (define nrows (add1 (- lastrow firstrow)))

  (when (= id 0)
   (printf "NAS Parallel Benchmarks 3.3 -- CG Benchmark")
   (printf "Size: ~a" na)
   (printf "Iterations: ~a" niter)
   (printf "Number of active processes: ~a" nprocs)
   (printf "Number of nonzeroes per row: ~a" nonzero)
   (printf "Eigenvalue shift: ~a" shift)
   )

  ;(define naa aa)
  ;(define nzz nz)

  (define-values (log2nprocs npcols nprows) (setup-proc-info num-procs num-proc-rows num-proc-cols))

  (setup-sumbatrix-info l2npcols reduce-exch-proc
                        reduce-exch-starts
                        reduce-send-starts
                        reduce-send-lengths
                        reduce-recv-starts
                        reduce-recv-lengths)

  (for ([i T_LAST]) (timer-clear i))

  (define tran 314159265.0)
  (define amult 1220703125.0)
  (define zeta (randlc tran amult))

  ;Set up partition's sparse random matrix for given class size 
  (makea naa nzz a colidx rowstr nonzer firstrow lastrow firstcol lastcol
         rcond arow acol aelt v iv shift)


  (for* ([j (in-range 1 (fx+ (fx- lastrow firstrow) 2))]
         [k (in-range (fxvr rowstr j) (fxvr rowstr (fx+1 j)))])
    (fxvs! colidx k (fx+1 (fx- (fxvr colidx k) firstcol))))

  (define (cg-body its silent)
    (for ([i (in-range 1 (fx+ (fx/ na num-proc-rows) 2))]) (fl! x i 1.0))
    (for/fold ([zeta 0.0])
              ([it its])
      (conj-grad colidx rowstr x z a p q r w rnomr l2npcols
                 reduce-exch-proc reduce-send-starts reduce-send-lengths
                 reduce-recv-starts reduce-recv-lengths)


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
            (match-define (fl-vector norm-temp20 norm-temp21) 
              (with-timer T_NCOMM
                (rmpi-send comm i (flvector norm-temp10 norm-temp11))
                (rmpi-recv comm i)))
            (values (fl+ norm-temp10 norm-temp20)
                    (fl+ norm-temp11 norm-temp21)))])
           
           (define norm_temp11 (fl/ 1.00 (flsqrt norm-temp11)))

           (when (and (not silent) (= id 0))
             (define zeta (fl+ shift (fl/ 1.0 norm_temp11)))
             (when (and)(= it 0)
               (printf "   iteration           ||r||                 zeta\n")
               (printf " ~a ~a ~a\n" it rnorm zeta))
             )

           (for ([j (in-range 1 (fx+ (fx- lastcol firstcol) 2))])
             (fl! x j (fx* norm_temp11 (flr z j))))))))

  (cg-body 1 #t) 
  (for ([i T_LAST]) (timer-clear i))
  (rmpi-barrier comm)
  (with-timer 1
    (cg-body niter #f))

  (define t (rmpi-reduce comm 0 max (read-timer 1)))

  (when (= id 0)
    (printf " Benchmark Complete\n")
    (if (fl< (flabs (fl/ (fl- zet zeta-verify-value) zeta-verify-value)) 0.0000000001)
        (printf " VERIFICATION SUCCESSFUL\n")
        (printf " VERIFICATION FAILED\n"))
    (printf " Zeta                ~a" zeta)
    (printf " The correct zeta is ~a" zeta)

    (define mflops
      (if (not (= tmax 0.0))
          (/ (/ (* (*  2 niter na)
                  (+ 3.0 (* nonzer (+ nonzer 1))
                     (* 25.0 (+ 5 (* nonzer (+ nonzer 1))))
                     3.0))
                tmax)
             1000000.0)
          0.0))

    (print-results "CG" class na 0 0 niter nnodes-compiled nprocs tmax mflops
                   "          floating point" verfied npbversion compiletime
                   cs1 cs2 cs3 cs4 cs5 cs6 cs7))


  (define t1 (for/flvector #:length T_LAST ([i T_LAST]) (read-timer i)))
  (fl! t1 T_CONFIG (fl- (flr t1 T_CONFIG) (flr t1 T_RCOMM)))
  (fl! t1 T_COMM   (fl+ (flr t1 T_RCOMM) (flr t1 T_NCOMM)))
  (fl! t1 T_COMP   (fl- (flr t1 T_TOTAL) (flr t1 T_COMM)))

  (define tsum (rmpi-reduce comm 0 + t1))
  (define tming (rmpi-reduce comm 0 min t1))
  (define tmaxg (rmpi-reduce comm 0 max t1))

  (when (= id 0)
    (printf " nprocs = ~a          minimum     maximum     average\n" cnt)
    (for ([i T_LAST]
          [d '(total conjg rcomm ncomm totcomp totcomm)])
      (printf " timer ~a (~a): ~a ~a ~a\n"
              (~r (fx+ i 1) #:min-width 2)
              (~a d #:width 8)
              (~r (flr tmin i) #:precision '(= 4) #:min-width 10)
              (~r (flr tmax i) #:precision '(= 4) #:min-width 10)
              (~r (/ (flr tsum i) cnt) #:precision '(= 4) #:min-width 10)
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
