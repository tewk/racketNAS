#lang racket/base

(provide main)
  
(require "bm-args.rkt") 
(require "bm-results.rkt") 
(require "rand-generator.rkt")
(require "timer.rkt")
(require "parallel-utils.rkt")
(require "debug.rkt")
(require (for-syntax racket/base))

(require (only-in racket/fixnum make-shared-fxvector))
(require (only-in racket/flonum make-shared-flvector))
#|
(require (only-in racket/fixnum fxvector-set! fxvector-ref fx+ fx- fx* fx=))
(require (only-in racket/flonum flvector-set! flvector-ref fl+ fl- fl* fl/ flsqrt))
(define vr vector-ref)
(define vs! vector-set!)
(define flvs! flvector-set!)
(define flvr flvector-ref)
(define fxvs! fxvector-set!)
(define fxvr fxvector-ref)
|#

(require (rename-in racket/unsafe/ops
                    [unsafe-fxvector-ref fxvr] 
                    [unsafe-fxvector-set! fxvs!]
                    [unsafe-flvector-ref flvr] 
                    [unsafe-flvector-set! flvs!]
                    [unsafe-flsqrt flsqrt]
                    [unsafe-fl+ fl+]
                    [unsafe-fl- fl-]
                    [unsafe-fl* fl*]
                    [unsafe-fl/ fl/]
                    [unsafe-fx+ fx+]
                    [unsafe-fx- fx-]
                    [unsafe-fx* fx*]
                    [unsafe-fx= fx=])
;         (rename-in racket/unsafe/ops [unsafe-fl+ fl+X])
        (only-in racket/flonum [fl+ fl+X])
)
#|
|#


;Constants 
(define amult 1220703125.0)
(define-syntax-rule (fx++ a) (fx+ a 1))

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values   1400  7  10 15 .1 8.5971775078648)]
    [(#\W) (values   7000  8  12 15 .1 10.362595087124)] 
    [(#\A) (values  14000 11  20 15 .1 17.130235054029)]
    [(#\B) (values  75000 13  60 75 .1 22.712745482631)]
    [(#\C) (values 150000 15 110 75 .1 28.973605592845)]
    [else (error "Unknown class")]))

(define (main . argv) 
  (let ([args (parse-cmd-line-args argv "Conjugate Gradient")]) 
    (run-benchmark args)))

(define cgitmax 25)

(define (run-benchmark args) 
  (let ([bmname "CG"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(na nonzer shift niter rcond zeta-verify-value ) (get-class-size CLASS)])
    (let* ([firstrow 1]
          [lastrow na]
          [firstcol 1]
          [lastcol na]
          [naa na]
          [nz (+ (* na (add1 nonzer) (add1 nonzer)) (* na (+ nonzer 2)))]
          [nzz nz]
          [colidx (make-shared-fxvector (+ nz 1) 0)]
          [rowstr (make-shared-fxvector (+ na 2) 0)]
          [iv (make-shared-fxvector (* 2 (+ na 2)) 0)]

          ;sparse matrix
          [arow (make-shared-fxvector (+ nz 1) 0)]
          [acol (make-shared-fxvector (+ nz 1) 0)]
          [aelt (make-shared-flvector (+ nz 1) 0.0)]
          [a (make-shared-flvector (+ nz 1) 0.0)]

          [v (make-shared-flvector (+ na 2) 0.0)]
          [p (make-shared-flvector (+ na 3) 0.0)]
          [q (make-shared-flvector (+ na 3) 0.0)]
          [r (make-shared-flvector (+ na 3) 0.0)]
          [x (make-shared-flvector (+ na 3) 0.0)]
          [z (make-shared-flvector (+ na 3) 0.0)]
          [rhomaster (make-shared-flvector num-threads 0.0)]
          [dmaster (make-shared-flvector num-threads 0.0)]
          [rnormmaster (make-shared-flvector num-threads 0.0)]
          [tnorm1master (make-shared-flvector num-threads 0.0)]
          [tnorm2master (make-shared-flvector num-threads 0.0)]
          [presults (make-shared-flvector 7 0.0)])
      (define ncols (add1 (- lastcol firstcol)))
      (define nrows (add1 (- lastrow firstrow)))

      (randlc amult)

      (makea firstrow lastrow firstcol lastcol naa nzz a colidx rowstr nonzer rcond arow acol aelt v iv shift)

      (for* ([j (in-range 1 (add1 (add1 (- lastrow firstrow))))]
             [k (in-range (fxvr rowstr j) (fxvr rowstr (add1 j)))])
        (fxvs! colidx k (add1 (- (fxvr colidx k) firstcol))))

      (let ([zeta
              (CGspawn (if serial 0 num-threads) cg-body
                niter nrows ncols naa shift colidx rowstr x z a p q r
                rhomaster dmaster rnormmaster tnorm1master tnorm2master presults)])

        (timer-stop 1)
        (let ([verified (verify zeta zeta-verify-value)])
          (print-banner "Conjugate Gradient" args) 
          (printf "Size = ~a niter = ~a~n" na niter) 
          (print-verification-status CLASS verified bmname)
          (let* ([time (/ (read-timer 1) 1000)]
                 [results (new-BMResults bmname CLASS na 0 0  niter time 
                                         (get-mflops na nonzer time niter)
                                         "floating point" 
                                         (if verified 1 0)
                                         serial 
                                         num-threads 
                                         -1)]) 
              (print-results results))))))))

(define (cg-body cg  
          niter nrows ncols naa shift colidx rowstr x z a p q r rhomaster 
          dmaster rnormmaster tnorm1master tnorm2master presults)
  (define serial (fx= 0 (CGnp cg)))

  (CG-n0-only cg 
    (for ([i (in-range 1 (fx+ naa 2))]) (flvs! x i 1.0)))

  (if serial
    (conj-grad nrows ncols naa shift colidx rowstr x z a p q r)
    (parallel-conj-grad cg -1 nrows ncols naa shift colidx rowstr x z a p q r
      rhomaster dmaster rnormmaster tnorm1master tnorm2master presults))

  (CG-n0-only cg 
    (for ([i (in-range 1 (fx+ naa 2))]) (flvs! x i 1.0))
    (timer-start 1))

  (begin0
    (for/last ([it (in-range niter)])
      (let-values ([(rnorm zeta)
        (if serial
          (conj-grad nrows ncols naa shift colidx rowstr x z a p q r)
          (values 
            (flvr presults 5)
            (parallel-conj-grad cg it nrows ncols naa shift colidx rowstr x z a p q r
              rhomaster dmaster rnormmaster tnorm1master tnorm2master presults)))])
      ;(printf "    ~a       ~a          ~a~n" it rnorm zeta)
      zeta))
    (CG-n0-only cg 
      (timer-stop 1))))

;//---------------------------------------------------------------------
;//       generate the test problem for benchmark 6
;//       makea generates a sparse matrix with a
;//       prescribed sparsity distribution
;//
;//       parameter    type        usage
;//
;//       input
;//
;//       n            i           number of cols/rows of matrix
;//       nz           i           nonzeros as declared array size
;//       rcond        r*8         condition number
;//       shift        r*8         main diagonal shift
;//
;//       output
;//
;//       a            r*8         array for nonzeros
;//       colidx       i           col indices
;//       rowstr       i           row pointers
;//
;//       workspace
;//
;//       iv, arow, acol i
;//       v, aelt        r*8
;//---------------------------------------------------------------------
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
      (sprnvc iouter n nn1 nonzer v iv colidx 0 colidx n)
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

(define (sprnvc iouter n nn1 nz v iv nzloc nzloc-offset mark mark-offset)
  (let loop ([nzv 0]
             [nzrow 0])
    (if ( nzv . >= . nz )
      (for ([ii (in-range 1 (add1 nzrow))])
        (fxvs! mark (+ (fxvr nzloc (+ ii nzloc-offset)) mark-offset) 0))
      (let* ([vecelt (randlc amult)]
             [vecloc (randlc amult)]
             [idx (inexact->exact (floor (add1 (* vecloc nn1))))])
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
  (if (not (for/fold ([set #f]) ([k (in-range 1 (add1 nzv))])
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

(define-syntax-rule (fxvs!++ v i) (fxvs! v i (fx++ (fxvr v i))))
(define-syntax-rule (flvs!+ v i val) (flvs! v i (fl+ (flvr v i) val)))
(define-syntax-rule (flvs!- v i val) (flvs! v i (fl- (flvr v i) val)))

(define (sparse firstrow lastrow a colidx rowstr n arow acol aelt x mark mark-offset nzloc nzloc-offset nnza)
  (define nrows (add1 (- lastrow firstrow)))
  (for ([j (in-range 1 (add1 n))])
    (fxvs! rowstr j 0)
    (fxvs! mark (+ j mark-offset) 0))
  (fxvs! rowstr (add1 n) 0)

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
    (let ([curval (fxvr rowstr j)])
      (fxvs! rowstr (add1 j) (fxvr rowstr j))))
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


(define (parallel-conj-grad cg it nrows ncols naa shift colidx rowstr x z a p q r rhomaster dmaster rnormmaster tnorm1master tnorm2master presults)
  (define-values (ALPHA BETA RHO TNORM1 TNORM2 RNORM ZETA) (values 0 1 2 3 4 5 6))
  (define id (CGid cg))
  (define np (CGnp cg))

    ;step 3
    (flvs! rhomaster id
      (CGfor/fold cg ([rho 0.0]) ([j (in-range 1 (fx++ naa))])
        (flvs! q j 0.0)
        (flvs! z j 0.0)
        (let ([xj (flvr x j)])
          (flvs! r j xj)
          (flvs! p j xj)
          (fl+X rho (fl* xj xj)))))
    (CG-n0-only cg (flvs! presults RHO (for/fold ([rho 0.0]) ([m (in-range np)]) (fl+ rho (flvr rhomaster m)))))


    (for ([cgit (in-range cgitmax)])
      ;step 0
      (CGfor cg ([j (in-range 1 (fx++ nrows))])
        (flvs! q j (for/fold ([sum 0.0]) ([k (in-range (fxvr rowstr j) (fxvr rowstr (fx++ j)))])
                    (fl+ sum (fl* (flvr a k) (flvr p (fxvr colidx k)))))))

      (flvs! dmaster id
        (CGfor/fold cg ([d 0.0]) ([j (in-range 1 (fx++ nrows))])
                    (fl+X d (fl* (flvr p j) (flvr q j)))))

      (CG-n0-only cg
        (flvs! presults ALPHA
          (fl/
            (flvr presults RHO)
            (for/fold ([d 0.0]) ([m (in-range np)]) (fl+ d (flvr dmaster m))))))
     
      ;step 1 
      (let ([alpha (flvr presults ALPHA)])
        (flvs! rhomaster id
          (CGfor/fold cg ([rho 0.0]) ([j (in-range 1 (fx++ ncols))])
            (let ([nrj (fl- (flvr r j) (fl* alpha (flvr q j)))])
              (flvs!+ z j (fl* alpha (flvr p j)))
              (flvs! r j nrj) 
              (fl+X rho (fl* nrj nrj))))))

      (CG-n0-only cg
        (let ([rho (for/fold ([rho 0.0])  ([m (in-range np)]) (fl+ rho (flvr rhomaster m)))])
          (flvs! presults BETA
            (fl/
              rho
              (flvr presults RHO)))
          (flvs! presults RHO rho)))


      ;step 2
      (let ([beta (flvr presults BETA)])
        (CGfor cg ([j (in-range 1 (fx++ ncols))])
          (flvs! p j (fl+ (flvr r j) (fl* beta (flvr p j))))))

      (CG-B cg))

    ;step 4
    (let-values ([(tnorm1 tnorm2 rnorm)
      (CGfor/fold cg ([tnorm1 0.0]
                      [tnorm2 0.0]
                      [rnorm  0.0]) ([j (in-range 1 (fx++ ncols))])
        (let* ([rj (for/fold ([sum 0.0]) ([k (in-range (fxvr rowstr j) (fxvr rowstr (fx++ j)))])
                    (fl+ sum (fl* (flvr a k) (flvr z (fxvr colidx k)))))]
               [zj (flvr z j)]
               [xj (flvr x j)]
               [xj-rj (fl- xj rj)])
          (values (fl+X tnorm1 (fl* xj zj))
                  (fl+X tnorm2 (fl* zj zj))
                  (fl+ rnorm (fl* xj-rj xj-rj)))))])
      (flvs! tnorm1master id tnorm1)
      (flvs! tnorm2master id tnorm2)
      (flvs! rnormmaster id rnorm))

    (CG-n0-only cg
      (let-values ([(tnorm1 tnorm2 rnorm)
          (for/fold ([tnorm1 0.0]
                     [tnorm2 0.0]
                     [rnorm 0.0])  ([m (in-range np)]) 
            (values (fl+ tnorm1 (flvr tnorm1master m))
                    (fl+ tnorm2 (flvr tnorm2master m))
                    (fl+ rnorm (flvr rnormmaster m))))])
        (flvs! presults RNORM (flsqrt rnorm))
        (flvs! presults TNORM2 (fl/ 1.0 (flsqrt tnorm2)))
        (flvs! presults ZETA (fl+ (exact->inexact shift) (fl/ 1.0 tnorm1)))))

    (let ([tnorm2 (flvr presults TNORM2)])
      (CGfor cg ([j (in-range 1 (fx++ ncols))])
        (flvs! x j (fl* tnorm2 (flvr z j)))))

    (CG-B cg)

    (flvr presults ZETA))

;;ilog2 : int -> int
(define (ilog2 n) 
  (if (= n 1) 
      0
      (let loop ([lg 1]
                 [nn 2])
        (if (< nn n)
          (loop (+ lg 1) (* 2 nn))
          lg))))

(define (get-mflops na nonzer total-time niter)
  (if (not (= total-time 0.0))
    (* (* 2 na (+ 3.0 (* nonzer (add1 nonzer)) (* 25.0 (+ 5.0 (* nonzer (add1 nonzer)))) 3))
       (/ niter (* total-time 1000000.0)))
    0.0))

(define (verify zeta zeta-verify-value)
  (define epsilon 1.0E-10)
  (printf " Zeta is    ~a~n" zeta)
  (let ([dev (abs (- zeta zeta-verify-value))])
    (if (dev . <= . epsilon)
      (begin
        (printf " Deviatation is   ~a~n" dev)
        #t)
      (begin
        (printf " The correct zeta is   ~a~n" zeta-verify-value)
        #f))))
