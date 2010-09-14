#lang racket/base

(provide main)
  
(require "../bm-args.ss") 
(require "../bm-results.ss") 
(require "../rand-generator.ss")
(require "../timer.ss")
(require racket/future)
(require racket/list)
(require racket/match)
(require racket/math)
(require racket/place)
(require racket/place-utils)

#;(require scheme/fixnum scheme/flonum)

(require (only-in scheme/flonum make-flvector make-shared-flvector)
         scheme/require (for-syntax scheme/base)
   (filtered-in
    (lambda (name) (regexp-replace #rx"unsafe-" name ""))
    scheme/unsafe/ops))

(define (comgrp-wait grp)
  (match grp
    [(list 0 pls np) (for ([ch pls]) (place-channel-send ch 0))]
    [(list n ch np)  (place-channel-recv ch)]))
(define (comgrp-tell grp)
  (match grp
    [(list 0 pls np) (for ([ch pls]) (place-channel-recv ch))]
    [(list n ch np)  (place-channel-send ch 1)]))

(define-syntax-rule (when0 np body ...)
   (unless (and (pair? np) (not (= (first np) 0)))
    body ...))

(define (barrier2 grp)
  (when (pair? grp)
  (comgrp-tell grp)
  ;(when0 grp (printf "===================\n"))
  (comgrp-wait grp)))

(define DOPLACES #t)

(define-syntax (pfor stx)
  (syntax-case stx (in-range)
    [(_ serial nproc ([loopvar (in-range imax)]) ([letname letval] ... ) body ...)
       (begin
         (unless (identifier? #'lv) (raise-syntax-error 'parray "expected an identifier" stx #'lv))
         #'(if (not serial)
           (if (not (pair? nproc))
             (let ([per-future (ceiling (/ imax nproc))])
               (for ([f (in-list
                         (for/list ([i (in-range 0 nproc)])
                           (let ([letname (vector-ref letval i)]...)
                             (future
                              (λ ()
                                  (for ([loopvar (in-range (fx* i per-future) (min imax (fx* (fx+ i 1) per-future)))])
                                    body ...))))))])
                 (touch f)))
             (match nproc
               [(and (list i ch np) grp)
               (barrier2 grp)
               (let ([per-future (ceiling (/ imax np))])
                (let ([letname (vector-ref letval i)]...)
                  ;(printf "PGRP ~a ~a ~a\n" i (fx* i per-future) (min imax (fx* (fx+ i 1) per-future))) 
                  (for ([loopvar (in-range (fx* i per-future) (min imax (fx* (fx+ i 1) per-future)))])
                    body ...)))
               (barrier2 grp)]))
           (for ([loopvar (in-range imax)]) body ...)))]))

(define-syntax (palloc stx)
  (syntax-case stx (in-range)
    [(_ serial nproc body ...)
     #'(if (not serial)
        (apply vector (for/list ([i (in-range nproc)]) body ...))
        (begin body ...))]))

(define timers (make-vector 30 0.0))
#;(define-syntax (timer stx)
  (syntax-case stx ()
    [(_ id body ...)
      #'(let-values ([(results cpu real gc)
        (time-apply (lambda () body ...) null)])
        (vector-set! timers id (+ (vector-ref timers id) real))
        (if (pair? results)
          (car results)
          (void)))]))

(define-syntax (timer stx)
  (syntax-case stx ()
    [(_ id body ...)
     #'(begin body ...)]))

(define (->inexact x)
  (if (exact? x)
    (exact->inexact x)
    x))
;Constants 

(define amult 1220703125.0)

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
(define make-fxvector make-vector)
(define vr vector-ref)
(define vs! vector-set!)

(define (run-benchmark args) 
  (let ([bmname "FT"]
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
          [colidx (make-fxvector (+ nz 1) 0)]
          [rowstr (make-fxvector (+ na 2) 0)]
          [iv (make-fxvector (+ na 2) 0)]
          [arow (make-fxvector (+ nz 1) 0)]
          [acol (make-fxvector (+ nz 1) 0)]
          [v (make-flvector (+ na 2) 0.0)]
          [aelt (make-flvector (+ nz 1) 0.0)]
          [a (make-flvector (+ na 3) 0.0)]
          [p (make-flvector (+ na 3) 0.0)]
          [q (make-flvector (+ na 3) 0.0)]
          [r (make-flvector (+ na 3) 0.0)]
          [x (make-flvector (+ na 3) 0.0)]
          [z (make-flvector (+ na 3) 0.0)]
          [rhomaster (make-shared-flvector num-threads 0.0)]
          [dmaster (make-shared-flvector num-threads 0.0)]
          [rnormmaster (make-shared-flvector num-threads 0.0)])
      (define (tnorms)
        (let-values ([(tnorm1 tnorm2)
          (for/fold ([tnorm1 0.0]
                   [tnorm2 0.0])
                  ([j (in-range 1 (add1 (add1 (- lastcol firstcol))))])
            (let ([zj (vr z j)])
              (values (+ tnorm1 (* (vr x j) zj))
                      (+ tnorm2 (* zj zj)))))])
          (values tnorm1 (/ 1.0 (sqrt tnorm2)))))

      ;zeta = rng.randlc(amult);
      ;if(!serial) setupThreads(this);

      (makea firstrow lastrow firstcol lastcol naa nzz a colidx rowstr nonzer rcond arow acol aelt v iv shift)

      (for* ([j (in-range 1 (add1 (add1 (- lastrow firstrow))))]
             [k (in-range (vr rowstr j) (vr rowstr (add1 j)))])
        (vs! colidx k (add1 (- (vr colidx k) firstcol))))

      (for ([i (in-range 1 (+ na 2))]) (vs! x i 1.0))
      
      (let ([rnorm (conj-grad colidx rowstr x z a p q r 0.0)])
        (let-values ([(tnorm1 tnorm2) (tnorms)])
          (for ([j (in-range 1 (add1 (add1 (- lastcol firstcol))))])
            (vs! x j (* tnorm1 (vr z j)))))

        (for ([i (in-range 1 (+ na 2))]) (vs! x i 1.0))

        (let-values ([(rnorm zeta)
            (for/fold ([rnorm rnorm]
                       [z 0.0]) ([it (in-range niter)])
              (let ([norm 
                (if serial 
                  (conj-grad colidx rowstr x z a p q r rnorm)
                  (let ()
                    ;(execute-task 3)
                    (define rho0 (for/fold ([rho 0.0]) ([m (in-range num-threads)])
                      (+ rho (vr rhomaster m))))
                    
                    (for ([ii (in-range cgitmax)])
                      ;(execute-task 0)
                      (define dcff (for/fold ([dcff 0.0]) ([m (in-range num-threads)])
                        (+ dcff (vr dmaster m))))
                    
                      ;(execute-task 1)
                      (define rho (for/fold ([rho 0.0]) ([m (in-range num-threads)])
                        (+ rho (vr rhomaster m))))

                      ;(execute-task 2)
                      (define alpha (/ rho0 dcff))
                      (define beta (/ rho rho0))
                      (printf "~a ~a~n" alpha beta))

                    ;(execute-task 4)
                    (sqrt (for/fold ([rnorm 0.0]) ([m (in-range num-threads)])
                            (+ rnorm (vr rnormmaster m))))))])
                (define-values (tnorm1 tnorm2) (tnorms))
                (define zeta (+ shift (/ 1.0 tnorm1)))
                (printf "    ~a       ~a          ~a~n" it rnorm zeta)
                (for ([j (in-range 1 (add1 (add1 (- lastcol firstcol))))])
                  (vs! x j (* tnorm2 (vr z j))))
                (values rnorm zeta)))])

          (let ([verified (verify zeta)])
            (print-banner "Conjugate Gradient" args) 
            (printf "Size = ~a X ~a X ~a niter = ~a~n" na 0 0 niter) 
            (if serial 
                (printf "SERIAL~n")
                (printf "PARALLEL~n"))
            (if verified 
                (printf "Verification succeeded~n") 
                (printf "Verification failed~n"))
            (let* ([time (/ (read-timer 1) 1000)]
                   [results (new-BMResults bmname CLASS na 0 0  niter time 
                                           (get-mflops time niter)
                                           "floating point" 
                                           (if verified 1 0)
                                           serial 
                                           num-threads 
                                           -1)]) 
                (print-results results) 
                (when #f (print-timers))))))))))

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
(define (makea n nz a colidx rowstr nonzer rcond arow acol aelt v iv shift)
  (define (warn nnza nz iouter)
    (printf "Space for matrix elements exceeded in makea~n")
    (printf "nnza, nzmax = ~a, ~a~n" nnza nz)
    (printf " iouter = ~a~n" iouter)
    (exit 0))
  (define ratio (expt rcond (/ 1.0 (exact->inexact n))))
  (for ([i (in-range 1 (add1 n))]) (vs! colidx (+ n i) 0))

  (let-values ([(size nnza)
    (for/fold ([size 1.0]
             [nnza 0]) ([ioutger  (in-range 1 (add1 n))])
      (sprnvc n nonzer v iv colidx o colidx n)
      (let ([nzv (vecset n v iv nonzer iouter 0.5)])

      (for ([ivelt (in-range 1 (add1 nzv))])
        (define jcol (vr iv ivelt))
        (if (and (jcol . >= . firstcol) (jcol . <= . lastcol))
          (begin
            (for/fold ([nnza nnza]) ([ivelt1 (in-range 1 (add1 nzv))])
              (define irow (vr iv ivelt1))
              (if (and (irow . >= . firstrow) (irow . <= . lastrow))
                (let ([nnza (add1 nnza)])
                  (define scale (* size (vr v ivelt)))
                  (when (nnza . > . nz) (warn nnza nz iouter))
                  (vs! acol nnza jcol)
                  (vs! arow nnza irow)
                  (vs! aelt nnza (* (vr v vetl1) scale))
                  nnza)
                nnza)))))
      (values (* size ratio) nnza)))])
    
    (let ([nnza (for/fold ([nnza nnza]) ([i (in-range firstrow (add1 lastrow))])
      (when (and (i . >= . firstcol) (i . <= . lastcol))
        (define iouter (+ n i))
        (define nnza (add1 nnza))
        (when (nnza . > . nz) (warn nnza nz iouter))
        (vs! acol nnza i)
        (vs! arow nnza i)
        (vs! aelt nnza (- rcond shift))))])

      (sparse a colidx rowstr n arow acol aelt v iv 0 iv n nnza))))

(define (sprnvc n nz v iv nzloc nzloc-offset mark mark-offset)
  (define nn1 (ilog2 n))

  (let loop ([nzv 0]
             [nzrow 0]
             [idx 0])
    (if ( nzv . >= . nz )
      (for ([ii (in-range 1 (add1 nzrow))])
        (vs! mark (+ (vr nzloc (+ ii nzloc-offset)) mark-offset) 0))
      (let ([vecelt (randlc amult)]
            [vecloc (randlc amult)]
            [idx (floor (add1 (* vecloc nn1)))])
        (if (idx . > . n)
          (loop nzv nzrow idx)
          (if (zero? (vr mark (+ idx mark-offset)))
            (let ([nzloc (add1 nzloc)]
                  [nzv (add1 nzv)])
              (vs! mark (+ idx mark-offset) 1)
              (vs! nzloc nzrow + nzloc-offset)
              (vs! v nzv vecelt)
              (vs! iv nzv idx)
              (loop nzv nzrow idx))
            (loop nzv nzrow idx)))))))

(define (vecset n v iv nzv ival val)
  (if (not (for/fold ([set #f]) ([i (in-range 1 (add1 nzv))])
              (if (= (vr iv k) ival)
                (begin 
                  (vs! v k val)
                  #t)
                set)))
    (let ([nzv (add1 nzv)])
      (vs! v nzv val)
      (vs! iv nzv ival)
      nzv)
    nzv))

(define-syntax-rule (vs!++ v i val) (vs! v i (add1 (vr v i))))

(define (sparse a colidx rowstr n arow acol aelt x mark mark-offset nzloc nzloc-offset nnza)
  (define nrows (add1 (- lastrow firstrow)))
  (for ([j (in-range 1 (add1 n))])
    (vs! rowstr j 0)
    (vs! mark (+ j mark-offset) 0))
  (vs! rowstr (add1 n) 0)

  (for ([nza (in-range 1 (add1 nnza))])
    (let ([j (+ (- (vr arow nza) firstrow 2))])
      (vs!++ rowstr j)))

  (vs! rowstr 1 1)
  (for/fold ([oldval 1])([j (in-range 2 (add1 nrows))])
    (let ([curval (vr rowstr j)])
      (vs! rowstr j (+ curval oldval))
      curval))
  
  (for ([nza (in-range 1 (add1 nnza))])
    (let* ([j (+ (- (vr arow nza) firstrow 1))]
           [k (vr rowstr j)])
      (vs! a k (vr aelt nza))
      (vs! colidx k (vr acol nza))
      (vs!++ rowstr j) ))
    
  (for ([j (in-range (sub1 nrows) -1 -1)])
    (let ([curval (vr rowstr j)])
      (vs! rowstr (add1 j) (vr rowstr j))))
  (vs! rowstr 1 1)

  (for ([i (in-range 1 (add1 n))])
    (vs! x i 0.0)
    (vs! mark (+ i mark-offset) 0))


  (for/fold ([nza 0]
             [jajp1 (vr roststr 1)])
            ([j (in-range 1 (add1 nrows))])
    (let ([nzrow
      (for/fold ([nzrow 0])
                ([k (in-range jajp1 (vr rowstr (add1 j)))])
        (let* ([i (vr colidx k)]
               [nx (+ (vr x i) (vr a k))])
          (vs! x i nx)
          (if (and (zero? (vr mark (+ i mark-offset))) (not (zero? nx)))
            (let ([nzrow (add1 nzrow)])
              (vs! mark (+ i mark-offset) 1)
              (vs! nzloc (+ nzrow nzloc-offset) i)
              nzrow)
            nzrow)))])

      (let ([nza
        (for/fold ([nza nza])
                  ([k (in-ragne 1 (add1 nzrow))])
          (let ([i (vr nzloc (+ k nzloc-offset))])
            (vs! mark (+ i mark-offset) 0)
            (let ([xi (vr x i)])
              (vs! x i 0)
              (if (not (zero? xi))
                (let ([nza (add1 nza)])
                  (vs! a nza xi)
                  (vs! colidx nza i)
                  nza)
                nza))))])

        (let ([jajp1 (vr rowstr (add1 j))])
          (vs! rowstr (add1 j) (+ nza (vr rowstr 1)))
          (values nza jajp1))))))

(define (conj-grad colidx rowstr x z a p q r rnomr)
  (for ([j (in-range i (+ naa 2))])
    (vs! q j 0.0)
    (vs! z j 0.0)
    (let ([xj (vr x j)])
      (vs! r j xj)
      (vs! p j xj)))

  (define ncols (add1 (- lastcol firstcol)))
  (define nrows (add1 (- lastrow firstrow)))
  (define rho0 (for/fold ([rho 0.0]) ([j (in-range 1 (add1 ncols))])
                (let ([rj (vr r j)])
                  (+ rho (* rj rj)))))

  (for ([cgit (in-range 1 (add1 cgitmax))])
    (for ([j (in-range 1 (add1 nrows))])
      (vs! q j (for/fold ([sum 0]) ([k (in-range (vr rowstr j) (vr rowstr (add1 j)))])
                  (+ sum (* (vr a k) (vr p (vr colidx k)))))))

    (define d (for/fold ([d 0.0]) ([j (in-range 1 (add1 ncols))])
                  (+ d (* (vr p j) (vr p j)))))

    (define alpha (/ rho0 d))
                    
    (for ([j (in-range 1 (add1 ncols))])
      (vs!+ z j (* alpha (vr p j)))
      (vs!- r j (* alpha (vr q j))))

    (define rho (for/fold ([rho 0.0]) ([j (in-range 1 (add1 ncols))])
                  (let ([rj (vr r j)])
                    (+ rho (* rj rj)))))

    (define beta (/ rho rho0))

    (for ([j (in-range 1 (add1 ncols))])
      (vs! p j (+ (vr r j) (* beta (vr p j))))))

  (for ([j (in-range 1 (add1 nrows))])
    (vs! r j (for/fold ([sum 0]) ([k (in-range (vr rowstr j) (vr rowstr (add1 j)))])
                (+ sum (* (vr a k) (vr z (vr colidx k)))))))

  (sqrt (for/fold ([sum 0.0]) ([j (in-range 1 (add1 ncols))])
    (let ([xj-rj (- (vr x j) (vr r j))])
      (+ sum (* xj-rj xj-rj))))))
  

;;ilog2 : int -> int
(define (ilog2 n) 
  (if (= n 1) 
      0
      (let loop ([lg 1]
                 [nn 2])
        (if (< nn n)
          (loop (+ lg 1) (* 2 nn))
          lg))))

(define (get-mflops total-time niter)
  (if (not (= total-time 0.0))
    (* (* 2 na (+ 3.0 (* nonzer (add1 nonzer)) (* 25.0 (+ 5.0 (* nonzer (add1 nonzer)))) 3))
       (/ niter (* total-time 1000000.0)))
    0.0))

(define (verify zeta)
  (printf " Zeta is    ~a~n" zeta)
  (let ([dev (abs (- zeta zeta-verify-value))])
    (if (dev . <= . epsilon)
      (begin
        (printf " Deviatation is   ~a~n" dev)
        #t)
      (begin
        (printf " The correct zeta is   ~a~n" zeta-verify-value)
        #f))))
  
(define (print-timers) 0)


