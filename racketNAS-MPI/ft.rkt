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

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values 64 64 64 6)]
    [(#\W) (values 128 128 32 6)] 
    [(#\A) (values 256 256 128 6)]
    [(#\B) (values 512 256 256 20)]
    [(#\C) (values 512 512 512 20)]
    [else (error "Unknown class")]))

(define-syntax mk-dims-fns
  (syntax-rules ()
    [(mk-dims-fns dimsa)
      (begin
        (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
        (define (dims! x y v) (fxvector-set! dimsa (fx+ (fx* x 3) y) v))
        )]))

(define (setup id nx ny nz cnt xstart ystart zstart xend yend zend dimsa)
  (mk-dims-fns dimsa)
    (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
    (define (dims! x y v) (fxvector-set! dimsa (fx+ (fx* x 3) y) v))
  (define-values (np1 np2 layout-type)
    (cond
      [(= cnt 1) (values 1 1 '0D)]
      [(< cnt nz) (values 1 cnt '1D)]
      [else (values nz (quotient cnt nz) '2D)]))
  (case layout-type
    [('0D)
     (for ([i (in-range 1 4)])
       (dims! 0 i nx)
       (dims! 1 i ny)
       (dims! 2 i nz))]
    [('1D)
     (dims! 0 0 nx)
     (dims! 1 0 ny)
     (dims! 2 0 (quotient nz np2))
     (dims! 0 1 nx)
     (dims! 1 1 ny)
     (dims! 2 1 (quotient nz np2))
     (dims! 0 2 nz)
     (dims! 1 2 nx)
     (dims! 2 2 (quotient ny np2))
     ]
    [('2D)
     (dims! 0 0 nx)
     (dims! 1 0 (quotient ny np1))
     (dims! 2 0 (quotient nz np2))
     (dims! 0 1 ny)
     (dims! 1 1 (quotient nx np1))
     (dims! 2 1 (quotient nz np2))
     (dims! 0 2 nz)
     (dims! 1 2 (quotient nx np1))
     (dims! 2 2 (quotient ny np2))
     ])
  (define me1 (quotient id np2))
  (define me2 (modulo id np2))
  (define nz/np2 (quotient nz/np2))
  (define ny/np2 (quotient ny np2))
  (define nx/np1 (quotient nx np1))
  (define ny/np1 (quotient ny np1))
  (define (start_ xx) (fx+ 1 (fx* me2 xx)))
  (define (end_ xx) (fx* (fx+ me2 1) xx))
  (define sz2 (start_ nz/np2))
  (define ez2 (end_ nz/np2))
  (define sy2 (start_ ny/np2))
  (define ey2 (end_ ny/np2))
  (define sx1 (start_ nx/np1))
  (define ex1 (end_ nx/np1))
  (define sy1 (start_ ny/np1))
  (define ey1 (end_ ny/np1))
  (case layout-type
    [('0D)
     (for ([i (in-range 1 4)])
       (fx! xstart i 1)
       (fx! ystart i 1)
       (fx! zstart i 1)
       (fx! xend i nx)
       (fx! yend i ny)
       (fx! zend i nz))]
    [('1D)
       (fx! xstart 1 1)
       (fx! ystart 1 1)
       (fx! zstart 1 sz2)
       (fx! xend 1 nx)
       (fx! yend 1 ny)
       (fx! zend 1 ez2)
       (fx! xstart 2 1)
       (fx! ystart 2 1)
       (fx! zstart 2 sz2)
       (fx! xend 2 nx)
       (fx! yend 2 ny)
       (fx! zend 2 ez2)
       (fx! xstart 3 1)
       (fx! ystart 3 sy2)
       (fx! zstart 3 1)
       (fx! xend 3 nx)
       (fx! yend 3 ey2)
       (fx! zend 3 nz)]
    [('2D)
       (fx! xstart 1 1)
       (fx! ystart 1 sy1)
       (fx! zstart 1 sz2)
       (fx! xend 1 nx)
       (fx! yend 1 ey1)
       (fx! zend 1 ez2)
       (fx! xstart 2 sx1)
       (fx! ystart 2 1)
       (fx! zstart 2 sz2)
       (fx! xend 2 ex1)
       (fx! yend 2 ny)
       (fx! zend 2 ez2)
       (fx! xstart 3 sx1)
       (fx! ystart 3 sy2)
       (fx! zstart 3 1)
       (fx! xend 3 ex1)
       (fx! yend 3 ey2)
       (fx! zend 3 nz)]
     )

  (define-values (fftblock fftblockpad)
    (case layout-type
      [('2D)
       (let* ([fftblock (min 16 (dimsr 1 0) (dimsr 1 1) (dimsr 1 2))]
              [fftblockpad (min 18 (+ fftblock 3))])
         (values fftblock fftblockpad))]
      [else (values 16 18)]))

  (values np1 np2 me1 me2 layout-type fftblock fftblockpad))

;FIXME  idx order d3 d2 d1:
(define (evolve u0 u1 twiddle d1 d2 d3)
  (define kidxT (fx* d2 d1))
  (define jidxT d1)
  (define iidxT 1)
  (define kidxK (fx* d2 d1 2))
  (define jidxK (fx* d1 2))
  (define iidxK 2)
  (for/fold ([kidx 0] [kidx2 0]) ([k d3])
    (for/fold ([jidx kidx] [jidx2 kidx2]) ([j d2])
      (for/fold ([iidx jidx] [iidx2 jidx2]) ([i d1])
        (define iidx+1 (fx+ iidx 1))
        (define twiddle-val (flr twiddle iidx2))
        (define u00 (fl* (flr u0 iidx) twiddle-val))
        (define u01 (fl* (flr u0 iidx+1) twiddle-val))
        (fl! u0 iidx u00)
        (fl! u0 iidx+1 u01)
        (fl! u1 iidx u00)
        (fl! u1 iidx+1 u01)
        (values (fx+ iidx iidxK) (fx+ iidx2 iidxT)))
      (values (fx+ jidx jidxK) (fx+ jidx2 jidxT)))
    (values (fx+ kidx kidxK) (fx+ kidx2 kidxT))))

(define (compute-initial-conditions nx ny d1 d2 d3 u0 zstartr)
  (let ([start (randlc/a 314159265.0 
                         (ipow46 1220703125.0 (fx* 2 nx) (fx+ (fx* (fx- (zstartr 0) 1) ny) 
                                                              (fx- (zstartr 0) 1))))]
        [an (ipow46 1220703125.0 (fx* 2 nx) ny)])
    (define iidxK 2)
    (for/fold ([kidx 0] [KV start]) ([k d3])
      (for/fold ([jidx kidx] [JV KV]) ([j d2])
        (for/fold ([iidx jidx] [IV JV]) ([i d1])
          (define iidx+1 (fx+ iidx 1))
          (define-values (NV v1) (krandlc/2 IV an))
          (define-values (NNV v2) (krandlc/2 NV an))
          (fl! u0 iidx v1)
          (fl! u0 iidx+1 v2)
          (values (fx+ iidx iidxK) NNV))))))

(define (ipow46 a exp-1 exp-2)
  (if (or (= 0 exp-1) (= 0 exp-2))
      1
      (let loop ([n (quotient exp-1 2)]
                 [n% (modulo exp-1 2)]
                 [q a])
        (if (= n% 0)
            (loop (quotient n 2)
                  (modulo n 2)
                  (randlc/a q q))
            (let ([n (fx* n exp-2)])
              (let loop2 ([n (quotient n 2)]
                          [n% (modulo n 2)]
                          [q q]
                          [r 1])
                (if (> n 1)
                    (loop2 (quotient n 2)
                           (modulo n 2)
                           (if (= n% 0) (randlc/a q q) q)
                           (if (= n% 0) r (randlc/a r q)))
                    (randlc/a r q))))))))

(define (compute-indexmap layout-type twiddle d1 d2 d3 nx ny nz dimsa xstart3 ystart3 zstart3)
  (define alpha .000001)
  (define ap (fl* -4.0 alpha pi pi))
  (define (XX xx Xstart3 NX) 
    (define NX/2 (fx/ NX 2))
    (fx- (modulo (fx+ xx Xstart3 -2 NX/2) NX) NX/2))
  (define-syntax-rule (worker RDIM1 RDIM2 RDIM3 II JJ KK RI RJ RK)
    (for ([II RDIM1])
      (define ii (XX II xstart3 nx))
      (define ii2 (fx* ii ii))
      (for ([JJ RDIM2])
        (define jj (XX JJ ystart3 ny))
        (define ij2 (fx+ (fx* jj jj) ii2))
        (for ([KK RDIM3])
          (define kk (XX KK zstart3 nz))
          (define IDX (fx+ (fx* RI d2 d3) (fx* RJ d3) RK))
          (fl! twiddle IDX (exp (fl* ap (exact->inexact (fx+ (fx* kk kk) ij2)))))))))

  (match-define (vector d00 d01 d02 d10 d11 d12 d20 d21 d22) dimsa)
  (case layout-type
    [(D0) (worker d02 d12 d22 i j k i j k)]
    [(D1 D2) (worker d12 d22 d02 i j k k i j)]
    [else (error (format "unknown layout type ~a" layout-type))])
  )

(define-values (T_TOTAL T_SETUP T_FFT T_EVOLVE T_CHECKSUM T_FFTLOW T_FFTCOPY
                T_TRANSPOSE 
                T_TRANSXZLOC T_TRANSXZGLO T_TRANSXZFIN
                T_TRANSXYLOC T_TRANSXYGLO T_TRANSXYFIN
                T_SYNCH T_INIT T_TOTALCOMP T_TOTALCOMM T_MAX
                )
               (values 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18))

(define (print-timers id comm np)
  (define t1 (for/flvector #:length T_MAX ([i T_MAX])
    (timer-read i)))

  (fl! t1 T_TOTALCOMM (+ (flr t1 T_TRANSXZGLO)
                         (flr t1 T_TRANSXYGLO)
                         (flr t1 T_SYNCH)))

  (fl! t1 T_TOTALCOMP (- (flr t1 T_TOTAL)
                         (flr t1 T_TOTALCOMM)))

  (define tsum (rmpi-reduce comm 0 + t1))
  (define tming (rmpi-reduce comm 0 min t1))
  (define tmaxg (rmpi-reduce comm 0 max t1))

  (when (= id 0)
    (define descrs '(
    "          total " 
    "          setup " 
    "            fft " 
    "         evolve " 
    "       checksum " 
    "         fftlow " 
    "        fftcopy " 
    "      transpose " 
    " transpose1_loc " 
    " transpose1_glo " 
    " transpose1_fin " 
    " transpose2_loc " 
    " transpose2_glo " 
    " transpose2_fin " 
    "           sync "
    "           init "
    "        totcomp "
    "        totcomm "))

    (printf "nprocs ~a minimum maximum average\n" np)
    (for ([i (in-naturals)]
          [descr descrs]
          [ts (in-vector tsum)]
          [tm (in-vector tming)]
          [tma (in-vector tmaxg)])
      (printf "timer ~a (~a) : ~a ~a ~a\n"
        (~r #:min-wideth 2 i)
        descr
        (~r #:precision 4 #:min-width 14 ts)
        (~r #:precision 4 #:min-width 14 tm)
        (~r #:precision 4 #:min-width 14 tma)))))

(define (fft dir x1 x2 scratch dimsa layout-type)
  (match-define (vector d00 d01 d02 d10 d11 d12 d20 d21 d22) dimsa)
  (case layout-type
    [(D0)
      (cond [(= dir 1)
             (cffts1 1 d00 d10 d20 x1 x1 scratch)
             (cffts2 1 d01 d11 d21 x1 x1 scratch)
             (cffts3 1 d02 d12 d22 x1 x2 scratch)]
            [else
             (cffts3 -1 d02 d12 d22 x1 x1 scratch)
             (cffts2 -1 d01 d11 d21 x1 x1 scratch)
             (cffts1 -1 d00 d10 d20 x1 x2 scratch)])]
    [(D1)
      (cond [(= dir 1)
             (cffts1 1 d00 d10 d20 x1 x1 scratch)
             (cffts2 1 d01 d11 d21 x1 x1 scratch)
             (transpose-xy-z 2 3 x1 x2)
             (cffts1 1 d02 d12 d22 x2 x2 scratch)]
            [else
             (cffts1 -1 d02 d12 d22 x1 x1 scratch)
             (transpose-x-yz 3 2 x1 x2)
             (cffts1 -1 d00 d10 d20 x2 x2 scratch)
             (cffts2 -1 d01 d11 d21 x2 x2 scratch)])]
    [(D3)
      (cond [(= dir 1)
             (cffts1 1 d00 d10 d20 x1 x1 scratch)
             (transpose-x-y 1 2 x1 x2)
             (cffts1 1 d01 d11 d21 x2 x2 scratch)
             (transpose-x-z 2 3 x2 x1)
             (cffts1 1 d02 d12 d22 x1 x2 scratch)]
            [else
             (cffts1 -1 d02 d12 d22 x1 x1 scratch)
             (transpose-x-z 3 2 x1 x2)
             (cffts1 -1 d01 d11 d21 x2 x2 scratch)
             (transpose-x-y 2 1 x2 x1)
             (cffts1 -1 d00 d10 d20 x1 x2 scratch)])]
   ))

(define-syntax-rule (MK-cffts NAME i j k d1 d2 d3 fftblock DX K D_32 JJII D_21 J I1 I2 I+ J+ X1 X2)
  (define (NAME is d1 d2 d3 fftblock x xout y)
    (define logdX (ilog2 DX))

    (for ([K D_32])
      (for ([JJII (in-range 0 D_21 fftblock)])
        (with-timer T_FFTCOPY
          (for ([J I1])
            (for ([i I2])
              (fl! y (fx+ (fx* X1 DX 2 2) (fx* X2 2 2) (fx* 1 2))
                   (flr x (fx+ (fx* I+ d2 d3 2) (fx* J+ d3 2) (fx* k 2))))
              (fl! y (fx+ (fx* X1 DX 2 2) (fx* X2 2 2) (fx* 1 2) 1)
                   (flr x (fx+ (fx* I+ d2 d3 2) (fx* J+ d3 2) (fx* k 2) 1)))
               )))

        (with-timer T_FFTLOW
          (cfftz is logdX DX y (flr y (fx+ (fx* 2 2)))))

        (with-timer T_FFTCOPY
          (for ([J I1])
            (for ([i I2])
              (fl! xout (fx+ (fx* I+ d2 d3 2) (fx* J+ d3 2) (fx* k 2))
                   (flr y (fx+ (fx* X1 DX 2 2) (fx* X2 2 2) (fx* 1 2))))
              (fl! xout (fx+ (fx* I+ d2 d3 2) (fx* J+ d3 2) (fx* k 2) 1)
                   (flr y (fx+ (fx* X1 DX 2 2) (fx* X2 2 2) (fx* 1 2) 1)))
               )))
    
    ))))
;         NAME   i j k DX K D_32 JJ D_21 J I1       I2         I+       J+ X1 X2)
;                                II
;                      DX K D_32 JJ D_21 J I1       I2       JII+       J+    LOGD  DX   I1             I2 I+ J+    X1 X2
(MK-cffts cffts1 i j k d1 d2 d3 fftblock d1 k d3 jj d2   j fftblock d1       i        (+ j jj) j i)
(MK-cffts cffts2 i j k d1 d2 d3 fftblock d2 k d3 ii d1   j d2       fftblock (+ i ii) j        i j)
(MK-cffts cffts3 i j k d1 d2 d3 fftblock d3 j d2 ii d2   k d3       fftblock (+ i ii) j        i k)

; DX K D_32 JJ D_21 J I1         I2       JI     I+ J+    LOGD  DX   I1             I2 I+ J+    X1 X2
; d1 k d3   jj d2   j fftblock i       d1 j  i 1 i j+jj k logd1 d1 j fftblock       d1 i j+jj k j i
; d2 k d3   ii d1   j d2       i fftblock i  j 1 i+ii j k logd2 d2 j d2       fftblock i+ii j k i j
; d3 j d2   ii d1   k d3       i fftblock i  k 1 i+ii j k logd3 d3 j d3       fftblock i+ii j k i k
;
(define (->inexact x)
  (if (exact? x) (exact->inexact x) x))

(define (fft-init n u)
  (define REAL 0)
  (define IMAG 1)
  (let ([nu n]
        [m (ilog2 n)] 
        [eps 1.0E-16] 
        [ku 1] 
        [ln 1]) 
    (fl! u 0 (->inexact m))
    (for/fold ([ku 2] [ln 1]) ([j (in-range 1 (+ m 1))]) 
      (let ([t (/ pi ln)]) 
        (for ([i (in-range ln)])
          (let ([ti (* i t)] 
                [idx (* (+ i ku) 2)]) 
            (fl! u (+ REAL idx) (->inexact (cos ti)))
            (fl! u (+ IMAG idx) (->inexact (sin ti)))
            (when (< (abs (fl! u (+ REAL idx))) eps) 
              (error (format "smaller thatn eps ~a ~a\n" (fl! u (+ REAL idx)) eps))
              (fl! u (+ REAL idx) 0.0))
            (when (< (abs (fl! u (+ IMAG idx))) eps) 
              (error (format "smaller thatn eps ~a ~a\n" (fl! u (+ IMAG idx)) eps))
              (fl! u (+ IMAG idx) 0.0))
)))
      (values (+ ku ln) (* 2 ln)))))

(define (cfftz is m n u x y fftblock fftblockpad)
  (define mx (flr u 0))
  (when (or
          (and (not (= is 1))
               (not (= is -1)))
          (< m 1)
          (> m mx))
    (printf "CFFTZ: Either U has not been initialized, or else one of the input parmeters is invalid ~a ~a ~a\n" is m mx))

  (let/ec k
    (for ([l (in-range 0 m 2)])
      (fftz2 is l m n fftblock fftblockpad u x y)
      (when (= l m ) (k))
      (fftz2 is (+ l 1) m n fftblock fftblockpad u x y)))

  (for* ([j n]
         [i fftblock])
    (define idx (fx+ (fx* i n 2) (fx* j 2)))
    (fl! x idx (flr y idx))
    (fl! x (fx+ idx 1) (flr y (fx+ idx 1)))))

(define (fftz2 is l m n ny ny1 u x y)
  (define nl (/ n 2))
  (define lk (expt 2 (fx- l 1)))
  (define li (expt 2 (fx- m l)))
  (define lj (fx* 2 lk))
  (define ku (fx+ li 1))

  (for ([i li])
    (define i11 (fx+ (fx* i lk ) 1))
    (define i12 (fx+ i11 nl))
    (define i21 (fx+ (fx* i lj ) 1))
    (define i22 (fx+ i21 lk))
    (define u10 (flr u (fx* (fx+ ku i) 2)))

    (define u11 (let ([v (flr u (fx+ (fx* (fx+ ku i) 2) 1))])
                  (if (> is 1) (fx- 0 v) v)))

    (for* ([k lk]
           [j ny])
      (define x110 (flr x (fx+ (fx* j 2) (fx* (fx+ i11 k) 2))))
      (define x111 (flr x (fx+ (fx* j 2) (fx* (fx+ i11 k) 2)) 1))
      (define x210 (flr x (fx+ (fx* j 2) (fx* (fx+ i12 k) 2))))
      (define x211(flr x (fx+ (fx* j 2) (fx* (fx+ i12 k) 2)) 1))
      (fl! y (fx+ (fx* j 2) (fx* (fx+ i21 k) 2)) (fl+ x110 x210))
      (fl! y (fx+ (fx* j 2) (fx* (fx+ i21 k) 2) 1) (fx+ x111 x211))
      (fl! y (fx+ (fx* j 2) (fx* (fx+ i21 k) 2)) (fx* u10 (fx- x110 x210)))
      (fl! y (fx+ (fx* j 2) (fx* (fx+ i21 k) 2) 1) (fx* u11 (fx- x111 x211)))
   )))

(define (ilog2 n) 
  (if (= n 1) 
      0
      (let loop ([lg 1]
                 [nn 2])
        (if (fx< nn n)
          (loop (fx+ lg 1) (fx* 2 nn))
          lg))))

(define (transpose-x-yz l1 l2 xin xout dimsa)
  (mk-dims-fns dimsa)
    (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
  (define d0_l1 (dimsr 0 l1))
  (define d1_l1*d2_l1 (fx* (dimsr 1 l1) (dimsr 2 l1)))
  (transpose2-local d0_l1 d1_l1*d2_l1 xin xout)
  (transpose2-global xout xin)
  (transpose2-finish d0_l1 d1_l1*d2_l1 xin xout))

(define (transpose-xy-z l1 l2 xin xout dimsa)
    (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
  (define d0_l1*d1_l1 (fx* (dimsr 0 l1) (dimsr 1 l1)))
  (define d2_l1 (dimsr 2 l1))
  (transpose2-local d0_l1*d1_l1 d2_l1 xin xout)
  (transpose2-global xout xin)
  (transpose2-finish d0_l1*d1_l1 d2_l1 xin xout))

(define (transpose2-local n1 n2 xin xout transblock transblockpad)
  (define z (make-flvector (fx* transblockpad transblock 2)))
  (with-timer T_TRANSXZLOC
    (if (or (< n1 transblock) (< n2 transblock))
        (if (>= n1 n2)
            (for* ([j n2]
                   [i n1])
              (fl! xout (fx+ (fx* j n1 2) (fx* i 2))
                   (flr xin (fx+ (fx* i n2 2) (fx* j 2)))))
            (for* ([i n1]
                   [j n2])
              (fl! xout (fx+ (fx* j n1 2) (fx* i 2))
                   (flr xin (fx+ (fx* i n2 2) (fx* j 2))))))
        (for* ([j (in-range 0 n2 transblock)]
               [i (in-range 0 n1 transblock)])
          (for* ([jj transblock]
                 [ii transblock])
            (fl! z (fx+ (fx* jj transblock 2) (fx* ii 2))
                 (flr xin (fx+ (fx* (fx+ i ii) n2 2) (fx* (fx+ j jj) 2))))
            (fl! z (fx+ (fx* jj transblock 2) (fx* ii 2) 1)
                 (flr xin (fx+ (fx* (fx+ i ii) n2 2) (fx* (fx+ j jj) 2) 1))))
          (for* ([ii transblock]
                 [jj transblock])
            (fl! xout (fx+ (fx* (fx+ j jj) n1 2) (fx* (fx+ i ii) 2))
                 (flr z (fx+ (fx* jj transblock 2) (fx* ii 2))))
            (fl! xout (fx+ (fx* (fx+ j jj) n2 2) (fx* (fx+ i ii) 2) 1)
                 (fl! z (fx+ (fx* jj transblock 2) (fx* ii 2) 1))))
              ))))

(define (transpose2-global xin xout commslice1)
  (with-timer T_TRANSXZGLO
    (rmpi-alltoall commslice1 xin)))

(define (transpose2-finish n1 n2 xin xout np np2) 
  (with-timer T_TRANSXZFIN
    (for ([p np2])
      (define ioff (* p n2))
      (for* ([j (fx/ n1 np2)]
             [i n2])
        (fl! xout (fx+ (fx* (fx+ i ioff) (fx/ n1 np2) 2) (fx* j 2))
             (flr xin (fx+ (fx* i (fx* np2 2)) (fx* j np 2) (fx* p 2))))
        (fl! xout (fx+ (fx* (fx+ i ioff) (fx/ n1 np2) 2) (fx* j 2) 1)
             (flr xin (fx+ (fx* i (fx* np2 2)) (fx* j np 2) (fx* p 2) 1)))
        ))))

(define (transpose-x-z l1 l2 xin xout dimsa)
    (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
  (define d0-l1 (dimsr 0 l1))
  (define d1-l1 (dimsr 1 l1))
  (define d2-l1 (dimsr 2 l1))
  (define d0-l2 (dimsr 0 l2))
  (define d1-l2 (dimsr 1 l2))
  (define d2-l2 (dimsr 2 l2))
  (transpose-x-z-local d0-l1 d1-l1 d2-l1 xin xout)
  (transpose-x-z-global d0-l1 d1-l1 d2-l1 xin xout)
  (transpose-x-z-finish d0-l2 d1-l2 d2-l2 xin xout)
  )

(define (transpose-x-z-local d1 d2 d3 xin xout transblock transblockpad maxdim)
  (define buf (make-flvector (fx* transblockpad maxdim 2)))
  (with-timer T_TRANSXZLOC
    (cond
      [(or (< d1 32) (= d3 1))
       (for* ([j d2]
              [k d3]
              [i d1])
         (fl! xout (fx+ (fx* k d2 d1 2) (fx* j d1 2) (fx* i 2)) 
              (flr xin (fx+ (fx* i d2 d3 2) (fx* j d3 2) (fx* k 2))))
         (fl! xout (fx+ (fx* k d2 d1 2) (fx* j d1 2) (fx* i 2) 1) 
              (flr xin (fx+ (fx* i d2 d3 2) (fx* j d3 2) (fx* k 2) 1)))
         )]
      [else
        (define block3 (if (> d3 transblock) transblock d3))
        (define block1 (if (> (fx* d1 block3) (fx* transblock transblock))
                           (fl/ (fx* transblock transblock) block3)
                           d1))
        (for* ([j d2]
               [kk (in-range 0 (fx- d3 block3) block3)]
               [ii (in-range 0 (fx- d1 block1) block1)])
              (for ([k block3])
                (define k1 (fx+ k kk))
                (for ([i block1])
                  (define idx1 (fx+ (fx* k maxdim 2) (fx* i 2)))
                  (define idx2 (fx+ (fx* (fx+ i ii) d2 d3 2) (fx* j d3 2) (fx* k1 2)))
                  (fl! buf idx1 (flr xin idx2))
                  (fl! buf (fx+ idx1 1) (flr xin (fx+ idx2 1)))))
              (for ([i block1])
                (define i1 (fx+ i ii))
                (for ([k block3])
                  (define idx1 (fx+ (fx* k maxdim 2) (fx* i 2)))
                  (define idx2 (fx+ (fx* (fx+ k kk) d2 d1 2) (fx* j d1 2) (fx* i1 2)))
                  (fl! xout idx2 (flr buf idx1))
                  (fl! xout (fx+ idx2 1) (flr buf (fx+ idx1 1)))))
              )])))

(define (transpose-x-z-global xin xout commslice1)
  (with-timer T_TRANSXZGLO
    (rmpi-alltoall commslice1 xin)))

(define (transpose-x-z-finish n1 n2 xin xout np2 d1 d2 d3)
  (with-timer T_TRANSXZFIN
    (for ([p np2])
      (define ioff (fx/ (fx* p d1) np2))
      (for* ([k d3]
             [j d2]
             [i (fx/ d1 np2)])
        (define idx1 (fx+ (fx* i d2 d3 np2 2) (fx* j d3 np2 2) (fx* k np2 2) (fx* p 2)))
        (define idx2 (fx+ (fx* (fx+ i ioff) d2 d3 2) (fx* j d3 2) (fx* k 2)))
        (fl! xout idx2 (flr xin idx1))
        (fl! xout (fx+ idx2 1) (flr xin (fx+ idx1 1)))))))

(define (transpose-x-y l1 l2 xin xout dimsa)
    (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
  (define d0-l1 (dimsr 0 l1))
  (define d1-l1 (dimsr 1 l1))
  (define d2-l1 (dimsr 2 l1))
  (define d0-l2 (dimsr 0 l2))
  (define d1-l2 (dimsr 1 l2))
  (define d2-l2 (dimsr 2 l2))
  (transpose-x-y-local d0-l1 d1-l1 d2-l1 xin xout)
  (transpose-x-y-global d0-l1 d1-l1 d2-l1 xin xout)
  (transpose-x-y-finish d0-l2 d1-l2 d2-l2 xin xout)
  )

(define (transpose-x-y-local d1 d2 d3 xin xout)
  (with-timer T_TRANSXYLOC
    (for* ([k d3]
           [i d1]
           [j d2])
      (define idx1 (fx+ (fx* i d2 d3 2) (fx* j d3 2) (fx* k 2)))
      (define idx2 (fx+ (fx* j d3 d1 2) (fx* k d1 2) (fx* i 2)))
      (fl! xout idx2 (flr xin idx1))
      (fl! xout (fx+ idx2 1) (flr xin (fx+ idx1 1)))
      )))

(define (transpose-x-y-global xin xout commslice2)
  (with-timer T_TRANSXYGLO
    (rmpi-alltoall commslice2 xin)))

(define (transpose-x-y-finish n1 n2 xin xout np1 d1 d2 d3)
  (with-timer T_TRANSXYFIN
    (for ([p np1])
      (define ioff (fx/ (fx* p d1) np1))
      (for* ([k d3]
             [j d2]
             [i (fx/ d1 np1)])
        (define idx1 (fx+ (fx* i d2 d2 np1 2) (fx* k d2 np1 2) (fx* j np1 2) (fx* p 2)))
        (define idx2 (fx+ (fx* (fx+ i ioff) d2 d3 2) (fx* j d3 2) (fx* k 2)))
        (fl! xout idx2 (flr xin idx1))
        (fl! xout (fx+ idx2 1) (flr xin (fx+ idx1 1)))))))


(define (checksum id i sum u1 d1 d2 d3 nx ny nz xstart1 xend1 ystart1 yend1 zstart1 zend1 u ntotal-f comm)
  (define isize3 2)
  [define jsize3 (* isize3 (+ d1 1))] 
  [define ksize3 (* jsize3 d2)]
  [define d123 (* d1 (* d2 d3))]
  (let-values ([(csumr csumi)
    (for/fold ([csumr 0.0]
               [csumi 0.0])
              ([i (in-range 1 1025)]) 
      (define q (fxmodulo i nx))
      (cond
        [(and (>= q xstart1) (<= q xend1))
          (define r (fxmodulo (fx* 3 i) ny))
          (cond
            [(and (>= r ystart1) (<= r yend1))
              (define s (fxmodulo (fx* 5 i) nz))
              (cond
                [(and (>= s zstart1) (<= s zend1))
                 (define idx (fx+ (fx+ (fx* (fx- q xstart1) isize3) 
                                       (fx* (fx- r ystart1) jsize3))
                                       (fx* (fx- s zstart1) ksize3)))
                 (values
                   (fl+ csumr (flr u idx))
                   (fl+ csumi (flr u (fx+ idx 1))))]
                [else (values csumr csumi)])]
            [else (values csumr csumi)])]
        [else (values csumr csumi)]))])
    
    (define chk (flvector
                  (fl/ csumr ntotal-f)
                  (fl/ csumi ntotal-f)))

    (define allchk
      (with-timer T_SYNCH
        (rmpi-reduce comm 0 + chk)))

    (match-define (vector chkr chki) allchk)

    (when (= id 0)
      (printf " T = ~a Checksum = ~a ~a\n" i chkr chki))
    (when (> i 0)
      (fl! sum (fx* i 2) chkr)
      (fl! sum (fx+ (fx* i 2) 1) chki))))

(define (synchup comm)
  (with-timer T_SYNCH
    (rmpi-barrier comm)))

(define (verify class niter-default cksum)
  (define cexpd
    (case class
       ;Class S reference values
    [(#\S) 
         #(
         554.6087004964 484.5363331978;IMAG 0
         554.6385409189 486.5304269511;IMAG 1  
         554.6148406171 488.3910722336;IMAG 2
         554.5423607415 490.1273169046;IMAG 3
         554.4255039624 491.7475857993;IMAG 4
         554.2683411902 493.2597244941;IMAG 5
        )]

    [(#\W)
          #( 
            567.3612178944 529.3246849175;IMAG 0
            563.1436885271 528.2149986629;IMAG 1
            559.4024089970 527.0996558037;IMAG 2
            556.0698047020 526.0027904925;IMAG 3
            553.0898991250 524.9400845633;IMAG 4
            550.4159734538 523.9212247086;IMAG 5
          )]

    [(#\A) 
          #( 
            504.6735008193 511.4047905510;IMAG 0
            505.9412319734 509.8809666433;IMAG 1
            506.9376896287 509.8144042213;IMAG 2
            507.7892868474 510.1336130759;IMAG 3
            508.5233095391 510.4914655194;IMAG 4
            509.1487099959 510.7917842803;IMAG 5
          )]
    [(#\B) 
     ;Class B reference values 
          #(
            517.7643571579 507.7803458597;IMAG 0
            515.4521291263 508.8249431599;IMAG 1
            514.6409228649 509.6208912659;IMAG 2
            514.2378756213 510.1023387619;IMAG 3
            513.9626667737 510.3976610617;IMAG 4
            513.7423460082 510.5948019802;IMAG 5
            513.5547056878 510.7404165783;IMAG 6
            513.3910925466 510.8576573661;IMAG 7
            513.2470705390 510.9577278523;IMAG 8
            513.1197729984 511.0460304483;IMAG 9
            513.0070319283 511.1252433800;IMAG 10
            512.9070537032 511.1968077718;IMAG 11
            512.8182883502 511.2616233064;IMAG 12
            512.7393733383 511.3203605551;IMAG 13
            512.6691062020 511.3735928093;IMAG 14
            512.6064276004 511.4218460548;IMAG 15
            512.5504076570 511.4656139760;IMAG 16
            512.5002331720 511.5053595966;IMAG 17
            512.4551951846 511.5415130407;IMAG 18
            512.4146770029 511.5744692211;IMAG 19
          )]
    [(#\C) 
     ;Class C reference values 
          #(
            519.5078707457  514.9019699238;IMAG 0   
            515.5422171134  512.7578201997;IMAG 2
            514.4678022222  512.2251847514;IMAG 3
            514.0150594328  512.1090289018;IMAG 4
            513.7550426810  512.1143685824;IMAG 5
            513.5811056728  512.1496764568;IMAG 6
            513.4569343165  512.1870921893;IMAG 7
            513.3651975661  512.2193250322;IMAG 8
            513.2955192805  512.2454735794;IMAG 9
            513.2410471738  512.2663649603;IMAG 10
            513.1971141679  512.2830879827;IMAG 11
            513.1605205716  512.2965869718;IMAG 12
            513.1290734194  512.3075927445;IMAG 13
            513.1012720314  512.3166486553;IMAG 14
            513.0760908195  512.3241541685;IMAG 15
            513.0528295923  512.3304037599;IMAG 16
            513.0310107773  512.3356167976;IMAG 17
            513.0103090133  512.3399592211;IMAG 18
            512.9905029333  512.3435588985;IMAG 19
            512.9714421109  512.3465164008;IMAG 20
          )]
    [(#\D)
        #(
         512.2230065252D 511.8534037109D
         512.0463975765D 511.7061181082D
         511.9865766760D 511.7096364601D
         511.9518799488D 511.7373863950D
         511.9269088223D 511.7680347632D
         511.9082416858D 511.7967875532D
         511.8943814638D 511.8225281841D
         511.8842385057D 511.8451629348D
         511.8769435632D 511.8649119387D
         511.8718203448D 511.8820803844D
         511.8683569061D 511.8969781011D
         511.8661708593D 511.9098918835D
         511.8649768950D 511.9210777066D
         511.8645605626D 511.9307604484D
         511.8647586618D 511.9391362671D
         511.8654451572D 511.9463757241D
         511.8665212451D 511.9526269238D
         511.8679083821D 511.9580184108D
         511.8695433664D 511.9626617538D
         511.8713748264D 511.9666538138D
         511.8733606701D 511.9700787219D
         511.8754661974D 511.9730095953D
         511.8776626738D 511.9755100241D
         511.8799262314D 511.9776353561D
         511.8822370068D 511.9794338060D
)]
    [(#\E)
        #(
         512.1601045346 511.7395998266D
         512.0905403678 511.8614716182D
         512.0623229306 511.9074203747D
         512.0438418997 511.9345900733D
         512.0311521872 511.9551325550D
         512.0226088809 511.9720179919D
         512.0169296534 511.9861371665D
         512.0131225172 511.9979364402D
         512.0104767108 512.0077674092D
         512.0085127969 512.0159443121D
         512.0069224127 512.0227453670D
         512.0055158164 512.0284096041D
         512.0041820159 512.0331373793D
         512.0028605402 512.0370938679D
         512.0015223011 512.0404138831D
         512.0001570022 512.0432068837D
         511.9987650555 512.0455615860D
         511.9973525091 512.0475499442D
         511.9959279472 512.0492304629D
         511.9945006558 512.0506508902D
         511.9930795911 512.0518503782D
         511.9916728462 512.0528612016D
         511.9902874185 512.0537101195D
         511.9889291565 512.0544194514D
         511.9876028049 512.0550079284D
)]

          ))
  (define epsilon 1.0E-12)
  (if (niter-default . <= . 0) 
    -1
    (for/fold ([verified #t]) ([it (in-range niter-default)]) 
      (let* ([it2 (* it 2)]
             [rit2 (+ 0 it2)]
             [iit2 (+ 1 it2)]
             [cexpdr (vector-ref cexpd rit2)]
             [cexpdi (vector-ref cexpd iit2)]
             [csumr (/ (- (flr cksum rit2) cexpdr) cexpdr)] 
             [csumi (/ (- (flr cksum iit2) cexpdi) cexpdi)])
        (if (or (<= (abs csumr) epsilon)
                (<= (abs csumi) epsilon)) 
            (and verified #t)
            (begin 
              (printf "Verification failure: ~n") 
              (printf "epsilon: ~a~n" epsilon)
              (printf "csumr: ~a~n" csumr) 
              (printf "csumi: ~a~n" csumi)
              (and verified #f)))))))

(define/provide (ft-place ch)
  (define-values (comm args tc) (rmpi-init ch))
  (match-define (list class bmname num-threads) args)
  (define id (rmpi-id comm))
  (define cnt (rmpi-cnt comm))
  (define-values (nx ny nz niter) (get-class-size class))
  (define dimsa (make-fxvector 9))
  (define xstart (make-fxvector 3))
  (define ystart (make-fxvector 3))
  (define zstart (make-fxvector 3))
  (define xend (make-fxvector 3))
  (define yend (make-fxvector 3))
  (define zend (make-fxvector 3))
  (define dims (make-fxvector (fx* 3 3)))
  (define ntdivnp (fx* (fx/ (fx* nx ny) cnt) nz))
  (define u (make-flvector (fx* 2 nx)))
  (define u0 (make-flvector (fx* 2 ntdivnp)))
  (define u1 (make-flvector (fx* 2 ntdivnp)))
  (define u2 (make-flvector (fx* 2 ntdivnp)))
  (define twiddle (make-flvector ntdivnp))

  (timer-start T_INIT)
  (define (dimsr x y)   (fxvector-ref  dimsa (fx+ (fx* x 3) y)))
  (define ntotal-f (* nx ny nz 1.0))
  (define-values (np1 np2 me1 me2 layout-type fftblock fftblockpad) (setup xstart ystart zstart xend yend zend dimsa))
  (compute-indexmap twiddle (dimsr 0 2) (dimsr 1 2) (dimsr 2 2))
  (compute-initial-conditions u1 (dimsr 0 0) (dimsr 1 0) (dimsr 2 0))
  (fft-init (dimsr 0 0))
  (fft 1 u1 u0)
  (timer-stop T_INIT)

  (rmpi-finish comm tc)
  
  (for ([i T_MAX]) (timer-clear i))

  (rmpi-barrier comm)

  (with-timer! T_TOTAL
    (with-timer T_SETUP
      (compute-indexmap twiddle (dimsr 0 2) (dimsr 1 2) (dimsr 2 2))
      (compute-initial-conditions u1 (dimsr 0 0) (dimsr 1 0) (dimsr 2 0))
      (fft-init (dimsr 0 0)))
    (with-timer T_FFT
      (fft 1 u1 u0))

    (for ([i niter])
      (with-timer T_EVOLVE
        (evolve u0 u1 twiddle (dimsr 0 0) (dimsr 1 0) (dims 2 0)))
      (with-timer T_FFT
        (fft -1 u1 u2))
      (with-timer T_CHECKSUM
        (checksum i u2 (dimsr 0 0) (dims 1 0) (dims 2 0)))))

  (define verified (verify nx ny nz niter class))

  (define total-time (timer-read T_TOTAL))
  (define mflops
    (cond 
      [(= total-time 0.0) 0.0]
      [else
        (/
          (+
            (* .000001 ntotal-f (+ 14.8157 (* 7.19641 (log ntotal-f))))
            (* (+ 5.23518 (* 7.21113 (log ntotal-f))) niter))
          total-time)]))

  (when (= id 0)
    (print-results
      "FT"
      class
      nx
      ny
      nz
      niter
      cnt
      cnt
      total-time
      mflops
      "          floating point"
      verified
      "3.3"
      ""
      ""
      ""
      ""
      ""
      ""
      ""
      ""
      )
    )

  (print-timers)

  (rmpi-finish tc)
  )

(define (main . argv) 
  (let* ([args (parse-cmd-line-args argv "Fourier Transform")]
         [class (BMArgs-class args)]
         [serial (BMArgs-serial args)]
         [num-threads (BMArgs-num-threads args)]
         [bmname "FT"])

    (print-banner "Fourier Transform" args) 

    (define place-args (list class bmname num-threads))

    (rmpi-launch
      (rmpi-build-default-config
        #:mpi-module (quote-module-path)
        #:mpi-func 'ft-place
        #:mpi-args place-args)
      (rmpi-make-localhost-config num-threads 6341 'ft))))

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

(module+ test (void))
