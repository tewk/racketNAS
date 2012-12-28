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

(define (ispow2 i)
  (cond
    [(< i 0) #f]
    [(= i 0) 1]
    [else
      (let loop ([i i]
                 [x 1])
        (cond 
          [(zero? i) x]
          [else (loop (- i 1) (arithmetic-shift x 1))]))]))

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values  32  32  32  4 5 )] 
    [(#\W) (values  64  64  64 40 6 )]
    [(#\A) (values 256 256 256  4 8 )]
    [(#\B) (values 256 256 256 20 8 )] 
    [(#\C) (values 512 512 512 20 9 )]
    [(#\D) (values 1024 1024 1024 50 10 )]
    [(#\E) (values 2048 2048 2048 50 11 )]
    [else (error "Unknown class")]))

(define/provide (mg-place ch)
  (define-values (comm args tc) (rmpi-init ch))
  (match-define (list class bmname num-threads) args)
  (define id (rmpi-id comm))
  (define nprocs 256)
  ;(define nprocs (rmpi-cnt comm))
  (define maxlevel 11)

  (define-values  (T_BENCH T_INIT T_PSINV T_RESID T_RPRJ3 T_INTERP T_NORM2U3 T_COMM3 T_RCOMM T_TOTCOMP T_TOTCOMM T_LAST)
                  (values 0 1 2 3 4 5 6 7 8 9 10 11))
  (define t-recs '( total init psinv resid rprj3 interp norm2u3 comm3 rcom totcomp totcomm))

  (define-values (nx_default ny_default nz_default nit_default lt_default) (get-class-size class))
  (for ([i T_LAST]) (timer-clear i))
  (rmpi-barrier comm)
  (with-timer T_INIT

  (define nx (make-vector maxlevel 0))
  (define ny (make-vector maxlevel 0))
  (define nz (make-vector maxlevel 0))
  (define ir (make-vector maxlevel 0))
  (define m1 (make-vector maxlevel 0))
  (define m2 (make-vector maxlevel 0))
  (define m3 (make-vector maxlevel 0))
  (define lt lt_default)
  (define log2-size (ilog2 nx_default))
  (define log2-nprocs (ilog2 nprocs))
  (define lm (fx- log2-size (fx/ log2-nprocs 3)))
  (define ndim1 lm)
  (define ndim2 (fx- log2-size (fx/ (fx+ log2-nprocs 2) 3)))
  (define ndim3 (fx- log2-size (fx/ (fx+ log2-nprocs 1) 3)))

  (define lb 1)
  (define k lt)

  (printf "LM ndim1 ndim2 ndim3 ~a ~a ~a ~a\n" lm ndim1 ndim2 ndim3)

  (setup )

  #|
  ;(define lt1 (fx-- lt)) 
  (define nit nit_default)
  (define nm (+ (arithmetic-shift 1 lm) 2))
  (define nv (* (+ 2 (arithmetic-shift 1 ndim1)) 
         (+ 2 (arithmetic-shift 1 ndim2)) 
         (+ 2 (arithmetic-shift 1 ndim3))))
  (define nr (floor (/ (* 8 (+ nv (* nm nm) (* 5 nm) (* 7 lm))) 7)))
  (define nm2 (* 2 nm nm))
  (define r (make-flvector nr 0.0))
  (define v (make-flvector nv 0.0))
  (define u (make-flvector nr 0.0))
  (define a (flvector (/ -8.0 3.0) 0.0 (/ 1.0 6.0) (/ 1.0 12.0)))
  (define c (case class
       [(#\A #\S #\W) 
          (flvector (/ -3.0 8.0) (/ 1.0 32.0) (/ -1.0 64.0) 0.0)]
       [else
          (flvector (/ -3.0 17.0) (/ 1.0 33.0) (/ -1.0 61.0) 0.0)]))

   (vs! nx lt1 nx_default)
   (vs! ny lt1 ny_default)
   (vs! nz lt1 nz_default)
|#

   )
  (rmpi-finish comm tc))

(define (setup lt)
  (define ng (make-fxvector (fx* 4 maxlevel)))
  (define (cidx i j)
    (fx+ (fx* i maxlevel) j))

  (for* ([j (in-range 0 3)]
         [d (in-range 0 3)])
    (fx! msg-type (fx+ (fx* d 3) j) (fx* 100 (fx+ j 3 (fx* 10 d)))))

  (fx! ng (fx+ (fx* 1 maxlevel) lt) (fxr nx lt))
  (fx! ng (fx+ (fx* 2 maxlevel) lt) (fxr ny lt))
  (fx! ng (fx+ (fx* 3 maxlevel) lt) (fxr nz lt))
  (for ([ax (in-range 1 4)])
    (fx! next ax 1)
    (for ([k 3])
      (fx! ng (fx+ (fx* ax 4) k) (fx/ (fxr ng (fx+ (fx* ax 4) (fx+1 k))) 2))))

  (for ([k (in-range lt 0 -1)])
    (fx! nx k (fxr ng (fx+ (fx* 1 maxlevel) k)))
    (fx! ny k (fxr ng (fx+ (fx* 2 maxlevel) k)))
    (fx! nz k (fxr ng (fx+ (fx* 3 maxlevel) k))))

  (define log-p (/ (log (+ nprocs 0.0001))
                   (log 2.0)))
  (define dx (quotient log_p 3))
  (define dy (quotient (- log-p dx) 2))

  (fx! pi 1 (arithmetic-shift 1 dx))
  (fx! idi 1 (modulo id (fxr pi 1)))
  (fx! pi 2 (arithmetic-shift 1 dy))
  (fx! idi 2 (modulo (fx/ id (fxr pi 1)) (fxr pi 2)))
  (fx! pi 3 (arithmetic-shift 1 dy))
  (fx! idi 3 (modulo (fx/ id (fxr pi 1)) (fxr pi 2)))

  (for ([k (in-range lt 0 -1)])
    (vector-set dead k #f)
    
    (for ([ax (in-range 1 3)])
      (v! take-ex (cidx ax k) #f)
      (v! give-ex (cidx ax k) #f)
      
      (fx! mi (cidx ax k) (fx+ 2
                               (fx- (fx/ (fx* (fx+ (fxr idi ax) 1) (fxv ng (cidx ax k))) (fxv pi ax))
                                    (fx/ (fx* (fx+ (fxr idi ax) 0) (fxv ng (cidx ax k))) (fxv pi ax))
                                    )))
      
      (fx! mi (cidx ax k) (fx+ 2
                               (fx- (fx/ (fx* (fx+ (fxv next ax) (fx+ (fxr idi ax) 1)) (fxv ng (cidx ax k))) (fxv pi ax))
                                    (fx/ (fx* (fx+ (fxv next ax) (fx+ (fxr idi ax) 0)) (fxv ng (cidx ax k))) (fxv pi ax))
                                    )))
      (when (or (= 2 (fxr mip (cidx ax k)))
                (= 2 (fxr mi (cidx ax k))))
        (fx! next ax (fx* (fvr next ax) 2)))

      (when ((fx+1 k) . < . lt)
        (when (and (= 2 (fxr mip (cidx ax k)))
                   (= 3 (fxr mi  (cidx ax k))))
          (v! give-ex (cidx ax (fx+1 k)) #t))
        (when (and (= 3 (fxr mip (cidx ax k)))
                   (= 2 (fxr mi  (cidx ax k))))
          (v! take-ex (cidx ax (fx+1 k)) #t))))

    (when (or (= (fxr mi (cidx 1 k)) 2)
              (= (fxr mi (cidx 2 k)) 2)
              (= (fxr mi (cidx 3 k)) 2))
      (v! dead k #t))

    (fx! m1 k (fxb mi (cidx 1 k)))
    (fx! m2 k (fxb mi (cidx 2 k)))
    (fx! m3 k (fxb mi (cidx 3 k)))

    (for ([ax (in-range 1 4)])
      (fx! idin (idin-idx ax 2) (modulo (fx+ (fxr idi ax) (fx+ (fxr next ax) (fxr pi ax))) (fxr pi ax)))
      (fx! idin (idin-idx ax 0) (modulo (fx- (fxr idi ax) (fx+ (fxr next ax) (fxr pi ax))) (fxr pi ax))))

    (for ([dir (list 0 2)])
      (fx! nbr (nbr-idx 1 dir k) (+ (fx* (fxr (idin-idx 1 dir)) (fxr pi 1))
                                    (fx* (fxr idi 2) (fxr pi 2))
                                    (fx* (fxr idi 3) (fxr pi 3))))

      (fx! nbr (nbr-idx 2 dir k) (+ (fx* (fxr idi 1) (fxr pi 1))
                                    (fx* (fxr (idin-idx 2 dir)) (fxr pi 2))
                                    (fx* (fxr idi 3) (fxr pi 3))))

      (fx! nbr (nbr-idx 3 dir k) (+ (fx* (fxr idi 1) (fxr pi 1))
                                    (fx* (fxr idi 2) (fxr pi 2))
                                    (fx* (fxr (idin-idx 3 dir)) (fxr pi 3))))))
  (define k lt)
  (define is1 (- (+ 2 (fxr ng (cidx 1 k)))
                 (fx/ (fx* (- (fxr pi 1) 0 (fxr idi 1)) (fxr ng (cidx 1 lt))) (fxr pi 1))))
  (define ie1 (- (+ 1 (fxr ng (cidx 1 k)))
                 (fx/ (fx* (- (fxr pi 1) 1 (fxr idi 1)) (fxr ng (cidx 1 lt))) (fxr pi 1))))
  (define n1 (+ 3 (- ie1 is1)))
     
  (define is1 (- (+ 2 (fxr ng (cidx 2 k)))
                 (fx/ (fx* (- (fxr pi 2) 0 (fxr idi 2)) (fxr ng (cidx 2 lt))) (fxr pi 2))))
  (define ie1 (- (+ 1 (fxr ng (cidx 2 k)))
                 (fx/ (fx* (- (fxr pi 2) 1 (fxr idi 2)) (fxr ng (cidx 2 lt))) (fxr pi 2))))
  (define n2 (+ 3 (- ie2 is2)))

  (define is1 (- (+ 2 (fxr ng (cidx 1 k)))
                 (fx/ (fx* (- (fxr pi 3) 0 (fxr idi 3)) (fxr ng (cidx 3 lt))) (fxr pi 3))))
  (define ie1 (- (+ 1 (fxr ng (cidx 1 k)))
                 (fx/ (fx* (- (fxr pi 3) 1 (fxr idi 3)) (fxr ng (cidx 3 lt))) (fxr pi 3))))
  (define n3 (+ 3 (- ie3 is3)))


  (fx! ir lt 1)
  (for ([j (in-range lt 0 -1)])
    (define j+1 (fx+1 j))
    (fx! ir j (fx+ (fxr ir j+1) (fx* (fxr m1 j+1) (fxr m2 j+1) (fxr m3 j+1))))))

(define (mg3P u v r a c n1 n2 n3 k)
  (for ([k (in-range lt lb -1)])
    (define j (fx-1 k))
    (rprj3 (flr r (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k)
           (flr r (fxr ir j)) (fxr m1 j) (fxr m2 j) (fxr m3 j) k))
  (let ([k lb])
    (zero3 (flr u (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k))
    (psinv (flr r (fxr ir k)) (flr u (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k) c k))
  
  (for ([k (in-range (fx+1 lb) lt)])
    (zero3 (flr u (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k))
    (interp (flr u (fxr ir j)) (fxr m1 j) (fxr m2 j) (fxr m3 j)
            (flr u (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k) k)
    (resid (flr u (fxr ir k)) (flr r (fxr ir k))
            (flr r (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k) a k)
    (psinv (flr r (fxr ir k)) (flr u (fxr ir k)) (fxr m1 k) (fxr m2 k) (fxr m3 k) c k))

  (let ([j (fx-1 lt)]
        [k lt])
    (interp (flr u (fxr ir j)) (fxr m1 j) (fxr m2 j) (fxr m3 j)
            u n1 n2 n3 k)
    (resid u v r n1 n2 n3 a k)
    (psinv r u n1 n2 n3 c k)))

(define (psinv r u n1 n2 n3 c k)
  (define (idx i j k)
    (+ (fx* i n2 n3) (fx* j n3) k))
  (define (rr v i j k)
    (flr v (idx i j k)))
  (define uidx idx)

  (with-timer T_PSINV
    (for ([i3 (in-range 2 n3)])
      (for ([i2 (in-range 2 n2)])
        (for ([i1 (in-range 1 (fx+1 n1))])
          (fl! r1 i1 (+ (rr r i1 (fx-1 i2) i3)
                        (rr r i1 (fx+1 i2) i3)
                        (rr r i1 i2 (fx-1 i3))
                        (rr r i1 i2 (fx+1 i3))))
          (fl! r2 i1 (+ (rr r i1 (fx-1 i2) (fx-1 i3))
                        (rr r i1 (fx+1 i2) (fx-1 i3))
                        (rr r i1 (fx-1 i2) (fx+1 i3))
                        (rr r i1 (fx+1 i2) (fx+1 i3)))))
        (for ([i1 (in-range 2 n1)])
          (fl! u (uidx i1 i2 i3)
               (+ (rr u i1 i2 i3)
                  (fl* (flr c 0) (rr r i1 i2 i3))
                  (fl* (flr c 1) (+ (rr (fx-1 il) i2 i3)
                                    (rr (fx+1 i1) i2 i3)
                                    (flr r1 i1)))
                  (fl* (flr c 3) (+ (flr r2 il)
                                    (flr r1 (fx-1 i1))
                                    (flr r1 (fx+1 il)))))))
                        )))
  (comm3 u n1 n2 n3 k))

(define (resid u v r n1 n2 n3 a k)
  (define (idx i j k)
    (+ (fx* i n2 n3) (fx* j n3) k))
  (define (rr v i j k)
    (flr v (idx i j k)))
  (define uidx idx)

  (with-timer T_PSINV
    (for ([i3 (in-range 2 n3)])
      (for ([i2 (in-range 2 n2)])
        (for ([i1 (in-range 1 (fx+1 n1))])
          (fl! u1 i1 (+ (rr u i1 (fx-1 i2) i3)
                        (rr u i1 (fx+1 i2) i3)
                        (rr u i1 i2 (fx-1 i3))
                        (rr u i1 i2 (fx+1 i3))))
          (fl! u2 i1 (+ (rr u i1 (fx-1 i2) (fx-1 i3))
                        (rr u i1 (fx+1 i2) (fx-1 i3))
                        (rr u i1 (fx-1 i2) (fx+1 i3))
                        (rr u i1 (fx+1 i2) (fx+1 i3)))))
        (for ([i1 (in-range 2 n1)])
          (fl! r (uidx i1 i2 i3)
               (- (rr v i1 i2 i3)
                  (fl* (flr a 0) (rr u i1 i2 i3))
                  (fl* (flr a 1) (+ (flv u2 i1)
                                    (flv u1 (fx-1 i1))
                                    (flv u1 (fx+1 i1))))

                  (fl* (flr a 3) (+ (flr u2 (fx-1 i1))
                                    (flr u2 (fx+1 il)))))))
                        )))
  (comm3 r n1 n2 n3 k))

(define (rprj3 r m1k m2k m3k s m1j m2j m3j k)
  (define (sidx i j k)
    (+ (fx* i m2j m3j) (fx* j m3j) k))
  (define (idx i j k)
    (+ (fx* i m2k m3k) (fx* j m3k) k))
  (define (rr v i j k)
    (flr v (idx i j k)))

  (with-timer T_RPRJ3
    (define d1 (if (= m1k 3) 2 1))
    (define d1 (if (= m2k 3) 2 1))
    (define d1 (if (= m3k 3) 2 1))

    (for ([j3 (in-range 2 m3j)])
      (define i3 (fx- (fx* 2 j2) d2))
      (for ([j2 (in-range 2 m2j)])
        (define i2 (f- (fx* 2 j2) d2))
        (for ([j1 (in-range 2 (fx+1 m1j))])
          (define i1 (fx- (fx* 2 j1) d1))
          (fl! x1 (fx-1 i1) (+ (rr r (fx-1 i1) (fx-1 i2) i3)
                               (rr r (fx-1 i1) (fx+1 i2) i3)
                               (rr r (fx-1 i1) i2 (fx-1 i3))
                               (rr r (fx-1 i1) i2 (fx+1 i3))))
          (fl! y1 (fx-1 i1) (+ (rr r (fx-1 i1) (fx-1 i2) (fx-1 i3))
                               (rr r (fx-1 i1) (fx-1 i2) (fx+1 i3))
                               (rr r (fx-1 i1) (fx+1 i2) (fx-1 i3))
                               (rr r (fx-1 i1) (fx+1 i2) (fx+1 i3)))))
        (for ([j1 (in-range 2 m1j)])
          (define i1 (rx- (fx* 2 j1) d1))
          (fl! y1 (fx-1 i1) (+ (rr r i1 (fx-1 i2) (fx-1 i3))
                               (rr r i1 (fx-1 i2) (fx+1 i3))
                               (rr r i1 (fx+1 i2) (fx-1 i3))
                               (rr r i1 (fx+1 i2) (fx+1 i3))))
          (fl! x1 (fx-1 i1) (+ (rr r i1 (fx-1 i2) i3)
                               (rr r i1 (fx+1 i2) i3)
                               (rr r i1 i2 (fx-1 i3))
                               (rr r i1 i2 (fx+1 i3))))
          (fl! s (sidx j1 j2 j3)
               (+ (fl* 0.5 (rr r i1 i2 i3))
                  (fl* 0.25 (fl+ (rr r (fx-1 i1) i2 i3) (rr r (fx+1 i1) i2 i3) x2))
                  (fl* 0.125 (fl+ (flr x1 (fx-1 i1)) (flr x1 (fx+1 i1)) y2))
                  (fl* 0.0625 (fl+ (flr y1 (fx-1 i1)) (flr y1 (fx+1 i1))))))))))
  (define j (fx-1 k))
  (comm3 s m1j m2j m3j j))

(define (interp z mm1 mm2 mm3 u n1 n2 n3 k)
  (with-timer T_INTERP
    (cond
      [(and (!= n1 3) (!= n2 3) (!= n3 3))
       (for* ([i3 (in-range 1 mm3)]
              [i2 (in-range 1 mm2)])
         (for ([i1 (in-range 1 (fx+1 mm1))])
           (fl! z1 i1 (fl+ (zr i1 (fx+1 i2) i3) (zr i1 i2 i3)))
           (fl! z2 i1 (fl+ (zr i1 i2 (fx+1 i3)) (zr i1 i2 i3)))
           (fl! z3 i1 (fl+ (zr i1 (fx+1 i2) (fx+1 i3)) (zr i1 i2 (fx+1 i3)) (flr z1 i1))))

         (for ([i1 (in-range 1 mm1)])
           (define ui1 (uidx (fx-1 (fx* 2 i1)) (fx-1 (fx* 2 i2)) (fx-1 (fx* 2 i3))))
           (define ui2 (uidx (fx* 2 i1) (fx-1 (fx* 2 i2)) (fx-1 (fx* 2 i3))))
           (define zv (zr i1 i2 i3))
           (fl! u ui1 (fl+ (flr u ui1) zr))
           (fl! u ui2 (fl+ (flr u ui2) (fl* 0.5 (fl* (zr (fx+1 i1) i2 i3) zr)))))

         (for ([i1 (in-range 1 mm1)])
           (define ui1 (uidx (fx-1 (fx* 2 i1)) (fx* 2 i2) (fx-1 (fx* 2 i3))))
           (define ui2 (uidx (fx* 2 i1) (fx* 2 i2) (fx-1 (fx* 2 i3))))
           (fl! u ui1 (fl+ (flr u ui1) (fl* 0.5 (flr z1 i1))))
           (fl! u ui2 (fl+ (flr u ui2) (fl* 0.25 (fl+ (flr z1 i1) (flr z1 (fx+1 i1)))))))

         (for ([i1 (in-range 1 mm1)])
           (define ui1 (uidx (fx-1 (fx* 2 i1)) (fx-1 (fx* 2 i2)) (fx* 2 i3)))
           (define ui2 (uidx (fx* 2 i1) (fx-1 (fx* 2 i2)) (fx* 2 i3)))
           (fl! u ui1 (fl+ (flr u ui1) (fl* 0.5 (flr z2 i1))))
           (fl! u ui2 (fl+ (flr u ui2) (fl* 0.25 (fl+ (flr z2 i1) (flr z2 (fx+1 i1)))))))

         (for ([i1 (in-range 1 mm1)])
           (define ui1 (uidx (fx-1 (fx* 2 i1)) (fx* 2 i2) (fx-1 (fx* 2 i3))))
           (define ui2 (uidx (fx* 2 i1) (fx* 2 i2) (fx* 2 i3)))
           (fl! u ui1 (fl+ (flr u ui1) (fl* 0.25 (flr z3 i1))))
           (fl! u ui2 (fl+ (flr u ui2) (fl* 0.125 (fl+ (flr z3 i1) (flr z3 (fx+1 i1))))))))]
      [else
       (define-values (d1 t1) (if (= n1 3) (values 2 1) (values 1 0)))
       (define-values (d2 t2) (if (= n2 3) (values 2 1) (values 1 0)))
       (define-values (d3 t3) (if (= n3 3) (values 2 1) (values 1 0)))
       (for ([i3 (in-range 1 mm3)])
         (for ([i2 (in-range 1 mm2)])
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) d1) (fx- (fx* 2 i2) d2) (fx- (fx* 2 i3) d3)))
             (fl! u ui1 (fl+ (flr ui1) (zr i1 i2 i3))))
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) t1) (fx- (fx* 2 i2) d2) (fx- (fx* 2 i3) d3))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.5 (fl+ (zr (fx+1 i1) i2 i3) (zr i1 i2 i3))))))))
         (for ([i2 (in-range 1 mm2)])
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) d1) (fx- )(fx* 2 i2) t2) (fx- (fx* 2 i3) d3))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.5 (fl+ (zr i1 (fx+1 i2) i3) (zr i1 i2 i3))))))
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) t1) (fx- (fx* 2 i2) t2) (fx- (fx* 2 i3) d3)))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.25 (fl+ (zr (fx+1 i1) (fx+1 i2) i3) 
                                                     (zr (fx+1 i1) i2 i3) 
                                                     (zr i1 (fx+1 i2) i3)
                                                     (zr i1 i2 i3))))))))

       (for ([i3 (in-range 1 mm3)])
         (for ([i2 (in-range 1 mm2)])
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) d1) (fx- (fx* 2 i2) d2) (fx- (fx* 2 i3) t3)))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.5 (fl+ (zr i1 i2 (fx+1 i3)) (zr i1 i2 i3))))))
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) t1) (fx- (fx* 2 i2) d2) (fx- (fx* 2 i3) t3)))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.25 (fl+ (zr (fx+1 i1) i2 (fx+1 i3)) 
                                                      (zr i1 i2 (fx+1 i3)) 
                                                      (zr (fx+1 i1) i2 i3) 
                                                      (zr i1 i2 i3)))))))
         (for ([i2 (in-range 1 mm2)])
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) d1) (fx- )(fx* 2 i2) t2) (fx- (fx* 2 i3) t3))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.25 (fl+ (zr i1 (fx+1 i2) (fx+1 i3))
                                                      (zr i1 i2 (fx+1 i3))
                                                      (zr i1 (fx+1 i2) i3)
                                                      (zr i1 i2 i3)
                                                      )))))
           (for ([i1 (in-range 1 mm1)])
             (define ui1 (uidx (fx- (fx* 2 i1) t1) (fx- (fx* 2 i2) t2) (fx- (fx* 2 i3) t3)))
             (fl! u ui1 (fl+ (flr ui1) (fl* 0.125 (fl+ (zr (fx+1 i1) (fx+1 i2) (fx+1 i3)) 
                                                       (zr (fx+1 i1) i2 (fx+1 i3)) 
                                                       (zr i1 (fx+1 i2) (fx+1 i3)) 
                                                       (zr i1 i2 (fx+1 i3)) 
                                                       (zr (fx+1 i1) (fx+1 i2) i3) 
                                                       (zr (fx+1 i1) i2 i3) 
                                                       (zr i1 (fx+1 i2) i3) 
                                                       (zr i1 i2 i3))))))))]))
  (comm3_ex u n1 n2 n3 k))

(define (norm2u3 r n1 n2 n3 nx0 ny0 nz0)
  (define-values (s rnmu)
    (with-timer T_NROM2U3
      (define dn (* 1.0 nx0 ny0 nz0))
      (for*/fold ([s 0.0]
                  [rnmu 0.0])
                 ([i3 (in-range 2 n3)]
                  [i2 (in-range 2 n2)]
                  [i1 (in-range 2 n1)])
        (define rv (flr r (+ (* i1 n2 n3) (fx* i2 n3) i3)))
        (values
          (fl+ s (expt rv 2))
          (max rnmu (flabs rv))))))
  (with-timer T_RCOMM
    (values
      (rmpi-all-reduce comm rnmu +)
      (fl-sqrt (fl/ (rmpi-allreduce comm s +)
                    dn)))))

(define (comm3 u n1 n2 n3 kk)
  (cond
    [(not (vector-ref dead kk))
     (for ([axis (in-range 1 4)])
       (cond 
         [(!= nprocs 1)
           (ready axis 0 kk)
           (ready axis 2 kk)
           
           (give3 axis 2 u n1 n2 n3 kk)
           (give3 axis 0 u n1 n2 n3 kk)
           
           (take3 axis 0 u n1 n2 n3)
           (take3 axis 2 u n1 n2 n3)]
         [else
           (comm1p axis u n1 n2 n3 kk)]))]
    [else
      (zero u n1 n2 n3)]))

(define (comm3_ex u n1 n2 n3 kk)
  (for ([axis (in-range 1 4)])
    (cond
      [(!= nprocs 1)
       (when (take-ex axis kk)
         (ready axis 0 kk)
         (ready axis 2 kk)
         (take3-ex axis 0 u n1 n2 n3)
         (take3-ex axis 2 u n1 n2 n3))

       (when (give-ex axis kk)
         (give-ex axis 2 u n1 n2 n3 kk)
         (give-ex axis 0 u n1 n2 n3 kk))]
      [else
        (comm1p-ex axis u n1 n2 n3 kk)])))

(define (ready axis dir k)
  (define buff-id (+ dir 3))
  (define buff-len nm2)
  
  (for ([i (in-range 1 (fx+1 nm2))])
    (fl! (buff (fx+ (fx* i XX) buff-id)) 0.0))
  
  (with-timer
    (fx! ()))
  )

(define (give3 axis dir u n1 n2 n3 k)
  (define (do-send IO II uu XX YY ZZ)
    (define UU (if (= dir -1) 2 (fx-1 uu)))
    (for ([io (in-range 2 IO)])
      (for ([ii (in-range 2 II)]
            [buff-len (in-naturals 1)])
        (define idx (+ (* XX n2 n3) (* YY n3) ZZ))
        (fl! buff (fx+ (fx* buf-len 55) buf-id) (flr u idx))))
    (rmpi-send comm (fxr nbr (+ (* axis 3 K) (fx* dir K) k)) buff))


  (case axis
    [(1) (do-send n3 n2 n1 UU ii io)]
    [(2) (do-send n3 n1 n2 ii UU io)]
    [(3) (do-send n2 n1 n3 ii io UU)]
    ))

(define (take3 axis dir u n1 n2 n3)
  (define buff (with-timer T_COMM3 (rmpi-recv comm RECV)))

  (define-syntax-rule (do-recv IO II uu XX YY ZZ)
    (define UU (if (= dir -1) uu 1))
    (for ([io (in-range 2 IO)])
      (for ([ii (in-range 2 II)]
            [idx (in-naturals 1)])
        (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
        (fl! u IDX (flr buff (fx+ (fx* idx 55) buf-id))))))

  (case axis
    [(1) (do-recv n3 n2 n1 UU ii io)]
    [(2) (do-recv n3 n1 n2 ii UU io)]
    [(3) (do-recv n2 n1 n3 ii io UU)]
    ))

(define (give3-ex axis dir u n1 n2 n3 k)
  (define buff-id (+ dir 2))
  
  (define-syntax-rule (do-send IO II RO RM RI III XX YY ZZ)
    (rmpi-send comm (fxr nbr (+ (* axis 3 K) (fx* dir K) k))
      (cond
        [(= dir -1)
           (define UU 2)
           (for*/flvector ([io IO]
                           [ii II])
             (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
             (flr u IDX))]
        [else
           (define UU III)
           (for*/flvector ([io RO]
                           [im RM]
                           [ii RI])
             (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
             (flr u IDX))])))

  (case axis
    [(1) (do-send n3 n2 (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n2)) (in-range (fx-1 n1) (fx+1 n1)) i1 UU i2 i3)]
    [(2) (do-send n3 n1 (in-range 1 (fx+1 n3)) (in-range (fx-1 n2) (fx+1 n2)) (in-range 1 (fx+1 n1)) i2 i1 UU i3)]
    [(3) (do-send n2 n1 (in-range (fx-1 n3) (fx+1 n3)) (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) i3 i1 i2 UU)]))

(define (take3-ex axis dir u n1 n2 n3)
  (define-syntax-rule (do-recv IO II RO RM RI AA BB XX YY ZZ)
    (define buff (mpi-recv comm PEER))

    (cond
      [(= dir -1)
       (define UU AA)
       (for*/fold ([i 0]) ([io IO]
                           [ii II])
         (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
         (fl! u IDX (flr buff i))
         (fx+1 i))]
      [else
       (define UU BB)
       (for*/fold ([i 0]) ([io RO]
                           [im RM]
                           [ii RM])
         (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
         (fl! u IDX (flr buff i))
         (fx+1 i))]))

  (case axis
    [(1) (do-recv n3 n2 (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n2)) (in-range 1 3) n1 i1 UU i2 i3)]
    [(2) (do-recv n3 n1 (in-range 1 (fx+1 n3)) (in-range 1 3) (in-range 1 (fx+1 n1)) n2 i2 i1 UU i3)]
    [(3) (do-recv n2 n1 (in-range 1 3) (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) n3 i3 i1 i2 UU)]))

(define (comm1p axis u n1 n2 n3 kk)
  (define-syntax-rule (copy-out IO II XX YY ZZ)
    (for*/flvector ([io IO]
                    [ii II])
      (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
      (flv u IDX))

  (define-syntax-rule (copy-in v IO II XX YY ZZ)
    (for*/fold ([i 0]) ([io IO]
                        [ii II])
      (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
      (fl! u IDX (flr v i)))

  (define dir1buff
    (case axis
      [(1) (copy-out n3 n2 (in-range 2 n3) (in-range 2 n2) (fx-1 n1) ii io)]
      [(2) (copy-out n3 n1 (in-range 2 n3) (in-range 1 (fx+1 n1)) ii (fx-1 n2) io)]
      [(3) (copy-out n2 n1 (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii io (fx-1 n3))]))

  (define dir2buff
    (case axis
      [(1) (copy-out n3 n2 (in-range 2 n3) (in-range 2 n2) 2 ii io)]
      [(2) (copy-out n3 n1 (in-range 2 n3) (in-range 1 (fx+1 n1)) ii 2 io)]
      [(3) (copy-out n2 n1 (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii io 2)]))

  (case axis
    [(1) (copy-in dir1buff n3 n2 (in-range 2 n3) (in-range 2 n2) n1 ii io)]
    [(2) (copy-in dir1buff n3 n1 (in-range 2 n3) (in-range 1 (fx+1 n1)) ii n2 io)]
    [(3) (copy-in dir1buff n2 n1 (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii io n3)]))

  (case axis
    [(1) (copy-in dir1buff n3 n2 (in-range 2 n3) (in-range 2 n2) 1 ii io)]
    [(2) (copy-in dir1buff n3 n1 (in-range 2 n3) (in-range 1 (fx+1 n1)) ii 1 io)]
    [(3) (copy-in dir1buff n2 n1 (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii io 1)])))

(define (comm1p-ex axis u n1 n2 n3 kk)
  (define-syntax-rule (copy-zero IO II XX YY ZZ)
    (for*/fold ([i 0]) ([io IO]
                        [ii II])
      (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
      (fl! u IDX 0.0)))

  (define-syntax-rule (copy-zero2 IO IM II XX YY ZZ)
    (for*/fold ([i 0]) ([io IO]
                        [im IM]
                        [ii II])
      (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
      (fl! u IDX 0.0)))

  (define-syntax-rule (copy-out IO IM II XX YY ZZ)
    (for*/flvector ([io IO]
                    [im IM]
                    [ii II])
      (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
      (flv u IDX)))

  (define-syntax-rule (copy-out IO II XX YY ZZ)
    (for*/flvector ([io IO]
                    [ii II])
      (define IDX (+ (* XX n2 n3) (* YY n3) ZZ))
      (flv u IDX))

   (take-ex axis kk)
    (case axis
      [(1) (copy-zero (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n2)) n1 ii io)]
      [(2) (copy-zero (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n1)) ii n2 io)]
      [(3) (copy-zero (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii io n3)])

    (case axis
      [(1) (copy-zero2 (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n2)) (in-range 1 3) ii im io)]
      [(2) (copy-zero2 (in-range 1 (fx+1 n3)) (in-range 1 3) (in-range 1 (fx+1 n1)) ii im io)]
      [(3) (copy-zero2 (in-range 1 3) (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii im io)]))

  (when (give-ex axis kk)
    (vector-set! buffs 3
      (case axis
        [(1) (copy-out (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n2)) (in-range (fx-1 n1) (fx+1 n1)) ii im io)]
        [(2) (copy-out (in-range 1 (fx+1 n3)) (in-range (fx-1 n2) (fx+1 n2)) (in-range 1 (fx+1 n1)) ii im io)]
        [(3) (copy-out (in-range (fx-1 n3) (fx+1 n3)) (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii im io)]))

    (vector-set! buffs 1
      (case axis
        [(1) (copy-out1 (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n2)) 2 ii io)]
        [(2) (copy-out1 (in-range 1 (fx+1 n3)) (in-range 1 (fx+1 n1)) ii 2 io)]
        [(3) (copy-out1 (in-range 1 (fx+1 n2)) (in-range 1 (fx+1 n1)) ii io 2)])))

  (match-define (vector v1 v2 v3 v4) buffs)
  (for ([i (in-range 1 (fx+1 nm2))])
    (fl! v4 i (flr v3 i))
    (fl! v2 i (flr v1 i))))


(define (zran3 z n1 n2 n3 nx ny k)
  (define a1 (power a nx 1 0))
  (define a2 (power a nx ny 0))
  
  (zero3 z n1 n2 n3)
  
  (define ai (power a nx (fx+ (fx- is2 2) (fx* ny (fx- is3 2))) (fx- is1 2)))
  (define d1 (fx+ (fx- ie1 is1) 1))
  (define e1 (fx+ (fx- ie1 is1) 2))
  (define e2 (fx+ (fx- ie2 is2) 2))
  (define e3 (fx+ (fx- ie3 is3) 2))
  (define x0 (randlc/a 314159265.0 ai))

  (for/fold ([x0 x0]) ([i3 (in-range 2 (fx+1 e3))])
    (define x1 x0)
    (for/fold ([x1 x0]) ([i2  (in-range 2 (fx+1 e2))])
      (vrandlc d1 xx a (fxl z GGGG))
      (randlc/a x1 a2))
    (randlc/a x0 a2))

  (for ([i (in-range 1 mm)])
    (define idx (fx+1 (fx* i 2)))
    (define idx2 (fx* i 2))
    (fl! ten idx 0.0)
    (fx! j1 idx 0)
    (fx! j2 idx 0)
    (fx! j3 idx 0)
    (fl! ten idx2 0.0)
    (fx! j1 idx2 0)
    (fx! j2 idx2 0)
    (fx! j3 idx2 0))

  (for ([i3 (in-range 2 n3)]
        [i2 (in-range 2 n2)]
        [i1 (in-range 2 n1)])
    (define IDX (+ (* i1 n2 n3) (* i2 n3) i3))
    (define zv (flr z IDX))
    (when (> zv (flr ten 3))
      (fl! ten 3 zr)
      (fx! j1 3 i1)
      (fx! j2 3 i2)
      (fx! j3 3 i3)
      (bubble ten j1 j2 j3 mm 1))
    (when (< zv (flr ten 2))
      (fl! ten 2 zr)
      (fx! j1 2 i1)
      (fx! j2 2 i2)
      (fx! j3 2 i3)
      (bubble ten j1 j2 j3 mm 0)))

  (rmpi-barrier comm)

  (define i1 mm)
  (define i0 mm)

  (define (jxr v i j)
    (fxr v (fx+ (fx* 2 i) j)))

  (define (jg! i j k v)
    (fx! jg (+ (* i (fx+1 mm) 2) (fx* j 2) k) v))

  (define (jgr i j k)
    (fxr jg (+ (* i (fx+1 mm) 2) (fx* j 2) k)))

  (define (zidx i1 i2 i3) (+ (* i1 n2 n3) (fx* i2 n3) i3))

  (for ([i (in-range mm 0 -1)])
    (define zv1 (flr z (zidx (jxr j1 i1 1) (jxr j2 i1 1) (jxr j3 i1 1))))
    (define best (rmpi-allreduce comm max zv1))
    (cond 
     [(= best zv1)
       (jg! 0 i 1 id)
       (jg! 1 i 1 (fx+ (fx- is1 2) (jxr j1 i1 1)))
       (jg! 2 i 1 (fx+ (fx- is2 2) (jxr j2 i1 1)))
       (jg! 3 i 1 (fx+ (fx- is3 2) (jxr j3 i1 1)))
       (set! i1 (fx-1 i1))]
     [else
       (jg! 0 i 1 0)
       (jg! 1 i 1 0)
       (jg! 2 i 1 0)
       (jg! 3 i 1 0)])
    (fl! ten (fx+ (fx* i 2) 1) best)
    (for ([x (in-fxvector (rmpi-allreduce comm max (for/fxvector ([j 4]) (jgr j i 1))))]
          [j (in-naturals)])
      (jg! j i 1 x))

    (define zv2 (flr z (zidx (jxr j1 i1 0) (jxr j2 i1 0) (jxr j3 i1 0))))
    (define best2 (rmpi-allreduce comm max zv2))
    (cond 
     [(= best2 zv2)
       (jg! 0 i 0 id)
       (jg! 1 i 0 (fx+ (fx- is1 2) (jxr j1 i1 0)))
       (jg! 2 i 0 (fx+ (fx- is2 2) (jxr j2 i1 0)))
       (jg! 3 i 0 (fx+ (fx- is3 2) (jxr j3 i1 0)))
       (set! i0 (fx-1 i0))]
     [else
       (jg! 0 i 0 0)
       (jg! 1 i 0 0)
       (jg! 2 i 0 0)
       (jg! 3 i 0 0)])
    (fl! ten (fx+ (fx* i 2) 0) best2)
    (for ([x (in-fxvector (rmpi-allreduce comm max (for/fxvector ([j 4]) (jgr j i 0))))]
          [j (in-naturals)])
      (jg! j i 0 x)))

    (define m1 (fx+1 i1))
    (define m0 (fx+1 i0))

    (rmpi-barrier comm)

    (for* ([i3 (in-range 1 (fx+1 n3))]
           [i2 (in-range 1 (fx+1 n2))]
           [i1 (in-range 1 (fx+1 n1))])
      (fl! z (zidx i1 i2 i3) 0.0))

    (for ([i (in-range mm (fx-1 m0) -1)])
      (fl! z (zidx (jxr j1 i 0) (jxr j2 i 0) (jxr j3 i 0)) -1.0))

    (for ([i (in-range mm (fx-1 m1) -1)])
      (fl! z (zidx (jxr j1 i 1) (jxr j2 i 1) (jxr j3 i 1)) 1.0))

    (comm z n1 n2 n3 k))


(define (power a n1 n2 n3)(zr (fx+1 i1) (fx+1 i2) i3) 

  (define n1j n1)
  (let loop ([aj a]
             [nj n3]
             [n2j n2]
             [power 1.0])
    (cond
      [(= nj 0) power]
      [else
        (let* ([t (> n2j 0)]
               [nj (if t (+ nj n1j) nj)]
               [n2j (if t (quotient n2j 2) n2j)])
          (loop
            (randlc/a aj aj)
            (quotient nj 2)
            nj2
            (if (= (modulo nj 2) 1) (randlc/a power aj) power))
          )]
      ))
   )

(define (bubble ten j1 j2 j3 m ind)
  (define op (if (= ind 1) > <))
  (define (ten! i j v) (fl! ten (fx+ (fx* i 2) j) v))
  (define (tenr i j) (flv ten (fx+ (fx* i 2) j)))
  (define (swap v i1 i2)
    (define i1x (fx+ (fx* i1 2) ind))
    (define i2x (fx+ (fx* i2 2) ind))
    (define t (flv v i1x))
    (fl! v i1x (flr v i2x))
    (fl! v i2x t))
  
  (let/ec k
    (for ([i (in-range 1 m)])
      (cond 
        [(op (tenv ten i ind) (tenv (fx+1 i) ind))
         (swap ten (fx+1 i) i)
         (swap j1 (fx+1 i) i)
         (swap j2 (fx+1 i) i)
         (swap j3 (fx+1 i) i)]
        [else (k)]))))

(define (zero3 z n1 n2 n3)
  (for ([i3 (in-range 1 (fx+1 n3))]
        [i2 (in-range 1 (fx+1 n2))]
        [i1 (in-range 1 (fx+1 n1))])
       (fl! v (fx+ (fx+ (fx* (fx* n2 n1) i3)
                        (fx* n1 i2))
                   i1) 0.0)))

(define (main . argv) 
  (let* ([args (parse-cmd-line-args argv "Multi-Grid")]
         [class (BMArgs-class args)]
         [serial (BMArgs-serial args)]
         [num-threads (BMArgs-num-threads args)]
         [bmname "CG"])

    (print-banner "Multi-Grid" args) 

    (define place-args (list class bmname num-threads))

    (rmpi-launch
      (rmpi-build-default-config
        #:mpi-module (quote-module-path)
        #:mpi-func 'mg-place
        #:mpi-args place-args)
      (rmpi-make-localhost-config num-threads 6341 'mg))))

(define-syntax-rule (v++! v idx)
  (vs! v idx (fx+ (vr v idx) 1)))

(define-syntax-rule (fx++! v idx)
  (fx! v idx (fx+ (fxr v idx) 1)))

(define-syntax-rule (fx--! v idx)
  (fx! v idx (fx- (fxr v idx) 1)))
