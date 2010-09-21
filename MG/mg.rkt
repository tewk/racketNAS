#lang racket/base

(provide main)
  
(require "../bm-args.ss") 
(require "../bm-results.ss") 
(require "../rand-generator.ss")
(require "../timer.ss")
(require "../parallel-utils.ss")
(require racket/match)
(require racket/math)
(require racket/place)
(require racket/place-utils)
(require (for-syntax scheme/base))

(require (only-in scheme/flonum make-flvector make-shared-flvector shared-flvector flvector-length))

#|
(require scheme/fixnum scheme/flonum)

(require (only-in scheme/flonum make-flvector make-shared-flvector)
         scheme/require (for-syntax scheme/base)
   (filtered-in
    (lambda (name) (regexp-replace #rx"unsafe-" name ""))
    scheme/unsafe/ops))
(define vr vector-ref)
(define vs! vector-set!)
(define flvs! flvector-set!)
(define flvr flvector-ref)
(require (rename-in scheme/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-flvector-ref flvr] 
                    [unsafe-flvector-set! flvs!]))
|#

#|
|#
(require (rename-in scheme/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-flvector-ref flvr] 
                    [unsafe-flvector-set! flvs!]
                    [unsafe-fl+ fl+]
                    [unsafe-fl- fl-]
                    [unsafe-fl* fl*]
                    [unsafe-fl/ fl/]
))
 
(define (comgrp-wait grp)
  (match grp
    [(list 0 pls np) (for ([ch pls]) (place-channel-send ch 0))]
    [(list n ch np)  (place-channel-recv ch)]))
(define (comgrp-tell grp)
  (match grp
    [(list 0 pls np) (for ([ch pls]) (place-channel-recv ch))]
    [(list n ch np)  (place-channel-send ch 1)]))

(define-syntax-rule (when0 np body ...)
   (unless (and (pair? np) (not (= (car np) 0)))
    body ...))

(define (barrier2 grp)
  (when (pair? grp)
  (comgrp-tell grp)
  (comgrp-wait grp)))

(define DOPLACES #t)

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
    [(#\S) (values  32  32  32  4 5 5 5 5 5)] 
    [(#\W) (values  64  64  64 40 6 6 6 6 6)]
    [(#\A) (values 256 256 256  4 8 8 8 8 8)]
    [(#\B) (values 256 256 256 20 8 8 8 8 8)] 
    [(#\C) (values 512 512 512 20 9 9 9 9 9)]
    [else (error "Unknown class")]))

(define (main . argv) 
  (let ([args (parse-cmd-line-args argv "Conjugate Gradient")]) 
    (run-benchmark args)))

(define cgitmax 25)
(define make-fxvector make-vector)

(define (run-benchmark args) 
  (define maxlevel 11)
  (let ([bmname "CG"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(nx_default ny_default nz_default nit_default lm lt_default ndim1 ndim2 ndim3) (get-class-size CLASS)])
    (let* (
          [nx (make-shared-flvector maxlevel 0.0)]
          [ny (make-shared-flvector maxlevel 0.0)]
          [nz (make-shared-flvector maxlevel 0.0)]
          [ir (make-shared-flvector maxlevel 0.0)]
          [m1 (make-shared-flvector maxlevel 0.0)]
          [m2 (make-shared-flvector maxlevel 0.0)]
          [m3 (make-shared-flvector maxlevel 0.0)]
          [lt lt_default]
          [nit nit_default]
          [nm (+ (arithmetic-shift 1 lm) 2)]
          [nv (* (+ 2 (arithmetic-shift 1 ndim1)) 
                 (+ 2 (arithmetic-shift 1 ndim2)) 
                 (+ 2 (arithmetic-shift 1 ndim3)))]
          [nr (/ (* 8 (+ nv (* nm nm) (* 5 nm) (* 7 lm))) 7)]
          [nm2 (* 2 nm nm)]
          [r (make-shared-flvector nr 0.0)]
          [v (make-shared-flvector nv 0.0)]
          [u (make-shared-flvector nr 0.0)]
          [a (shared-flvector maxlevel (/ -8.0 3.0) 0.0 (/ 1.0 6.0) (/ 1.0 12.0))]
          [c (case CLASS 
               [(#\A #\S #\W) 
                  (shared-flvector maxlevel (/ -3.0 8.0) (/ 1.0 32.0) (/ -1.0 64.0) 0.0)]
               [else
                  (shared-flvector maxlevel (/ -3.0 17.0) (/ 1.0 33.0) (/ -1.0 61.0) 0.0)])])
        (let-values ([(n1 n2 n3 is1 is2 is3 ie1 ie2 ie3) (setup)])
      


      (zero3 u 0 n1 n2 n3)
      (zran3 v n1 n2 n3 (vr nx (sub1 lt)) (vr ny (sub1 lt)))

      (resid uv r 0 n1 n2 n3)

      (timer-stop 1)
      (let ([verified (verify zeta zeta-verify-value)])
        (print-banner "Conjugate Gradient" args) 
        (printf "Size = ~a niter = ~a~n" na niter) 
        (if serial 
            (printf "SERIAL~n")
            (printf "PARALLEL~n"))
        (if verified 
            (printf "Verification succeeded~n") 
            (printf "Verification failed~n"))
        (let* ([time (/ (read-timer 1) 1000)]
               [results (new-BMResults bmname CLASS na 0 0  niter time 
                                       (get-mflops na nonzer time niter)
                                       "floating point" 
                                       (if verified 1 0)
                                       serial 
                                       num-threads 
                                       -1)]) 
            (print-results results) 
            (when #f (print-timers))))))))

(define-syntax-rule (vr3 v i1 i2 i3 n1 n2) (vr v (+ i1 (* n1 (+ i2 (* n2 i3))))))

(define (verify class rnm2)
  (define verify-value 
    (case class
      [(#\S) 0.530770700573E-4]
      [(#\W) 0.250391406439E-17]
      [(#\A) 0.2433365309E-5]
      [(#\B) 0.180056440132E-5]
      [(#\C) 0.570674826298E-6]))

  (printf " L2 Norm is ~a\n" rnm2)

  (let ([deviation (abs (- verify_value rnm2))])
    (if (deviation . < . 1.0E-8)
      (begin
        (printf " Deviation is   ~a\n" deviation)
        #t)
      (begin
        (printf " The correct L2 Norm is ~n\n" verify-value)
        #f))))

(define (get-mflops tm niter n1 n2 n3)
  (if (not (= total-time 0.0))
    (* (* 58.0 n1 n2 n3)
       (/ niter (* total-time 1000000.0)))
    0.0))

(define (get-input-pars)
  (define fn "mg.input")
  (if (file esists fn)
    (match (call-with-input-file fn read)
      [(list lt lnx lny lnz nit)
        (when (lt . > . maxlevel)
          (printf "lt=~a Maximum allowable=~a\n" lt maxleevel)
          (exit 0))
        (values nit lt lx lz)]
      [else 
        (printf "Error reading from file mg.input\n")
        (exit 0)])
    (printf "No input file mg.input, Using compiled defaults\n")))

(define (setup ir nx ny nz m1 m2 m3 )
  (define size1 3)
  (define size2 10)
  (define mi (make-fxvector (* size1 size2) 0.0))
  (define ng (make-fxvector (* size1 size2) 0.0))
  (define lt1 (sub1 lt))

  (for* ([ax (in-range size1)]
         [k (in-range (- lt 2) -1 -1)])
    (let ([k1 (add1 k)])
        (vs! nx k (/ (vr nx k1) 2))
        (vs! ny k (/ (vr ny k1) 2))
        (vs! nz k (/ (vr nz k1) 2))))
    
    (vs! ng (+ ax (* k size1)) (/ (vr ng (+ ax (* (add1 k) size1))) 2))

  (for ([k (in-range (- lt 2) -1 -1)])
      (let ([ksize1 (* k size1)])
        (vs! nx k (vr ng ksize1))
        (vs! ny k (vr ng (+ ksize1 1)))
        (vs! nz k (vr ng (+ ksize1 2)))))

  (for ([k (in-range (- lt -1) -1 -1)])
    (vs! m1 k (* 2 (vr ng (* k size1))))
    (vs! m2 k (* 2 (vr ng (+ (* k size1) 1))))
    (vs! m3 k (* 2 (vr ng (+ (* k size1) 2)))))

  (define k (sub1 lt))
  (define-syntax-rule (_x off)
    (vs! nsizies off (- (+ 3 (+ 1 (vr ng (+ off (* k size1)))) (- (+ 2 (vr nv (+ off (* k size1)))) (vr ng (+ off (* k size1))))))))

  (_x 0)
  (_x 1)
  (_x 2)

  (for ([j (in-range (- lt 2) -1 -1)])
    (let ([j1 (add1 j)])
    (vs! ir j (+ (vr ir j1) (* (vr m1 j1) (vr m2 j1) (vr m3 j1))))))))

(define (zero3 z off n1 n2 n3)
  (for* ([i3 (in-range n3)]
         [i2 (in-range n2)]
         [i1 (in-range n1)])
    (vs! z (+ off i1 (* n1 (+ i2 (* n2 i3)))) 0.0)))

(define (zran3 z n1 n2 n3 nx ny)
  (define j1 (make-fxvector (* 2 mm) 0))
  (define j2 (make-fxvector (* 2 mm) 0))
  (define j3 (make-fxvector (* 2 mm) 0))
  (define jg (make-fxvector (* 4 2 mm) 0))
  (define-syntax-rule (zidx i1 i2 i3) (+ i1 (* n1 (+ i2 (* n2 i3)))))
  (zero3 z 0 n1 n2 n3)

  (let* ([i (+ (- is1 2) (* nx (+ (- is2 2) (* ny (- is3 2)))))]
        [d1 (add1 (- ie1 is1))]
        [e1 (+ 2 (- ie1 is1))]
        [e2 (+ 2 (- ie2 is2))]
        [e3 (+ 2 (- ie3 is3))]
        [a (expt 5.0 13)]
        [rng (new-rand)]
        [a1 (power/r rng a nx)]
        [a2 (power/r rng a (* nx ny))]
        [ai (power/r rng a i)])
    (for ([x0 (randlc/r rng 3.14159265.0 ai)])  ([i3 (in-range 2 (add1 e3))])
      (for/fold ([x1 x0]) ([i2 (in-range 2 (add1 e2))])
        (let ([xx x1])
          (vranlc/r rng d1 xx a z (+ 1 (* n1 (+ (- i2 1) (* n2 (- i3 1))))))
          (ranlc/r rng x1 a1)))
      (randlc/r rng x0 a2)))

  (for ([i (in-range mm)])
    (let ([imm (+ i mm)])
      (vs! ten imm 0.0)
      (vs! j1 imm 0.0)
      (vs! j2 imm 0.0)
      (vs! j3 imm 0.0)
      (vs! ten i 0.0)
      (vs! j1 i 0.0)
      (vs! j2 i 0.0)
      (vs! j3 i 0.0)))
  
  (for* ([i3 (in-range 1 (sub1 n3))]
         [i2 (in-range 1 (sub1 n2))]
         [i1 (in-range 1 (sub1 n1))])
    (let ([zv (vr z (+ i1 (* n1 (+ i2 (* n2 i3)))))])
      (if (zv . > . (vr ten mm))
        (begin
        (vs! ten imm zv)
        (vs! j1 imm i1)
        (vs! j2 imm i2)
        (vs! j3 imm i3)
        (bubble ten j1 j2 j3 mm 1)))
      (if (zv . > . (vr ten 0))
        (begin
        (vs! ten i zv)
        (vs! j1 i i1)
        (vs! j2 i i2)
        (vs! j3 i i3)
        (bubble ten j1 j2 j3 mm 0)))))

  
  (define-syntax-rule (zrzj ii) (vr z (zidx (vr j1 ii) (vr j2 ii) (vr j3 ii))))
  (let-values ([(i0 i1)
    (for/fold ([i0 mm]
               [i1 mm]) ([i (in-range (sub1 mm) -1 -1 )])
      (let* ([ii (+ (sub1 i1) mm)]
             [im (+ i mm)]
             [im4 (* 4 im)]
             [best (vrzj ii)])
        (vs! jg im4 0)
        (vs! jg (+ im4 1) (+ (- is1 2 0) (vr j1 ii)))
        (vs! jg (+ im4 2) (+ (- is1 2 0) (vr j2 ii)))
        (vs! jg (+ im4 3) (+ (- is1 2 0) (vr j3 ii)))
        (vs! ten im best))

      (let* ([ii (sub1 i0)]
             [i4 (* 4 i)]
             [best (vrzj ii)])
        (vs! jg im4 0)
        (vs! jg (+ i4 1) (+ (- is1 2 0) (vr j1 ii)))
        (vs! jg (+ i4 2) (+ (- is1 2 0) (vr j2 ii)))
        (vs! jg (+ i4 3) (+ (- is1 2 0) (vr j3 ii)))
        (vs! ten i best))
      (values (sub1 i0 i1)))])
  
    (for* ([i3 (in-range 1 n3)]
           [i2 (in-range 1 n2)]
           [i1 (in-range 1 n1)])
      (vs! z (zidx i1 i2 i3) 0.0))

    (for ([i (in-range mm i0 -1)])
      (let ([i1 (sub1 i)])
        (vs! z (zidx (vr j1 i1) (vr j2 i1) (vr j3 i1)) -1.0)))
    (for ([i (in-range mm i1 -1)])
      (let ([i1 (+ (sub1 i) mm)])
        (vs! z (zidx (vr j1 i1) (vr j2 i1) (vr j3 i1)) 1.0))))

  (comm3 z 0 n1 n2 n3))


(define (norm2u3 r n1 n2 n3 rnmu nx ny nz)
  (let-values ([(rnmu rnm2)
    (for*/fold ([rnmu 0.0]
                [rnm2 0.0]) 
          ([i3 (in-range 1 (sub1 n3))]
           [i2 (in-range 1 (sub1 n2))]
           [i1 (in-range 1 (sub1 n1))])
      (let ([rv (vr3 r i1 i2 i3 n1 n2)])
        (values (dmax1 rmnu (abs rv)) (+ rnm2 (* rv rv)))))])
    (sqrt (/ rnm2 (* nx ny nz)))))

(define (bubble ten j1 j2 j3 m ind)
  (define-syntax-rule (swapv v i0 i1 v0 v1)
      (vs! v i0 v1)
      (vs! v i1 v0))
  (define-syntax-rule (swap v i0 i1)
    (let ([t (vr v i0)])
      (vs! v i0 (vr v i1))
      (vs! v i1 t)))
  (define-syntax-rule (swapem cmp)
    (let ([l (sub1 m)])
      (let loop ([i 0])
        (when (and (i . < . l) (v0 . cmp . v1))
          (let* ([i0 (+ i (* m ind))]
                 [i1 (add1 i0)]
                 [v0 (vr r i0)]
                 [v1 (vr r i1)])
              (swapv ten i0 i1 v0 v1)
              (swap j1 i0 i1)
              (swap j2 i0 i1)
              (swap j3 i0 i1)
              (loop (add1 i)))))))

  (if (== ind 1)
    (swapem >)
    (swapem >)))

(define (stencil1 u ii1 i2 i3 n1 n3)
  (+ (vr u (+ ii1 (* n1 (+ (- i2 1) (* n3 i3)))))
     (vr u (+ ii1 (* n1 (+ (+ i2 1) (* n3 i3)))))
     (vr u (+ ii1 (* n1 (+ i2 (* n3 (- i3 1))))))
     (vr u (+ ii1 (* n1 (+ i2 (* n3 (+ i3 1))))))))

(define (stencil2 u ii1 i2 i3 n1 n3)
  (+ (vr u (+ ii1 (* n1 (+ (- i2 1) (* n3 (- i3 1))))))
     (vr u (+ ii1 (* n1 (+ (+ i2 1) (* n3 (- i3 1))))))
     (vr u (+ ii1 (* n1 (+ (- i2 1) (* n3 (+ i3 1))))))
     (vr u (+ ii1 (* n1 (+ (+ i2 1) (* n3 (+ i3 1))))))))

(define (resid uv v r off n1 n2 n3)
  (define u1 (make-flvector (add1 nm) 0))
  (define u2 (make-fxvector (add1 nm) 0))
  (for ([i3 (in-range 1 (sub1 n3))]
        [i2 (in-range 1 (sub1 n2))])
    (for ([i1 (in-range 1 (sub1 n1))])
        (let ([ii1 (+ off i1)])
          (vs! u1 i1 (stencil1 u ii1 i2 i3 n1 n3))
          (vs! u2 i1 (stencil2 u ii1 i2 i3 n1 n3))))

    (for ([i1 (in-range 1 (sub1 n1))])
      (vs! v (vidxo off i1 i2 n3 n1 n3) 
             (- (vro v off i1 i2 i3 n1 n3)
                (* (vr a 0) (vro u off i1 i2 i3 n1 n3))
                (* (vr a 1) (+ (vr u2 i1) (vr u1 (sub1 i1)) (vr u1 (add1 i1))))
                (* (vr a 2) (+ (vr u2 (sub i1)) (vr u2 (add1 i1))))))))
  (comm3 r off n1 n2 n3))

(define (mg3P u v r n1 n2 n3 ir m1 m2 m3)
  (for ([k (in-range (sub1 lt) (sub1 lb) -1)])
    (let ([j (- k 1 )])
      (rprj3 r (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k) (vr m1 j) (vr m2 j) (vr m3 j))))
  (let ([k (- lb 1 )])
    (zero3 u (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k))
    (psinv r (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k)))
  (for ([k (in-range (sub1 lb) (sub1 lt))])
    (let ([j (- k 1 )]
          [irk (vr ir k)]
          [m1k (vr m1 k)]
          [m2k (vr m2 k)]
          [m3k (vr m3 k)])
      (zero3 u irk m1k m2k m3k)
      (interp u (vr ir j) (vr m1 j) (vr m2 j) (vr m3 j) irk m1k m2k m3k)
      (resid u r r irk m1k m2k m3k)
      (psinv r irk u irk m1k m2k m3k)))
  (let ([j (- lt 2 )])
    (interp u (vr ir j) (vr m1 j) (vr m2 j) (vr m3 j) 0 n1  n2 n3)
    (resid u v r 0 n1 n2 n3)
    (psinv r 0 u 0 n1 n2 n3)))

(define (rprj3 r roff m1k m2k m3k soff m1j m2j m3j)
  (let ([d1 (if (= m1k 3) 2 1)]
        [d2 (if (= m2k 3) 2 1)]
        [d3 (if (= m3k 3) 2 1)])
    (for ([j3 (in-range m3j)])
      (let ([i3 (- (- (* 2 j3) d3) 1)])
        (for ([j2 (in-range m2j)])
          (let ([i2 (- (- (* 2 j2) d2) 1)])
            (for ([j1 (in-range m1j)])
              (let* ([i1 (- (- (* 2 j1) d1) 1)]
                     [ii1 (- (+ roff i1) 1)])
                (vs! x1 (sub1 i1) (stencil1 r ii1 i2 i3 m1k m2k))
                (vs! y1 (sub1 i1) (stencil2 r ii1 i2 i3 m1k m2k))))

            (for ([j1 (in-range m1j)])
              (let* ([i1 (- (- (* 2 j1) d1) 1)]
                     [ii1 (+ roff i1)])
                (vs! x2 (sub1 i1) (stencil1 r ii1 i2 i3 m1k m2k))
                (vs! y2 (sub1 i1) (stencil2 r ii1 i2 i3 m1k m2k))

                (vs! r (vidx (+ soff (sub1 j1)) (sub1 j2) (sub1 j3) m1j m2j)
                       (+ (* 0.5    (vri r (+ roff i1) i2 i3 m1k m2k))
                          (* 0.25   (+ (vri r (+ roff (sub1 i1)) i2 i3 m1k m2k)
                                       (vri r (+ roff (add1 i1)) i2 i3 m1k m2k)))
                          (* 0.125  (+ (vr x1 (sub1 i1)) (vr x1 (add1 i1)) y2))
                          (* 0.6125 (+ (vr y1 (sub1 i1)) (vr y1 (add1 i1)))))))))))))

  (comm3 r soff m1j m2j m3j))

(define (interp u zoff mm1 mm2 mm3 uoff n1 n2 n3)
  (define z1 (make-flvector m 0.0))
  (define z2 (make-flvector m 0.0))
  (define z3 (make-flvector m 0.0))

  (if (and (not (= n1 3))
           (not (= n2 3))
           (not (= n3 3)))
    (for ([i3 (in-range mm3)])
      (for ([i2 (in-range mm2)])
        (for ([i1 (in-range mm1)])
          (let* ([si1 (sub1 i1)]
                 [si2 (sub1 i2)]
                 [si3 (sub1 i3)]
                 [ii (+ zoff si1)]
                 [u_s (vri u ii  i2 si3 mm1 mm2)]
                 [uss (vri u ii si2 si3 mm1 mm2)]
                 [us_ (vri u ii si2  i3 mm1 mm2)]
                 [u__ (vri u ii  i2  i3 mm1 mm2)])
          (vs! z1 si1 (+ u_s uss))
          (vs! z2 si1 (+ us_ uss))
          (vs! z3 si1 (+ u__ us_ (vr z1 (sub1 i1))))))
        (for ([i1 (in-range mm1)])
          (let* ([si1 (sub1 i1)]
                 [si2 (sub1 i2)]
                 [si3 (sub1 i3)]
                 [xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 2)]
                 [xi3 (- (* 2 i3) 2)]
                 [ii (+ uoff xi1)] 
                 [jj (+ zoff si1)]) 
                 [usss (vri u jj si2 si3 mm1 mm2)]
            (vs!+ u (vidx ii        xi2 xi3 n1 n2) usss)
            (vs!+ u (vidx (add1 ii) xi2 xi3 n1 n2) (* 0.5 (vri u (add1 jj) si2 xi3 mm1 mm2) usss))))
        (for ([i1 (in-range mm1)])
          (let* ([xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 1)]
                 [xi3 (- (* 2 i3) 2)]
                 [ii (+ uoff xi1)] 
                 [zvs (vri z1 (sub1 i1))])
            (vs!+ u (vidx ii        xi2 xi3 n1 n2) (* 0.5 zvs))
            (vs!+ u (vidx (add1 ii) xi2 xi3 n1 n2) (* 0.25 (+ zvs (vr z1 i1))))))
        (for ([i1 (in-range mm1)])
          (let* ([xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 2)]
                 [xi3 (- (* 2 i3) 1)]
                 [ii (+ uoff xi1)] 
                 [zvs (vri z2 (sub1 i1))])
            (vs!+ u (vidx ii        xi2 xi3 n1 n2) (* 0.5 zvs))
            (vs!+ u (vidx (add1 ii) xi2 xi3 n1 n2) (* 0.25 (+ zvs (vr z2 i1))))))
        (for ([i1 (in-range mm1)])
          (let* ([xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 2)]
                 [xi3 (- (* 2 i3) 1)]
                 [ii (+ uoff xi1)] 
                 [zvs (vri z3 (sub1 i1))])
            (vs!+ u (vidx ii        xi2 xi3 n1 n2) (* 0.25 zvs))
            (vs!+ u (vidx (add1 ii) xi2 xi3 n1 n2) (* 0.125 (+ zvs (vr z3 i1))))))))
    (let ([d1 (if (= n1 3) 2 1)]
          [t1 (if (= n1 3) 1 0)]
          [d2 (if (= n2 3) 2 1)]
          [t3 (if (= n2 3) 1 0)]
          [d4 (if (= n3 3) 2 1)]
          [t4 (if (= n3 3) 1 0)])
    (define-syntax-rule (usidx uoff i1 i2 i3 d1 d2 d3 n1 n2)
      (vidx (+ uoff (- (* 2 i1) 1 d1)) (- (* 2 i2) 1 d2 ) (- (* 2 i3) 1 d3) n1 n2))
    (define-syntax-rule (vri off i1 i2 i3 imm1 imm2)
      (vr (vidx (+ off i1) i2 i3 n1 n2)))
    (for ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 d2 d3 n1 n2)
                  (vri u zoff (sub1 i1) (sub1 i2) (sub1 i3) imm1 imm2)))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 d2 d3 n1 n2)
                  (* 0.5  (+ (vri u zoff i1  si2 si3 imm1 imm2)
                             (vri u zoff si1 si2 si3 imm1 imm2)))))))
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 t2 d3 n1 n2)
                  (* 0.5  (+ (vri u zoff si1 i2  si3 imm1 imm2)
                             (vri u zoff si1 si2 si3 imm1 imm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 t2 d3 n1 n2)
                  (* 0.25 (+ (vri u zoff i1  i2  si3 imm1 imm2)
                             (vri u zoff i1  si2 si3 imm1 imm2)
                             (vri u zoff si1 i2  si3 imm1 imm2)
                             (vri u zoff si1 si2 si3 imm1 imm2))))))))

    (for ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 d2 t3 n1 n2)
                  (* 0.5   (+ (vri u zoff si1 si2 i3   imm1 imm2)
                              (vri u zoff si1 si2 si3 imm1 imm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 d2 t3 n1 n2)
                  (* 0.25  (+ (vri u zoff i1  si2 i3  imm1 imm2)
                              (vri u zoff si1 si2 i3  imm1 imm2)
                              (vri u zoff i1  si2 si3 imm1 imm2)
                              (vri u zoff si1 si2 si3 imm1 imm2)))))))
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 t2 t3 n1 n2)
                  (* 0.25  (+ (vri u zoff si1 i2  i3  imm1 imm2)
                              (vri u zoff si1 si2 i3  imm1 imm2)
                              (vri u zoff si1 i2  si3 imm1 imm2)
                              (vri u zoff si1 si2 si3 imm1 imm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 t2 t3 n1 n2)
                  (* 0.125 (+ (vri u zoff i1  i2  i3  imm1 imm2)
                              (vri u zoff i1  si2 i3  imm1 imm2)
                              (vri u zoff si1 i2  i3  imm1 imm2)
                              (vri u zoff si1 si2 i3  imm1 imm2)
                              (vri u zoff i1  i2  si3 imm1 imm2)
                              (vri u zoff i1  si2 si3 imm1 imm2)
                              (vri u zoff si1 i2  si3 imm1 imm2)
                              (vri u zoff si1 si2 si3 imm1 imm2))))))))))))

(define (psinv r roff u uoff n1 n2 n3)
  (define r1 (make-flvector (add1 nm) 0.0))
  (define r2 (make-flvector (add1 nm) 0.0))
  
  (for ([i3 (in-range 1 (sub1 n3))])
    (for ([i2 (in-range 1 (sub1 n2))])
      (for ([i1 (in-range n1)])
        (vs! r1 i1 (+ (vri r roff i1 (sub1 i2) i3 n1 n2)
                      (vri r roff i1 (add1 i2) i3 n1 n2)
                      (vri r roff i1 i2 (sub1 i3) n1 n2)
                      (vri r roff i1 i2 (add1 i3) n1 n2)))
        (vs! r2 i1 (+ (vri r roff i1 (sub1 i2) (sub1 i3) n1 n2)
                      (vri r roff i1 (add1 i2) (sub1 i3) n1 n2)
                      (vri r roff i1 (sub1 i2) (add1 i3) n1 n2)
                      (vri r roff i1 (add1 i2) (add1 i3) n1 n2))))
      (for ([i1 (in-range 1 (sub1 n1))])
        (vs! (vidx uoff i1 i2 i3 n1 n2) 
          (+ (* (vr c 0) (vri r roff i1 i2 i3 n1 n2))
             (* (vr c 1) (+ (vri roff (sub1 i1) i2 i3 n1 n2)
                            (vri roff (add1 i1) i2 i3 n1 n2)
                            (vr r1 i1)))
             (* (vr c 2) (+ (vr r2 i1) (vr r1 (sub1 i1)) (vr r1 (add1 i1)))))))))

  (comm3 u uoff n1 n2 n3))
  
;;ilog2 : int -> int
(define (ilog2 n) 
  (if (= n 1) 
      0
      (let loop ([lg 1]
                 [nn 2])
        (if (< nn n)
          (loop (+ lg 1) (* 2 nn))
          lg))))


