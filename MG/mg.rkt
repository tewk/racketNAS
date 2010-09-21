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
        (let-values ([(is1 is2 is3 ie1 ie2 ie3 n1 n2 n3) (setup lt ir nx ny nz m1 m2 m3 )])
      


      (zero3 u 0 n1 n2 n3)
      (zran3 v n1 n2 n3 (vr nx (sub1 lt)) (vr ny (sub1 lt)) is1 is2 is3 ie1 ie2 ie3)

      (resid a u v r 0 n1 n2 n3 nm)

      (for ([it (in-range 1 (add1 nit))])
        (mg3P c a u v r n1 n2 n3 lt)
        (resid a u v r 0 n1 n2 n3 nm))

      (timer-stop 1)


      (let* ([lt1 (sub1 lt)]
             [verified (verify CLASS 
              (norm2u3 r n1 n2 n3 void (vr nx lt1) (vr ny lt1) (vr nz lt1)))])
        (print-banner "MG" args) 
        (printf "Size = ~a X ~a X ~a niter = ~a~n" (vr nx lt1) (vr ny lt1) (vr nz lt1)  nit) 
        (if serial 
            (printf "SERIAL~n")
            (printf "PARALLEL~n"))
        (if verified 
            (printf "Verification succeeded~n") 
            (printf "Verification failed~n"))
        (let* ([time (/ (read-timer 1) 1000)]
               [results (new-BMResults bmname CLASS (vr nx lt1) (vr ny lt1) (vr nz lt1)  nit time 
                                       (get-mflops time nit n1 n2 n3)
                                       "floating point" 
                                       (if verified 1 0)
                                       serial 
                                       num-threads 
                                       -1)]) 
            (print-results results) 
            (when #f (print-timers)))))))))

(define-syntax-rule (vidx3 i1 i2 i3 n1 n2) (+ i1 (* n1 (+ i2 (* n2 i3)))))
(define-syntax-rule (vr3 v i1 i2 i3 n1 n2) (vr v (vidx3 i1 i2 i3 n1 n2)))
(define-syntax-rule (vidx off i1 i2 i3 n1 n2) (+ off (vidx3 i1 i2 i3 n1 n2)))
(define-syntax-rule (vro3 v off i1 i2 i3 n1 n2) (vr v (+ off (vidx3 i1 i2 i3 n1 n2))))

(define (verify class rnm2)
  (define verify-value 
    (case class
      [(#\S) 0.530770700573E-4]
      [(#\W) 0.250391406439E-17]
      [(#\A) 0.2433365309E-5]
      [(#\B) 0.180056440132E-5]
      [(#\C) 0.570674826298E-6]))

  (printf " L2 Norm is ~a\n" rnm2)

  (let ([deviation (abs (- verify-value rnm2))])
    (if (deviation . < . 1.0E-8)
      (begin
        (printf " Deviation is   ~a\n" deviation)
        #t)
      (begin
        (printf " The correct L2 Norm is ~n\n" verify-value)
        #f))))

(define (get-mflops total-time niter n1 n2 n3)
  (if (not (= total-time 0.0))
    (* (* 58.0 n1 n2 n3)
       (/ niter (* total-time 1000000.0)))
    0.0))

(define (get-input-pars maxlevel)
  (define fn "mg.input")
  (if (file-exists? fn)
    (match (call-with-input-file fn read)
      [(list lt lnx lny lnz nit)
        (when (lt . > . maxlevel)
          (printf "lt=~a Maximum allowable=~a\n" lt maxlevel)
          (exit 0))
        (values nit lt lnx lnz)]
      [else 
        (printf "Error reading from file mg.input\n")
        (exit 0)])
    (printf "No input file mg.input, Using compiled defaults\n")))

(define (setup lt ir nx ny nz m1 m2 m3 )
  (define lt1 (sub1 lt))

  (for ([k (in-range (- lt 2) -1 -1)])
    (let ([k1 (add1 k)])
        (vs! nx k (/ (vr nx k1) 2))
        (vs! ny k (/ (vr ny k1) 2))
        (vs! nz k (/ (vr nz k1) 2))
        (printf "~a ~a ~a ~a\n" k (vr nx k) (vr ny k) (vr nz k))
    ))
    
  (for ([k (in-range (- lt 1) -1 -1)])
    (let ([k1 (add1 k)])
        (vs! m1 k (+ 2 (vr nx k1)))
        (vs! m2 k (+ 2 (vr ny k1))) 
        (vs! m3 k (+ 2 (vr nz k1)))
        (printf "~a ~a ~a ~a\n" k (vr m1 k) (vr m2 k) (vr m3 k))
    ))

  (for ([j (in-range (- lt 2) -1 -1)])
    (let ([j1 (add1 j)])
    (vs! ir j (+ (vr ir j1) (* (vr m1 j1) (vr m2 j1) (vr m3 j1))))))

  (values 2 2 2 (add1 (vr nx lt1)) (add1 (vr ny lt1)) (add1 (vr nz lt1))
                (+ 2 (vr nx lt1)) (+ 2 (vr ny lt1)) (+ 2 (vr nz lt1))))

(define (zero3 z off n1 n2 n3)
  (for* ([i3 (in-range n3)]
         [i2 (in-range n2)]
         [i1 (in-range n1)])
    (vs! z (+ off i1 (* n1 (+ i2 (* n2 i3)))) 0.0)))

(define (zran3 z n1 n2 n3 nx ny is1 is2 is3 ie1 ie2 ie3)
  (define mm 10)
  (define j1 (make-fxvector (* 2 mm) 0))
  (define j2 (make-fxvector (* 2 mm) 0))
  (define j3 (make-fxvector (* 2 mm) 0))
  (define jg (make-fxvector (* 4 2 mm) 0))
  (define ten (make-flvector (* 2 mm) 0.0))
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
    (for/fold ([x0 (randlc/a 3.14159265 ai)])  ([i3 (in-range 2 (add1 e3))])
      (for/fold ([x1 x0]) ([i2 (in-range 2 (add1 e2))])
        (let ([xx x1])
          (vranlc d1 xx a z (vidx3 1 (sub1 i2) (sub1 i3) n1 n2))
          (randlc/a x1 a1)))
      (randlc/a x0 a2)))

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
      (when (zv . > . (vr ten mm))
        (begin
        (vs! ten mm zv)
        (vs! j1 mm i1)
        (vs! j2 mm i2)
        (vs! j3 mm i3)
        (bubble ten j1 j2 j3 mm 1)))
      (when (zv . > . (vr ten 0))
        (begin
        (vs! ten 0 zv)
        (vs! j1 0 i1)
        (vs! j2 0 i2)
        (vs! j3 0 i3)
        (bubble ten j1 j2 j3 mm 0)))))

  
  (define-syntax-rule (vrzj ii) (vr z (zidx (vr j1 ii) (vr j2 ii) (vr j3 ii))))
  (let-values ([(i0 i1)
    (for/fold ([i0 mm]
               [i1 mm]) ([i (in-range (sub1 mm) -1 -1 )])
      (let* ([ii (+ (sub1 i1) mm)]
             [im (+ i mm)]
             [im4 (* 4 im)]
             [best (vrzj ii)])
        (vs! jg im4 0)
        (vs! jg (+ im4 1) (+ (- is1 2) (vr j1 ii)))
        (vs! jg (+ im4 2) (+ (- is1 2) (vr j2 ii)))
        (vs! jg (+ im4 3) (+ (- is1 2) (vr j3 ii)))
        (vs! ten im best))

      (let* ([ii (sub1 i0)]
             [i4 (* 4 i)]
             [best (vrzj ii)])
        (vs! jg i4 0)
        (vs! jg (+ i4 1) (+ (- is1 2) (vr j1 ii)))
        (vs! jg (+ i4 2) (+ (- is1 2) (vr j2 ii)))
        (vs! jg (+ i4 3) (+ (- is1 2) (vr j3 ii)))
        (vs! ten i best))
      (values (sub1 i0) (sub1 i1)))])
  
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
        (values (max rnmu (abs rv)) (+ rnm2 (* rv rv)))))])
    (sqrt (/ rnm2 (* nx ny nz)))))

(define (bubble ten j1 j2 j3 m ind)
  (define-syntax-rule (swapv v i0 i1 v0 v1)
    (begin
      (vs! v i0 v1)
      (vs! v i1 v0)))
  (define-syntax-rule (swap v i0 i1)
    (let ([t (vr v i0)])
      (vs! v i0 (vr v i1))
      (vs! v i1 t)))
  (define-syntax-rule (swapem cmp)
    (let ([l (sub1 m)])
      (let loop ([i 0])
        (when (i . < . l)
          (let* ([i0 (+ i (* m ind))]
                 [i1 (add1 i0)]
                 [v0 (vr ten i0)]
                 [v1 (vr ten i1)])
            (when (v0 . cmp . v1)
              (swapv ten i0 i1 v0 v1)
              (swap j1 i0 i1)
              (swap j2 i0 i1)
              (swap j3 i0 i1)
              (loop (add1 i))))))))

  (if (= ind 1)
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

(define (resid a u v r off n1 n2 n3 nm)
  (define u1 (make-flvector (add1 nm) 0))
  (define u2 (make-fxvector (add1 nm) 0))
  (for ([i3 (in-range 1 (sub1 n3))]
        [i2 (in-range 1 (sub1 n2))])
    (for ([i1 (in-range 1 (sub1 n1))])
        (let ([ii1 (+ off i1)])
          (vs! u1 i1 (stencil1 u ii1 i2 i3 n1 n3))
          (vs! u2 i1 (stencil2 u ii1 i2 i3 n1 n3))))

    (for ([i1 (in-range 1 (sub1 n1))])
      (vs! v (vidx off i1 i2 n3 n1 n3) 
             (- (vro3 v off i1 i2 i3 n1 n3)
                (* (vr a 0) (vro3 u off i1 i2 i3 n1 n3))
                (* (vr a 1) (+ (vr u2 i1) (vr u1 (sub1 i1)) (vr u1 (add1 i1))))
                (* (vr a 2) (+ (vr u2 (sub1 i1)) (vr u2 (add1 i1))))))))
  (comm3 r off n1 n2 n3))

(define (mg3P c a u v r n1 n2 n3 ir m1 m2 m3 nm lt)
  (define lb 1)
  (for ([k (in-range (sub1 lt) (sub1 lb) -1)])
    (let ([j (- k 1 )])
      (rprj3 r (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k) (vr m1 j) (vr m2 j) (vr m3 j) nm)))
  (let ([k (- lb 1 )])
    (zero3 u (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k))
    (psinv c r (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k) nm))
  (for ([k (in-range (sub1 lb) (sub1 lt))])
    (let ([j (- k 1 )]
          [irk (vr ir k)]
          [m1k (vr m1 k)]
          [m2k (vr m2 k)]
          [m3k (vr m3 k)])
      (zero3 u irk m1k m2k m3k)
      (interp u (vr ir j) (vr m1 j) (vr m2 j) (vr m3 j) irk m1k m2k m3k)
      (resid a u r r irk m1k m2k m3k nm)
      (psinv c r irk u irk m1k m2k m3k nm)))
  (let ([j (- lt 2 )])
    (interp u (vr ir j) (vr m1 j) (vr m2 j) (vr m3 j) 0 n1  n2 n3)
    (resid a u v r 0 n1 n2 n3 nm)
    (psinv c r 0 u 0 n1 n2 n3 nm)))

(define (rprj3 r roff m1k m2k m3k soff m1j m2j m3j nm)
  (define x1 (make-flvector nm 0.0))
  (define y1 (make-flvector nm 0.0))
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
                     [ii1 (+ roff i1)]
                     [x2 (stencil1 r ii1 i2 i3 m1k m2k)]
                     [y2 (stencil2 r ii1 i2 i3 m1k m2k)])
                (vs! r (vidx soff (sub1 j1) (sub1 j2) (sub1 j3) m1j m2j)
                       (+ (* 0.5    (vro3 r roff i1 i2 i3 m1k m2k))
                          (* 0.25   (+ (vro3 r roff (sub1 i1) i2 i3 m1k m2k)
                                       (vro3 r roff (add1 i1) i2 i3 m1k m2k))
                                       x2)
                          (* 0.125  (+ (vr x1 (sub1 i1)) (vr x1 (add1 i1)) y2))
                          (* 0.6125 (+ (vr y1 (sub1 i1)) (vr y1 (add1 i1)))))))))))))

  (comm3 r soff m1j m2j m3j))

(define (interp u zoff mm1 mm2 mm3 uoff n1 n2 n3)
  (define-syntax-rule (vs!+ u idx v) (vs! u idx (+ (vr u idx) v)))
  (define m 535)
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
                 [u_s (vr3 u ii  i2 si3 mm1 mm2)]
                 [uss (vr3 u ii si2 si3 mm1 mm2)]
                 [us_ (vr3 u ii si2  i3 mm1 mm2)]
                 [u__ (vr3 u ii  i2  i3 mm1 mm2)])
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
                 [jj (+ zoff si1)] 
                 [usss (vr3 u jj si2 si3 mm1 mm2)])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) usss)
            (vs!+ u (vidx3 (add1 ii) xi2 xi3 n1 n2) (* 0.5 (vr3 u (add1 jj) si2 xi3 mm1 mm2) usss))))
        (for ([i1 (in-range mm1)])
          (let* ([xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 1)]
                 [xi3 (- (* 2 i3) 2)]
                 [ii (+ uoff xi1)] 
                 [zvs (vr z1 (sub1 i1))])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) (* 0.5 zvs))
            (vs!+ u (vidx3 (add1 ii) xi2 xi3 n1 n2) (* 0.25 (+ zvs (vr z1 i1))))))
        (for ([i1 (in-range mm1)])
          (let* ([xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 2)]
                 [xi3 (- (* 2 i3) 1)]
                 [ii (+ uoff xi1)] 
                 [zvs (vr z2 (sub1 i1))])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) (* 0.5 zvs))
            (vs!+ u (vidx3 (add1 ii) xi2 xi3 n1 n2) (* 0.25 (+ zvs (vr z2 i1))))))
        (for ([i1 (in-range mm1)])
          (let* ([xi1 (- (* 2 i1) 2)]
                 [xi2 (- (* 2 i2) 2)]
                 [xi3 (- (* 2 i3) 1)]
                 [ii (+ uoff xi1)] 
                 [zvs (vr z3 (sub1 i1))])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) (* 0.25 zvs))
            (vs!+ u (vidx3 (add1 ii) xi2 xi3 n1 n2) (* 0.125 (+ zvs (vr z3 i1))))))))
    (let ([d1 (if (= n1 3) 2 1)]
          [t1 (if (= n1 3) 1 0)]
          [d2 (if (= n2 3) 2 1)]
          [t2 (if (= n2 3) 1 0)]
          [d3 (if (= n3 3) 2 1)]
          [t3 (if (= n3 3) 1 0)])
    (define-syntax-rule (u!idx uoff i1 i2 i3 d1 d2 d3 n1 n2)
      (vidx3 (+ uoff (- (* 2 i1) 1 d1)) (- (* 2 i2) 1 d2 ) (- (* 2 i3) 1 d3) n1 n2))
    (for ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 d2 d3 n1 n2)
                  (vro3 u zoff (sub1 i1) (sub1 i2) (sub1 i3) mm1 mm2)))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 d2 d3 n1 n2)
                  (* 0.5  (+ (vro3 u zoff i1  si2 si3 mm1 mm2)
                             (vro3 u zoff si1 si2 si3 mm1 mm2)))))))
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 t2 d3 n1 n2)
                  (* 0.5  (+ (vro3 u zoff si1 i2  si3 mm1 mm2)
                             (vro3 u zoff si1 si2 si3 mm1 mm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 t2 d3 n1 n2)
                  (* 0.25 (+ (vro3 u zoff i1  i2  si3 mm1 mm2)
                             (vro3 u zoff i1  si2 si3 mm1 mm2)
                             (vro3 u zoff si1 i2  si3 mm1 mm2)
                             (vro3 u zoff si1 si2 si3 mm1 mm2))))))))

    (for ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 d2 t3 n1 n2)
                  (* 0.5   (+ (vro3 u zoff si1 si2 i3   mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 d2 t3 n1 n2)
                  (* 0.25  (+ (vro3 u zoff i1  si2 i3  mm1 mm2)
                              (vro3 u zoff si1 si2 i3  mm1 mm2)
                              (vro3 u zoff i1  si2 si3 mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2)))))))
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 t2 t3 n1 n2)
                  (* 0.25  (+ (vro3 u zoff si1 i2  i3  mm1 mm2)
                              (vro3 u zoff si1 si2 i3  mm1 mm2)
                              (vro3 u zoff si1 i2  si3 mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (sub1 i1)] [si2 (sub1 i2)] [si3 (sub1 i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 t2 t3 n1 n2)
                  (* 0.125 (+ (vro3 u zoff i1  i2  i3  mm1 mm2)
                              (vro3 u zoff i1  si2 i3  mm1 mm2)
                              (vro3 u zoff si1 i2  i3  mm1 mm2)
                              (vro3 u zoff si1 si2 i3  mm1 mm2)
                              (vro3 u zoff i1  i2  si3 mm1 mm2)
                              (vro3 u zoff i1  si2 si3 mm1 mm2)
                              (vro3 u zoff si1 i2  si3 mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2))))))))))))

(define (psinv c r roff u uoff n1 n2 n3 nm)
  (define r1 (make-flvector (add1 nm) 0.0))
  (define r2 (make-flvector (add1 nm) 0.0))
  
  (for ([i3 (in-range 1 (sub1 n3))])
    (for ([i2 (in-range 1 (sub1 n2))])
      (for ([i1 (in-range n1)])
        (vs! r1 i1 (+ (vro3 r roff i1 (sub1 i2) i3 n1 n2)
                      (vro3 r roff i1 (add1 i2) i3 n1 n2)
                      (vro3 r roff i1 i2 (sub1 i3) n1 n2)
                      (vro3 r roff i1 i2 (add1 i3) n1 n2)))
        (vs! r2 i1 (+ (vro3 r roff i1 (sub1 i2) (sub1 i3) n1 n2)
                      (vro3 r roff i1 (add1 i2) (sub1 i3) n1 n2)
                      (vro3 r roff i1 (sub1 i2) (add1 i3) n1 n2)
                      (vro3 r roff i1 (add1 i2) (add1 i3) n1 n2))))
      (for ([i1 (in-range 1 (sub1 n1))])
        (vs! (vidx uoff i1 i2 i3 n1 n2) 
          (+ (* (vr c 0) (vro3 r roff i1 i2 i3 n1 n2))
             (* (vr c 1) (+ (vro3 r roff (sub1 i1) i2 i3 n1 n2)
                            (vro3 r roff (add1 i1) i2 i3 n1 n2)
                            (vr r1 i1)))
             (* (vr c 2) (+ (vr r2 i1) (vr r1 (sub1 i1)) (vr r1 (add1 i1)))))))))

  (comm3 u uoff n1 n2 n3))

(define (comm3 u off n1 n2 n3)
  (for* ([i3 (in-range 1 (sub1 n3))]
         [i2 (in-range 1 (sub1 n2))])
    (vs! u (vidx off 0 i2 i3 n1 n2) (vro3 u off (- n1 2) i2 i3 n1 n2))
    (vs! u (vidx off (sub1 n1) i2 i3 n1 n2) (vro3 u off 1 i2 i3 n1 n2)))

  (for* ([i3 (in-range 1 (sub1 n3))]
         [i1 (in-range n1)])
    (vs! u (vidx off i1 0 i3 n1 n2) (vro3 u off i1 (- n2 2) i3 n1 n2))
    (vs! u (vidx off i1 (sub1 n2) i3 n1 n2) (vro3 u off i1 1 i3 n1 n2)))

  (for* ([i2 (in-range n2)]
         [i1 (in-range n1)])
    (vs! u (vidx off i1 i2 0 n1 n2) (vro3 u off i1 i2 (- n3 2) n1 n2))
    (vs! u (vidx off i1 i2 (sub1 n3) n1 n2) (vro3 u off i1 i2 1 n1 n2))))
  
;;ilog2 : int -> int
(define (ilog2 n) 
  (if (= n 1) 
      0
      (let loop ([lg 1]
                 [nn 2])
        (if (< nn n)
          (loop (+ lg 1) (* 2 nn))
          lg))))

(define (print-timers) 0)

