#lang racket/base

(provide main)
  
(require "../bm-args.ss") 
(require "../bm-results.ss") 
(require "../rand-generator.ss")
(require "../timer.ss")
(require "../parallel-utils.ss")
(require racket/match)
(require racket/math)
;(require racket/place)
;(require racket/place-utils)
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

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values  problem_size problem_size problem_size (make-fxvector 3 12) 0.015 100)]
    [(#\W) (values  problem_size problem_size problem_size (make-fxvector 3 36) 0.0015 400)]
    [(#\A) (values  problem_size problem_size problem_size (make-fxvector 3 64) 0.0015 400)
    [(#\B) (values  problem_size problem_size problem_size (make-fxvector 3 102) 0.001 400)]
    [(#\C) (values  problem_size problem_size problem_size (make-fxvector 3 162) 0.00067 400)]
    [else (error "Unknown class")]))

 (define ce (flvector 
     2.0 1.0 2.0 2.0 5.0 
     0.0 0.0 2.0 2.0 4.0 
     0.0 0.0 0.0 0.0 3.0 
     4.0 0.0 0.0 0.0 2.0 
     5.0 1.0 0.0 0.0 0.1 
     3.0 2.0 2.0 2.0 0.4 
     0.5 3.0 3.0 3.0 0.3 
     0.02 0.01 0.04 0.03 0.05 
     0.01 0.03 0.03 0.05 0.04 
     0.03 0.02 0.05 0.04 0.03 
     0.5 0.4 0.3 0.2 0.1 
     0.4 0.3 0.5 0.1 0.3 
     0.3 0.5 0.4 0.3 0.2))

(define (main . argv) 
  (let ([args (parse-cmd-line-args argv "Conjugate Gradient")]) 
    (run-benchmark args)))

(define make-fxvector make-vector)

(define (run-benchmark args) 
  (define maxlevel 11)
  (let ([bmname "CG"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(IMAX JMAX KMAX grid_points dt_default niter_default) (get-class-size CLASS)])
    (let* (
          [isize1 5]
          [jsize1 (* 5 (add1 IMAX))]
          [ksize1 (* 5 (add1 IMAX) (add1 JMAX))]
          [jsize2 (add1 IMAX)]
          [ksize2 (* (add1 IMAX) (add1 JMAX))]
          [jsize4 5]
          [s1 (* ksize1 KMAX)]
          [s2 (* ksize2 KMAX)]
          [s3 (* 5 (add1 problem_size))]
          [u (make-shared-flvector s1 0.0)]
          [rhs (make-shared-flvector s1 0.0)]
          [forcing (make-shared-flvector s1 0.0)]
          [us (make-shared-flvector s2 0.0)]
          [vs (make-shared-flvector s2 0.0)]
          [ws (make-shared-flvector s2 0.0)]
          [qs (make-shared-flvector s2 0.0)]
          [rho-i (make-shared-flvector s2 0.0)]
          [speed (make-shared-flvector s2 0.0)]
          [square (make-shared-flvector s2 0.0)]
          [lhs (make-shared-flvector s3 0.0)]
          [lhsp (make-shared-flvector s3 0.0)]
          [lhsm (make-shared-flvector s3 0.0)]
          [cv (make-shared-flvector problem_size 0.0)]
          [rhon (make-shared-flvector problem_size 0.0)]
          [rhos (make-shared-flvector problem_size 0.0)]
          [rhoq (make-shared-flvector problem_size 0.0)]
          [cuf (make-shared-flvector problem_size 0.0)]
          [q (make-shared-flvector problem_size 0.0)]



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
        (if verified 
            (printf "Verification Successful~n") 
            (printf "Verification Failed~n"))
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

(define (initialize)
  (for* ([k (in-range (vr grid_points 2))]
         [j (in-range (vr grid_points 1))]
         [i (in-range (vr grid_points 0))])
    (let ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))])
      (vs! u (+ 0 idx) 1.0)
      (vs! u (+ 1 idx) 0.0)
      (vs! u (+ 2 idx) 0.0)
      (vs! u (+ 3 idx) 0.0)
      (vs! u (+ 4 idx) 1.0)))

  (let ([Pface (make-fl-vector (* 5 3 2) 0.0)])
    (for ([k (in-range (vr grid_points 2))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range (vr grid_points 1))])
          (let ([eta (* j dnym1)])
            (for ([i (in-range (vr grid_points 0))])
              (let ([xi (* i dnxm1)])
                (for ([ix (in-range 2)])
                  (exact_solution ix eta zeta Pface (+ 0 (* 0 5) (* ix 15))))
                (for ([ix (in-range 2)])
                  (exact_solution xi ix zeta Pface (+ 0 (* 1 5) (* ix 15))))
                (for ([ix (in-range 2)])
                  (exact_solution xi eta ix Pface (+ 0 (* 2 5) (* ix 15))))

                (for ([m (in-range 5)])
                  (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
                        [pxi   (+ (* xi   (vr Pface (+ m (* 0 5) (* 1 15))))
                                 (* (- 1.0   xi) (vr Pface (+ m (* 0 5) (* 0 15)))))]
                        [peta  (+ (* eta  (vr Pface (+ m (* 1 5) (* 1 15))))
                                 (* (- 1.0  eta) (vr Pface (+ m (* 1 5) (* 0 15)))))]
                        [pzeta (+ (* zeta (vr Pface (+ m (* 2 5) (* 1 15))))
                                 (* (- 1.0 zeta) (vr Pface (+ m (* 2 5) (* 0 15)))))])
                    (flvs! u idx  (+ (- (+ pxi peta pzeta) (* pxi peta) (* pxi pzeta) (* peta pzeta))
                                     (* pxi peta pzeta))))))))))))
  (let ([temp (make-flvector 5 0.0)]
        [temp2 (make-flvector 5 0.0)]
        [i2 (sub1 (vr grid_points 0))]
        [j2 (sub1 (vr grid_points 1))]
        [k2 (sub1 (vr grid_points 2))])
    (for ([k (in-range (vr grid_points 2))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range (vr grid_points 1))])
          (let ([eta (* j dnym1)])
            (exact_solution 0.0 eta zeta temp 0)
            (exact_solution 1.0 eta zeta temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* j jsize1) (* k ksize1))]
                     [idx2 (+ (* i2 isize1) idx)])
                (flvs! u idx (vr temp m))
                (flvs! u idx2 (vr temp2 m))))))
            
        (for ([i (in-range (vr grid_points 0))])
          (let ([xi (* i dnxm1)])
            (exact_solution xi 0.0 zeta temp 0)
            (exact_solution xi 1.0 zeta temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* i isize1) (* k ksize1))]
                     [idx2 (+ (* j2 jsize1) idx)])
                (flvs! u idx (vr temp m))
                (flvs! u idx2 (vr temp2 m))))))))
  
    (for ([j (in-range (vr grid_points 1))])
      (let ([eta (* j dnym1)])
        (for ([i (in-range (vr grid_points 0))])
          (let ([xi (* i dnxm1)])
            (exact_solution xi eta 0.0 temp 0)
            (exact_solution xi eta 1.0 temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* i isize1) (* j jsize1))]
                     [idx2 (+ (* k2 ksize1) idx)])
                (flvs! u idx (vr temp m))
                (flvs! u idx2 (vr temp2 m))))))))))

(define (add)
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
      (flvs!+ u idx (vr rhs idx)))))

(define (error-norm rms)
  (for ([m (in-range 5)]) (flvs! rms m 0.0))

  (let ([u-exact (make-fl-vector 5 0.0)])
    (for ([k (in-range (vr grid_points 2))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range (vr grid_points 1))])
          (let ([eta (* j dnym1)])
            (for ([i (in-range (vr grid_points 0))])
              (let ([xi (* i dnxm1)])
                (exact_solution xi eta zeta u-exact 0)
                (for ([m (in-range 5)])
                  (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
                         [add (- (vr u idx) (vr u-exact m))])
                    (flvs!+ rms m (* add add)))))))))))

  
  (for ([m (in-range 5)])
    (flvs! rms m (sqrt (for/fold ([x (flvr rms m)]) ([d (in-range 3)])
      (/ x (1 (vr grid_points d) 2)))))))

(define (error-norm rms)
  (for ([m (in-range 5)]) (flvs! rms m 0.0))

  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
           [add (vr rhs idx)])
      (flvs!+ rms m (* add add))))
  
  (for ([m (in-range 5)])
    (flvs! rms m (sqrt (for/fold ([x (flvr rms m)]) ([d (in-range 3)])
      (/ x (1 (vr grid_points d) 2)))))))

(define-syntax-rule (flvector-set!-all v val)
  (for ([i (in-range (flvector-length v))])
    (flvs! v i val)))

(define-syntax-rule (sqr v ) (* v v))

(define (exact_rhs)
#|
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
      (flvs!+ forcing idx 0.0)))
|#

  (flvector-set!-all forcing 0.0)

(define-syntax-rule (body ii giidx t_2 f1jsize3 dtemp jsize3 ue buf cuf q c1 c2 i j k hh 
  _con_0 _con_1 _con_2 __con3 __con4 __con5
  d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1)
          (for ([ii (in-range 0 (- (vr grid_points giidx) 1))])
              (let ([xi (* i dnxm1)]
                    [eta (* j dnym1)]
                    [zeta (* k dnzm1)])
                (exact_solution xi eta zeta dtemp 0)
                (let ([dtpp (/ 1.0 (vr dtemp 0))])
                  (for ([m (in-range 5)])
                    (let ([idx (+ ii (* m jsize3))]
                          [dtv (vr dtemp m)])
                    (flvs! ue idx dtv)
                    (flvs! buf idx (* dtpp dtv)))))
                (let* ([ij (+ ii jsize3)]
                       [i2j (+ ii (* 2 jsize3))]
                       [i3j (+ ii (* 3 jsize3))]
                       [bufij (vr buf ij)]
                       [bufi2j (vr buf i2j)]
                       [bufi3j (vr buf i3j)]
                       [ueij (vr ue ij)]
                       [uei2j (vr ue i2j)]
                       [uei3j (vr ue i3j)]
                       [bufij2 (sqr bufij)]
                       [bufi2j2 (sqr bufi2j)]
                       [bufi3j2 (sqr bufi3j)]
                (flvs! cuf ii (sqr (vr buf (+ ii (* hh jsize3)))))
                (flvs! buf ii (+ bufij2 bufi2j2 bufi3j2))
                (flvs! q ii (* 0.5 (+ (* bufij ueij) (* bufi2j uei2j) (* bufi3j uei3j))))))))

          (for ([ii (in-range (- (vr grid_points giidx) 2))])
            (let ([ip1 (add1 ii)]
                  [im1 (sub1 ii)]
                  [idx (+ (* i isize1) (* j jsize1) (* k ksize1))])

              (define-syntax-rule (vrit v s t) (vr v (+ s t)))
              (define-syntax-rule (cit c v t)
                (* c (+ (- vrit v ip1 t) (* 2.0 (vrit v ii t)) (vrit v im1 t))))
              (define-syntax-rule (dit d t) (cit d ue t))
              (define-syntax-rule (xit x t) (cit x buf t))

              (define-syntax-rule (t_2it l r) (- (* t_2 (-  l r))))
              (define-syntax-rule (t_2_mw f1jsize3 mjs3) (* (vrit ue s mjs3) (vrit buf s f1jsize3)))
              (define-syntax-rule (t_2_m i1)
                (t_2it (t_2_mw ip1 i1)
                       (t_2_2mw im1 i1))))

              (flvs!+ forcing idx 
                   (t_2it (vrit ue ip1 f1jsize3)
                          (vrit ue im1 f1jsize3))
                   (dit d_1t_1 0))

              (let ([x4jsize3 (* 4 jsize3)])
                (define-syntax-rule (tx2_1 s)
                  (+ (* (vrit ue s jsize3) (vrit buf s jsize)) (- (* c2 (vrit ue s x4jsize3)))))
                (flvs!+ forcing (+ idx 1) t_1 (xit _con_0 jsize3) (dit d_2t_1 jsize3)))

              (let ([x2jsize3 (* 2 jsize3)])
                (flvs!+ forcing (+ idx 2) t_2 (xit _con_1 x2jsize3) (dit d_3t_1 x2jsize3)))

              (let ([x3jsize3 (* 3 jsize3)])
                (flvs!+ forcing (+ idx 3) t_3 (xit _con_2 x3jsize3) (dit d_4t_1 x3jsize3)))
              
              (let ([x4jsize3 (* 4 jsize3) ])
                (define-syntax-rule (tx2_4 s)
                  (* (vrit buf s iijsize3) (- (* c1 (vrit ue s x4jsize3)) (*c2 (vr q s)))))
                (flvs!+ forcing (+ idx 4)
                     (tx2it (t_2_4 ip1) (t_2_4 im1))
                     (* 0.5 (cit __con3 buf 0))
                     (cit __con4 cuf 0)
                     (xit __con4 x4jsize3)
                     (dit d_5t_1 x4jsize3))))
;//---------------------------------------------------------------------
;//            Fourth-order dissipation
;//---------------------------------------------------------------------
          (for ([m (in-range 5)])
            (let* ([ii 1]
                   [imj (+ ii (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (* 5.0 (flvr ue imj)) 
                              (* -4.0 (flvr ue (+ imj 1))))
                              (flvr ue (+ imj 2))))))
            (let* ([ii 2]
                   [imj (+ ii (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (* -4.0 (flvr ue (+ imj 1)))
                              (*  6.0 (flvr ue imj))
                              (* -4.0 (flvr ue (+ imj 1)))
                              (flvr ue (+ imj 2))))))))
          (for* ([m (in-range 5)]
                 [ii (in-range 3 (- (vr grid_points giidx) 3))])
            (let* ([imj (+ ii (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (flvr ue (- imj 2)) 
                              (* -4.0 (flvr ue (- imj 1))))
                              (*  6.0 (flvr ue imj))
                              (* -4.0 (flvr ue (+ imj 1)))
                              (flvr ue (+ imj 2)))))))
          (for ([m (in-range 5)])
            (let* ([ii (- (vr grid_points giidx) 3)]
                   [imj (+ ii (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (flvr ue (- imj 2)) 
                              (* -4.0 (flvr ue (- imj 1)))
                              (*  6.0 (flvr ue imj))
                              (* -4.0 (flvr ue (+ imj 1))))))))
            (let* ([ii (- (vr grid_points giidx) 3)]
                   [imj (+ ii (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (flvr ue (- imj 2)) 
                              (* -4.0 (flvr ue (- imj 1)))
                              (*  5.0 (flvr ue imj)))))))))


  (let ([dtemp (make-fl-vector 5 0.0)])
    (for ([k (in-range 1 (- (vr grid_points 2) 2) )])
        (for ([j (in-range 1 (- (vr grid_points 1) 2))])
           (body2 i 0 tx2 jsize3 jsize3 ue buf cuf q c1 c2 i j k 1
             (t_2it (t_2_1 ip1) (tx2_1 im1))

              (t_2_mw (* 1 jsize3) (* 2 jsize3))
              (t_2_mw (* 1 jsize3) (* 3 jsize3))
              xxcon1 xxcon2 xxcon2 xxcon3 xxcon4 xxcon5
              dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)
)

        (for ([i (in-range 1 (- (vr grid_points 0) 2))])
            (body2 j 1 ty2 (* 2 jsize3) jsize3 ue buf cuf q c1 c2 i j k 2
              (t_2_mw (* 2 jsize3) (* 1 jsize3))

              (t_2_mw (* 2 jsize3) (* 3 jsize3))
              yycon2 yycon1 yycon2 yycon3 yycon4 yycon5
              dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)
))

    (for ([j (in-range 1 (- (vr grid_points 1) 2))])
        (for ([i (in-range 1 (- (vr grid_points 0) 2))])
            (body2 k 2 tz2 (* 3 jsize3) jsize3 ue buf cuf q c1 c2 i j k 3
              (t_2_mw (* 3 jsize3) (* 1 jsize3))
              (t_2_mw (* 3 jsize3) (* 2 jsize3))

              zzcon2 zzcon2 zzcon1 zzcon3 zzcon4 zzcon5
              dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)
)))

    (for* ([k (in-range 1 (- (vr grid_points 2) 2))]
           [j (in-range 1 (- (vr grid_points 1) 2))]
           [i (in-range 1 (- (vr grid_points 0) 2))]
           [m (in-range 5)])
      (flvs!* forcing (+ m (* i isize1) (* j jsize1) (* k ksize1)) -1.0))
)
#|
                (exact_solution xi eta zeta dtemp 0)
                (let ([dtpp (/ 1.0 (vr dtemp 0))])
                  (for ([m (in-range 5)])
                    (let ([idx (+ i (* m jsize3))]
                          [dtv (vr dtemp m)])
                    (flvs! ue idx dtv)
                    (flvs! buf idx (* dtpp dtv)))))
                (let* ([ij (+ i jsize3)]
                       [i2j (+ i (* 2 jsize3))]
                       [i3j (+ i (* 2 jsize3))]
                       [bufij (vr buf ij)]
                       [bufi2j (vr buf i2j)]
                       [bufi3j (vr buf i3j)]
                       [ueij (vr ue ij)]
                       [uei2j (vr ue i2j)]
                       [uei3j (vr ue i3j)]
                       [bufij2 (sqr bufij)]
                       [bufi2j2 (sqr bufi2j)]
                       [bufi3j2 (sqr bufi3j)]
                (flvs! cuf i bufij2)
                (flvs! buf i (+ bufij2 bufi2j2 bufi3j2))
                (flvs! q i (* 0.5 (+ (* bufij ueij) (* bufi2j uei2j) (* bufi3j uei3j)))))))))

          (for ([i (in-range (- (vr grid_points 0) 2))])
            (let ([ip1 (add1 i)]
                  [im1 (sub1 i)]
                  [idx (+ (* i isize1) (* j jsize1) (* k ksize1))])

              (define-syntax-rule (vrit v s ii) (vr v (+ s ii)))
              (define-syntax-rule (cit c v ii)
                (* c (+ (- vrit v ip1 ii) (* 2.0 (vrit v i ii)) (vrit v im1 ii))))
              (define-syntax-rule (dit d ii) (cit d ue ii))
              (define-syntax-rule (xit x ii) (cit x buf ii))

              (define-syntax-rule (tx2it l r) (- (* tx2 (-  l r))))
              (define-syntax-rule (tx2_2or3clause s i1) (* (vrit ue s i1) (vrit buf s jsize3)))
              (define-syntax-rule (tx2_2or3 i1)
                (tx2it (tx2_2or3clause ip1 i1)
                       (tx2_2or3clause im1 i1))))

              (flvs!+ forcing idx 
                   (tx2it (vrit ue ip1 jsize3)
                          (vrit ue im1 jsize3))
                   (dit dx1tx1 0))

              (let ([x4jsize3 (* 4 jsize3)])
                (define-syntax-rule (tx2_1 s)
                  (+ (* (vrit ue s jsize3) (vrit buf s jsize)) (- (* c2 (vrit ue s x4jsize3)))))
                (flvs!+ forcing (+ idx 1)
                   (tx2it (tx2_1 ip1) (tx2_1 im1))
                   (xit xxcon1 jsize3)
                   (dit dx2tx1 jsize3)))

              (let ([x2jsize3 (* 2 jsize3) ])
                (flvs!+ forcing (+ idx 2)
                     (tx2_2or3 x2jsize3)
                     (xit xxcon2 x2jsize3)
                     (dit dx3tx1 x2jsize3)))

              (let ([x3jsize3 (* 3 jsize3) ])
                (flvs!+ forcing (+ idx 3)
                     (tx2_2or3 x3jsize3)
                     (xit xxcon2 x3jsize3)
                     (dit dx4tx1 x3jsize3)))
              
              (let ([x4jsize3 (* 4 jsize3) ])
                (define-syntax-rule (tx2_4 s)
                  (* (vrit buf s jsize3) (- (* c1 (vrit ue s x4jsize3)) (*c2 (vr q s)))))
                (flvs!+ forcing (+ idx 4)
                     (tx2it (tx2_4 ip1) (tx2_4 im1))
                     (* 0.5 (cit xxcon3 buf 0))
                     (cit xxcon4 cuf 0)
                     (xit xxcon4 x4jsize3)
                     (dit dx5tx1 x4jsize3))))


;//---------------------------------------------------------------------
;//            Fourth-order dissipation                         
;//---------------------------------------------------------------------

          (for ([m (in-range 5)])
            (let* ([i 1]
                   [imj (+ i (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (* 5.0 (flvr ue imj)) 
                              (* -4.0 (flvr ue (+ imj 1))))
                              (flvr ue (+ imj 2))))))
            (let* ([i 2]
                   [imj (+ i (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (* -4.0 (flvr ue (+ imj 1)))
                              (*  6.0 (flvr ue imj))
                              (* -4.0 (flvr ue (+ imj 1)))
                              (flvr ue (+ imj 2))))))))
          (for* ([m (in-range 5)]
                 [i (in-range 3 (- (vr grid_points 0) 3))])
            (let* ([imj (+ i (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (flvr ue (- imj 2)) 
                              (* -4.0 (flvr ue (- imj 1))))
                              (*  6.0 (flvr ue imj))
                              (* -4.0 (flvr ue (+ imj 1)))
                              (flvr ue (+ imj 2)))))))
          (for ([m (in-range 5)])
            (let* ([i (- (vr grid_points 0) 3)]
                   [imj (+ i (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (flvr ue (- imj 2)) 
                              (* -4.0 (flvr ue (- imj 1)))
                              (*  6.0 (flvr ue imj))
                              (* -4.0 (flvr ue (+ imj 1))))))))
            (let* ([i (- (vr grid_points 0) 3)]
                   [imj (+ i (* m jsize3))]
                   [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!+ forcing idx 
                (- (* dssp (+ (flvr ue (- imj 2)) 
                              (* -4.0 (flvr ue (- imj 1)))
                              (*  5.0 (flvr ue imj))))))))
|#

(define (ninvr)
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [r1 (flvr rhs (+ 0 idx))]
           [r2 (flvr rhs (+ 1 idx))]
           [r3 (flvr rhs (+ 2 idx))]
           [r4 (flvr rhs (+ 3 idx))]
           [r5 (flvr rhs (+ 4 idx))]
           [t1 (* bt r3)]
           [t2 (* 0.5 (+ r4 r5))])
      (flvs! rhs (+ 0 idx) (- r2)
      (flvs! rhs (+ 1 idx) r1)
      (flvs! rhs (+ 2 idx) (* bt (- r4 r5)))
      (flvs! rhs (+ 3 idx) (- t2 t1))
      (flvs! rhs (+ 4 idx) (+ t1 t2)))))

(define (pinvr)
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [r1 (flvr rhs (+ 0 idx))]
           [r2 (flvr rhs (+ 1 idx))]
           [r3 (flvr rhs (+ 2 idx))]
           [r4 (flvr rhs (+ 3 idx))]
           [r5 (flvr rhs (+ 4 idx))]
           [t1 (* bt r3)]
           [t2 (* 0.5 (+ r4 r5))])
      (flvs! rhs (+ 0 idx) (* bt (- r4 r5)))
      (flvs! rhs (+ 1 idx) (- r3))
      (flvs! rhs (+ 2 idx) r2)
      (flvs! rhs (+ 3 idx) (- t2 t1))
      (flvs! rhs (+ 4 idx) (+ t1 t2)))))

(define compute_rhs)
  (for* ([k (in-range (vr grid_points 2))]
         [j (in-range (vr grid_points 1))]
         [i (in-range (vr grid_points 0))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [idx2 (+ i (* j jsize2) (* k ksize2))]
           [rho_inv (/ 1.0 idx)]
           [u1 (vs u (+ idx 1))]
           [u2 (vs u (+ idx 2))]
           [u3 (vs u (+ idx 3))]
           [u4 (vs u (+ idx 4))]
           [sq (* 0.5 (+ (sqr u1) (sqr u2) (sqr u3)) rho_inv)])
      (flvs! rho_i idx2 rho_inv)
      (flvs! us idx2 (* rho_inv u1))
      (flvs! vs idx2 (* rho_inv u2))
      (flvs! ws idx2 (* rho_inv u3))
      (flvs! squre idx2 sq)
      (flvs! qs idx2 (* rho_inv sq))
      (flvs! speed idx2 (sqrt (* c2c2 rho_inv (- u4 sq))))
  
      (for* ([m (in-range 5)])
        (flvs! rhs (+ m idx) (vr forcing (+ m idx))))))

(b us (+ i 1 (* j jsize1) (* k ksize1)) (+ i -1 (* j jsize1) (* k ksize1))

dx1tx1
(b vs (+ i (* (+ j 1) jsize1) (* k ksize1)) (+ i (* (- j 1) jsize1) (* k ksize1))
dy1ty1
(b ws (+ i (* j jsize1) (* (+ k 1) ksize1)) (+ i (* j jsize1) (* (- k 1) ksize1))
dz1tz1
(define-syntax-rule (b zs idxz+ idxz-
d_1t_1
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [idx2 (+ i (* j jsize1) (* k ksize1))]
           [idxi (+ i 1 (* j jsize1) (* k ksize1))]
           [idxj (+ i (* (+ j 1)  jsize1) (* k ksize1))]
           [idxk (+ i (* j jsize1) (* (+ k 1) ksize1))])
           [zijk (vr zs idx2)]
           [zp1 (vr zs idx2z+)]
           [zm1 (vr zs idx2z-)]
      (flvs!+ rhs idx
        (* d_1t_1 (+ (- (vr u idxz+) (* 2.0 (vr u idx))) (vr u idxz-)))
        (- (* t_2 (- (vr u (+ idxz+ 1)) (vr u (+ idxz- 1))))))

      (flvs!+ rhs (+ idx 1)
        (* d_2t_1 (+ (- (vr u (+ idxz+ 1)) (* 2.0 (vr u (+ idx 1)))) (vr u (+ idxz- 1))))
        (* t_2 (- (vr 
   
      (flvs!+ rhs (+ idx 2)
        (* d_3t_1 (+ (- (vr u (+ idxz+ 1)) (* 2.0 (vr u (+ idx 1)))) (vr u (+ idxz- 1))))
        (* t_2 (- (vr 
    
      (flvs!+ rhs (+ idx 3)
        (* d_4t_1 (+ (- (vr u (+ idxz+ 1)) (* 2.0 (vr u (+ idx 1)))) (vr u (+ idxz- 1))))
        (* t_2 (- (vr 

      (flvs!+ rhs (+ idx 4)
        (* d_5t_1 (+ (- (vr u (+ idxz+ 1)) (* 2.0 (vr u (+ idx 1)))) (vr u (+ idxz- 1))))
        (* __con3 (+ (- (vr u (+ idxz+ 1)) (* 2.0 (vr u (+ idx 1)))) (vr u (+ idxz- 1))))
        (* __con4 (+ (- (vr u (+ idxz+ 1)) (* 2.0 (vr u (+ idx 1)))) (vr u (+ idxz- 1))))
        (* t_2 (- (vr 

(define (txinvr)
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [idx2 (+ i (* j jsize2) (* k ksize2))]

           [ru1 (flvr rho_i idx2)]
           [uu (flvr us idx2)]
           [vv (flvr vs idx2)]
           [ww (flvr ws idx2)]
           [ac (flvr speed idx2)]
           [ac2inv (/ 1.0 (sqr ac))]

           [r1 (flvr rhs (+ 0 idx))]
           [r2 (flvr rhs (+ 1 idx))]
           [r3 (flvr rhs (+ 2 idx))]
           [r4 (flvr rhs (+ 3 idx))]
           [r5 (flvr rhs (+ 4 idx))]
           [t1 (* c2 ac2inv (+ (- (* (vr qs idx2) r1) (* uu r2) (* vv r3) (* ww r4)) r5))]
           [t2 (* bt ru1 (- (* uu r1) r2))]
           [t3 (* bt ru1 ac t1)])

      (flvs! rhs (+ 0 idx) (- r1 t1))
      (flvs! rhs (+ 1 idx) (- (* ru1 (- (* ww r1) r4))))
      (flvs! rhs (+ 2 idx) (- (* ru1 (- (* vv r1) r3))))
      (flvs! rhs (+ 3 idx) (- t2 t3))
      (flvs! rhs (+ 4 idx) (+ t1 t3)))))

(define (pinvr)
  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [idx2 (+ i (* j jsize2) (* k ksize2))]

           [xvel (flvr us idx2)]
           [uvel (flvr vs idx2)]
           [zvel (flvr ws idx2)]
           [ac (flvr speed idx2)]
           [ac2u (sqr ac)]

           [r1 (flvr rhs (+ 0 idx))]
           [r2 (flvr rhs (+ 1 idx))]
           [r3 (flvr rhs (+ 2 idx))]
           [r4 (flvr rhs (+ 3 idx))]
           [r5 (flvr rhs (+ 4 idx))]

           [uzik1 (vr i idx)]
           [btiz  (* bt usik1)]

           [t1 (* (/ btuz ac) (+ r4 r5))]
           [t2 (+ r3 t1)]
           [t3 (* btuz (- r4 r5))])
      (flvs! rhs (+ 0 idx) t2)
      (flvs! rhs (+ 1 idx) (+ (- (* uzik1 r2)) (* xvel t2)))
      (flvs! rhs (+ 2 idx) (+ (* uzik1 r1) (* yvel t2))) 
      (flvs! rhs (+ 3 idx) (+ (* zvel t2) t3))
      (flvs! rhs (+ 4 idx) (+ (* uzik1 (- (* yvel r1) (*xvel r2)))
                              (* (vr qs idx2) t2)
                              (* c2iv ac2u t1)
                              (* zvel t3))))))
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

(define (print-timers) 0)

