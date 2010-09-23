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

(define (compute_rhs)
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

(BODY i j k 1 us 
  (+ (* (+ i 1) isize1) (* j jsize1) (* k ksize1)) 
  (+ (* (- i 1) isize1) (* j jsize1) (* k ksize1))
  (+ (+ i 1) (* j jsize2) (* k ksize2)) 
  (+ (- i 1) (* j jsize2) (* k ksize2))
  dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
  (* xxcon2 con43) xxcon2 xxcon2 xxcon3 xxcon4 xxcon5
  tx2 t_2m 0 0)

(BODY i j k 2 vs 
  (+ (* i isize1) (* (+ j 1) jsize1) (* k ksize1)) 
  (+ (* i isize1) (* (- j 1) jsize1) (* k ksize1))
  (+ i (* (+ j 1) jsize2) (* k ksize2)) 
  (+ i (* (- j 1) jsize2) (* k ksize2))
  dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
  yycon2 (* yycon2 con43) yycon2 yycon3 yycon4 yycon5
  ty2 0 t_2m 0)

(BODY i j k 3 ws 
  (+ (* i isize1) (* j jsize1) (* (+ k 1) ksize1)) 
  (+ (* i iszie1) (* j jsize1) (* (- k 1) ksize1))
  (+ i (* j jsize2) (* (+ k 1) ksize2)) 
  (+ i (* j jsize2) (* (- k 1) ksize2))
  dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
  zzcon2 zzcon2 (* zzcon2 con43) zzcon3 zzcon4 zzcon5
  tz2 0 0 t_2m)

(define-syntax-rule (BODY i j k A zs 
  idxz+ idxz- idx2z+ idx2z-
  d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1
  __con21 __con22  __con23 __con3 __con4 __con5
  t_2 t_2m1 t_2m2 t_2m3)

  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))])
    (let* ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))]
           [idx2 (+ i (* j jsize1) (* k ksize1))])

      (define-syntax-rule (vrit V I1 I2) (vr V (+ I1 I2)))
      (define-syntaz-rule (cit C V A)
        (* C (+ (- (vrit V idxz+ A) (* 2.0 (vrit V idz A)) (vrit V idxz- A)))))
      (define-syntaz-rule (cit2 C V)
        (* C (+ (- (vr V idx2z+) (* 2.0 (vr V idx2)) (vr V idx2z-)))))
      (define-syntaz-rule (cit3 C F V)
        (* C (+ (- (F (vr V idx2z+)) (* 2.0 (F (vr V idx2))) (F (vr V idx2z-))))))
      (define-syntax-rule (t_2m)
        (* (+ (- (vr u (+ idxz+ 4)) (vr square idx2z+) (vr u (+ idxz- 4))) (vr square idx2z-)) c2))

      (define-syntax-rule (mid d__t_1 __con2X A t_2m_)
        (flvs!+ rhs (+ idx A)
          (cit d__t_1 u A)
          (cit2 __con2X zs)
          (- (* t_2 (+ (- (* (vr u (+ idxz+ A) (vr zs idx2+))) 
                          (* (vr u (+ idxz- A) (vr zs idx2-)))) t_2m_)))))
  
      (flvs!+ rhs idx
        (cit d_1t_1 u 0)
        (- (* t_2 (- (vr u (+ idxz+ A)) 
                     (vr u (+ idxz- A))))))

      (mid d_2t_1 __con21 1 t_2m1)
      (mid d_3t_1 __con22 2 t_2m2)
      (mid d_4t_1 __con23 3 t_2m3)

      (flvs!+ rhs (+ idx 4)
        (cit d_5t_1 u 4)
        (cit2 __con3 qs)
        (cit3 __con4 sqr zs)
        (* __con5 (+ (-
            (*     (vr u (+ idxz+ 4)) (vr rho_i idx2z+))
            (* 2.0 (vr u (+ idxz  4)) (vr rho_i idx2z )))
            (*     (vr u (+ idxz- 4)) (vr rho_i idx2z-))))
        (- (* t_2 (- (* (- (* c1 (vr u (+ idxz+ 4))) 
                           (* c2 (vr square idx2z+))) 
                        (vr zs idxz+))
                     (* (- (* c1 (vr u (+ idxz- 4))) 
                           (* c2 (vr square idx2z-))) 
                        (vr zs idxz-)))))))))

(define-syntax-rule (dissip1 ii i j k midx- midx midx+ midx+2)
  (let* ([ii 1]
         [idx (+ (* i isize1) (* j jsize1) (* k ksize1))])
    (for ([m (in-range 5)])
      (let ([midx (+ m idx)])
        (flvs! rhs midx (- (* dssp (+ (*  5.0 (flvr u midx))
                                      (* -4.0 (flvr u midx+)))
                                              (flvr u midx+2)))))))
  (let* ([ii 2]
         [idx (+ (* i isize1) (* j jsize1) (* k ksize1))])
    (for ([m (in-range 5)])
      (let ([midx (+ m idx)])
        (flvs! rhs midx (- (* dssp (+ (* -4.0 (flvr u midx-)) 
                                      (*  6.0 (flvr u midx))
                                      (* -4.0 (flvr u midx+))
                                              (flvr u midx+2)))))))))

(define-syntax-rule (dissip2 i j k)
    (for ([m (in-range 5)])
      (let ([midx (+ m idx)])
        (flvs! rhs midx (- (* dssp (+         (flvr u midx-2)
                                      (* -4.0 (flvr u midx-))
                                      (*  6.0 (flvr u midx))
                                      (*  4.0 (flvr u midx+))
                                              (flvr u midx+2))))))))

(define-syntax-rule (dissip3 ii ii_2 i j k)
  (let* ([ii (sub1 ii_2)]
         [idx (+ (* i isize1) (* j jsize1) (* k ksize1))])
    (for ([m (in-range 5)])
      (let ([midx (+ m idx)])
        (flvs! rhs midx (- (* dssp (+         (flvr u midx-2)
                                      (* -4.0 (flvr u midx-))
                                      (*  6.0 (flvr u midx))
                                      (* -4.0 (flvr u midx+)))))))))
  (let* ([ii ii_2]
         [idx (+ (* i isize1) (* j jsize1) (* k ksize1))])
    (for ([m (in-range 5)])
      (let ([midx (+ m idx)])
        (flvs! rhs midx (- (* dssp (+         (flvr u midx-2)
                                      (* -4.0 (flvr u midx-))
                                      (*  5.0 (flvr u midx))))))))))

(define-syntax-case (DISSIP kk jj ii nkk2 njj2 nii2) 
  (for ([kk (in-range 1 (add1 nkk2))])
    (for ([jj (in-range 1 (add1 njj2))])
      (for ([ii (in-range 1 (add1 nii2))])
        (BODY))
      (dissip1 ii i j k)
      (for ([i (in-range 3 (sub1 nii2))])
        (dissip2 ii i j k))
      (dissip3 ii nii2 i j k))))

(DISSIP k j i nz2 ny2 nx2) 
(DISSIP k i j nz2 nx2 ny2) 
(DISSIP i j k nx2 ny2 nz2) 

  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
      (flvs!* fhs idx dt))))
  

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
(define (verify class rnm2 no_time_steps)
  (define xcrdif (make-flvector 5 0.0))
  (define xcedif (make-flvector 5 0.0))
  (define xcr (make-flvector 5 0.0))
  (define xce (make-flvector 5 0.0))

    int m;
    int verified=-1;
    char clss = 'U';
;;;//---------------------------------------------------------------------
;;;//   compute the error norm and the residual norm, and exit if not printing
;;;//---------------------------------------------------------------------
    (error_norm xce)
    (compute_rhs)
    (rhs_norm xcr)

    (for ([m (in-range 5)]) 
      (flvs!/ xcr m dt)
      (flvs! xcrref m 1.0)
      (flvs! xceref m 1.0))

;;;//---------------------------------------------------------------------
;;;//    reference data for 12X12X12 grids after 100 time steps, with DT = 1.50d-02
;;;//---------------------------------------------------------------------
;;;//---------------------------------------------------------------------
;;;//    Reference values of RMS-norms of residual.
;;;//---------------------------------------------------------------------
;;;//---------------------------------------------------------------------
;;;//    Reference values of RMS-norms of solution error.
;;;//---------------------------------------------------------------------
  (define-values (xcrref xceref dtref)
  (case class
  [(#\S) (values
      (flvector
          2.7470315451339479E-2
          1.0360746705285417E-2
          1.6235745065095532E-2
          1.5840557224455615E-2
          3.4849040609362460E-2)
      (flvector
          2.7289258557377227E-5
          1.0364446640837285E-5
          1.6154798287166471E-5
          1.5750704994480102E-5
          3.4177666183390531E-5)
      0.015)]
  [(#\W) (values
      (flvector
      0.1893253733584E-2
      0.1717075447775E-3
      0.2778153350936E-3
      0.2887475409984E-3
      0.3143611161242E-2)
      (flvector
      0.7542088599534E-4
      0.6512852253086E-5
      0.1049092285688E-4
      0.1128838671535E-4
      0.1212845639773E-3)
      0.0015)]
  [(#\A) (values
      (flvector
      2.4799822399300195
      1.1276337964368832
      1.5028977888770491
      1.4217816211695179
      2.1292113035138280)
      (flvector
      1.0900140297820550E-4
      3.7343951769282091E-5
      5.0092785406541633E-5
      4.7671093939528255E-5
      1.3621613399213001E-4)
      0.0015)]
  [(#\B) (values
      (flvector
      0.6903293579998E+02
      0.3095134488084E+02
      0.4103336647017E+02
      0.3864769009604E+02
      0.5643482272596E+02)
      (flvector
      0.9810006190188E-02
      0.1022827905670E-02
      0.1720597911692E-02
      0.1694479428231E-02
      0.1847456263981E-01)
      0.001)]
  [(#\C) (values
      (flvector
      0.5881691581829E+03
      0.2454417603569E+03
      0.3293829191851E+03
      0.3081924971891E+03
      0.4597223799176E+03)
      (flvector
      0.2598120500183
      0.2590888922315E-01
      0.5132886416320E-01
      0.4806073419454E-01
      0.5483377491301) 
      0.00067)]))

    (for ([m (in-range 5)]) 
      (let ([xcrr (vr xcrref m)]
            [xcer (vr xceref m)])
        (flvs! xcrdif m (/ (abs (- (vr xcr m) xcrr) xcrr)))
        (flvs! xcedif m (/ (abs (- (vr xce m) xcer) xcer)))))

  (define  epsilon 1.0E-8)
  (if (not (= class #\U))
    (begin
      (printf " Verification being performed for class ~a\n" class)
      (printf " Accuracy setting for epsilon = ~a\n" epsilon)
      (if ((abs (- dt dtref)) . <= . epsilon)
        #t
        (begin
          (printf "DT does not match the reference value of ~a\m" dtref)
          #f))
      (printf " Comparison of RMS-norms of residual"))
    (begin
      (printf " Unknown CLASS")
      (printf " RMS-norms of residual")
    ))
#|
    verified=BMResults.printComparisonStatus(clss,verified,epsilon,
                                             xcr,xcrref,xcrdif);
    if (clss != 'U') {
      System.out.println(" Comparison of RMS-norms of solution error");
    }else{
      System.out.println(" RMS-norms of solution error");
    }
    verified=BMResults.printComparisonStatus(clss,verified,epsilon,
                                             xce,xceref,xcedif);

    BMResults.printVerificationStatus(clss,verified,BMName); 
    return verified;
  }
|#
)
(__solve k j i nz2 ny2 nx2 i j k ninvr
us rhon
dx2 dx5 dx1
)
(__solve k i j nz2 nx2 ny2 i j k pinvr
vs rhoq
dy3 dy5 dy1
)
(__solve j i k ny2 nx2 nz2 i j k tzetar
ws rhos
dz4 dz5 dz1
)
(define-syntax-rule (DISSIP kk jj ii nkk2 njj2 nii2 i j k inverter
_s rho_
d__ d_5 d_1
dtt_1 dtt_2
c2dtt_
)
  (for ([kk (in-range 1 (add1 nkk2))])
    (for ([jj (in-range 1 (add1 njj2))])
      (for ([ii (in-range (+ nii2 2))])
        (let* ([idx2 (+ i (* j jsize1) (* k ksize1))]
               [ri (* c3c4 (flvr rho_i idx2))])
          (flvs! cv ii (flvr _s idx2))
          (flvs! rho_ ii (dmax1 (+ d__ (* con43 ru1))
                                (+ d_5 (* c1c5  ru1))  
                                d_1))))
      (lhsinit (sub1 nii2))
      (for ([ii (in-range 1 (+ nii2 1))])
        (let ([iij (* ii jsize4)]
              [ii+ (add1 ii)]
              [ii- (sub1 ii)])
        (flvs! lhs (+ 0 iij) 0.0)
        (flvs! lhs (+ 1 iij) (- (- (* dtt_2 (flvr cv ii-)) (* dtt_1 (flvr rho_ ii-)))))
        (flvs! lhs (+ 2 iij) (+ 1.0 (* c2dtt_ (vr rho_ ii))))
        (flvs! lhs (+ 3 iij)    (- (* dtt_2 (flvr cv ii+)) (* dtt_1 (flvr rho_ ii+))))
        (flvs! lhs (+ 4 iij) 0.0)))

      (let ([ii 1]
            [iij (* ii jsize4)]
            [iij+ (* (add1 ii) jsize4)])
        (flvs!+ lhs (+ 2 iij) comz5)
        (flvs!- lhs (+ 3 iij) comz4)
        (flvs!+ lhs (+ 4 iij) comz1)

        (flvs!- lhs (+ 1 iij+) comz4)
        (flvs!+ lhs (+ 2 iij+) comz5)
        (flvs!- lhs (+ 3 iij+) comz4)
        (flvs!+ lhs (+ 4 iij+) comz1))
    
      (for ([ii (in-range 3 (sub1 nii2))])
        (let ([iij (* ii jsize4)])
          (flvs!+ lhs (+ 0 iij) comz1)
          (flvs!- lhs (+ 1 iij) comz4)
          (flvs!+ lhs (+ 2 iij) comz6)
          (flvs!- lhs (+ 3 iij) comz4)
          (flvs!+ lhs (+ 4 iij) comz1)))

      (let ([ii nii2]
            [iij (* ii jsize4)]
            [iij+ (* (add1 ii) jsize4)])
        (flvs!+ lhs (+ 0 iij) comz1)
        (flvs!- lhs (+ 1 iij) comz3)
        (flvs!+ lhs (+ 2 iij) comz6)
        (flvs!- lhs (+ 3 iij) comz4)

        (flvs!+ lhs (+ 0 iij+) comz1)
        (flvs!- lhs (+ 1 iij+) comz4)
        (flvs!+ lhs (+ 2 iij+) comz5))
      
      (for ([ii (in-range 1 (add1 nii2))])
        (let ([iij (* ii jsize4)]
              [d+ (* dtt_2 (vr speed idx2z+))]
              [d- (* dtt_2 (vr speed idx2z-))])
          (flvs! lhsp (+ 0 iij)    (flvr lhs (+ 0 iij)))
          (flvs! lhsp (+ 1 iij) (- (flvr lhs (+ 1 iij)) d-))
          (flvs! lhsp (+ 2 iij)    (flvr lhs (+ 2 iij)))
          (flvs! lhsp (+ 3 iij) (+ (flvr lhs (+ 3 iij)) d+))
          (flvs! lhsp (+ 4 iij)    (flvr lhs (+ 4 iij)))

          (flvs! lhsm (+ 0 iij)    (flvr lhs (+ 0 iij)))
          (flvs! lhsm (+ 1 iij) (+ (flvr lhs (+ 1 iij)) d-))
          (flvs! lhsm (+ 2 iij)    (flvr lhs (+ 2 iij)))
          (flvs! lhsm (+ 3 iij) (- (flvr lhs (+ 3 iij)) d+))
          (flvs! lhsm (+ 4 iij)    (flvr lhs (+ 4 iij)))))

      (for ([ii (in-range nii2)])
        (let* ([iij (* ii jsize4)]
               [ii1j (* (add1 ii) jsize4)]
               [ii2j (* (+ ii 2) jsize4)]
               [fac1 (/ 1.0 (flvr lhs (+ 2 iij)))])
          (flvs!* lhs (+ 3 iij) fac1)
          (flvs!* lhs (+ 4 iij) fac1)
          (for ([m (in-range 3)])
            (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!* rhs idx fac1)))
          (flvs!- lhs (+ 2 ii1j) (* (flvr lhs (+ 1 ii1j)) (flvr lhs (+ 3 iij))))
          (flvs!- lhs (+ 3 ii1j) (* (flvr lhs (+ 1 ii1j)) (flvr lhs (+ 4 iij))))
          (for ([m (in-range 3)])
            (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!- rhs idxz+ (* (flvr lhs (+ 1 ii1j)) (flvr lhs idx)))))
          (flvs!- lhs (+ 1 ii2j) (* (flvr lhs ii2j) (flvr lhs (+ 3 iij))))
          (flvs!- lhs (+ 2 ii2j) (* (flvr lhs ii2j) (flvr lhs (+ 4 iij))))
          (for ([m (in-range 3)])
            (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!- rhs idxz+2 (* (flvr lhs ii2j) (flvr lhs idx)))))))

      (let* ([ii nii2]
             [iij (* ii jsize4)]
             [ii1j (* (add1 nii2) jsize4)]
             [fac1 (/ 1.0 (vr lhs (+ 2 iij)))])
        (flvs!* lhs (+ 3 iij) fac1)  
        (flvs!* lhs (+ 4 iij) fac1)  

        (for ([m (in-range 3)])
          (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
            (flvs!* rhs idx fac1)))

        (flvs!- lhs (+ 2 ii1j) (* (flvr lhs (+ 1 ii1j)) (flvr lhs (+ 3 iij))))
        (flvs!- lhs (+ 3 ii1j) (* (flvr lhs (+ 1 ii1j)) (flvr lhs (+ 4 iij))))

        (for ([m (in-range 3)])
          (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
            (flvs!- rhs idxz+ (* (flvr lhs (+ 1 ii1j)) (flvr lhs idx)))))
        (let ([fac2 (/ 1.0 (flvr lhs (+ 2 ii1j)))])
            (for ([m (in-range 3)])
            (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!* rhs idxz+ fac2)))))

      (for ([ii (in-range (sub1 nii2))])
        (let* ([iij (* ii jsize4)]
               [ii1j (* (add1 ii) jsize4)]
               [ii2j (* (+ ii 2) jsize4)])
          (let* ([fac1 (/ 1.0 (flvr lhsp (+ 2 iij)))]
                 [m 3]
                 [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
            (flvs!* lhsp (+ 3 iij) fac1)
            (flvs!* lhsp (+ 4 iij) fac1)
            (flvs!* rhs midx (+ 4 iij) fac1)
            (flvs!- lhsp (+ 2 ii1j)    (+ (flvr lhsp (+ 1 ii1j)) (flvr lhsp (+ 3 iij))))
            (flvs!- lhsp (+ 3 ii1j)    (+ (flvr lhsp (+ 1 ii1j)) (flvr lhsp (+ 4 iij))))
            (flvs!- rhs midx (+ 4 iij) (* (flvr lhsp (+ 1 ii1j)) (flvr rhs midx)))
            (flvs!- lhsp (+ 1 ii2j)    (+ (flvr lhsp ii2j) (flvr lhsp (+ 3 iij))))
            (flvs!- lhsp (+ 2 ii2j)    (+ (flvr lhsp ii2j) (flvr lhsp (+ 4 iij))))
            (flvs!- rhs midx (+ 4 iij) (* (flvr lhsp ii2j) (flvr rhs midx))))

          (let* ([fac1 (/ 1.0 (flvr lhsm (+ 2 iij)))]
                 [m 4]
                 [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
            (flvs!* lhsm (+ 3 iij) fac1)
            (flvs!* lhsm (+ 4 iij) fac1)
            (flvs!* rhs midx (+ 4 iij) fac1)
            (flvs!- lhsm (+ 2 ii1j)    (+ (flvr lhsm (+ 1 ii1j)) (flvr lhsm (+ 3 iij))))
            (flvs!- lhsm (+ 3 ii1j)    (+ (flvr lhsm (+ 1 ii1j)) (flvr lhsm (+ 4 iij))))
            (flvs!- rhs midx (+ 4 iij) (* (flvr lhsm (+ 1 ii1j)) (flvr rhs midx)))
            (flvs!- lhsm (+ 1 ii2j)    (+ (flvr lhsm ii2j) (flvr lhsm (+ 3 iij))))
            (flvs!- lhsm (+ 2 ii2j)    (+ (flvr lhsm ii2j) (flvr lhsm (+ 4 iij))))
            (flvs!- rhs midx (+ 4 iij) (* (flvr lhsm ii2j) (flvr rhs midx))))))

      (let* ([ii nii2]
             [iij (* ii jsize4)]
             [ii1j (* (add1 ii) jsize4)])
        (let* ([fac1 (/ 1.0 (flvr lhsp (+ 2 iij)))]
               [m 3]
               [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
          (flvs!* lhsp (+ 3 iij) fac1)
          (flvs!* lhsp (+ 4 iij) fac1)
          (flvs!* rhs midx (+ 4 iij) fac1)
          (flvs!- lhsp (+ 2 ii1j)    (+ (flvr lhsp (+ 1 ii1j)) (flvr lhsp (+ 3 iij))))
          (flvs!- lhsp (+ 3 ii1j)    (+ (flvr lhsp (+ 1 ii1j)) (flvr lhsp (+ 4 iij))))
          (flvs!- rhs midx (+ 4 iij) (* (flvr lhsp (+ 1 ii1j)) (flvr rhs midx))))
        (let* ([fac1 (/ 1.0 (flvr lhsm (+ 2 iij)))]
               [m 4]
               [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
          (flvs!* lhsm (+ 3 iij) fac1)
          (flvs!* lhsm (+ 4 iij) fac1)
          (flvs!* rhs midx (+ 4 iij) fac1)
          (flvs!- lhsm (+ 2 ii1j)    (+ (flvr lhsm (+ 1 ii1j)) (flvr lhsm (+ 3 iij))))
          (flvs!- lhsm (+ 3 ii1j)    (+ (flvr lhsm (+ 1 ii1j)) (flvr lhsm (+ 4 iij))))
          (flvs!- rhs midx (+ 4 iij) (* (flvr lhsm (+ 1 ii1j)) (flvr rhs midx))))

        (flvs!/ rhs (+ 3 idxz+) (flvr lhsp (+ 2 ii1j)))
        (flvs!/ rhs (+ 4 idxz+) (flvr lhsm (+ 2 ii1j))))

      (let* ([ii nii2]
             [iij (* ii jsize4)]
             [ii1j (* (add1 ii) jsize4)]
             [idx (* i isize1) (* j jsize1) (* k ksize1)])
        (for ([m (in-range 3)])
          (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
            (flvs!- rhs idx (* (flvr lhs (+ 3 ii1j)) (flvr lhs idx)))))

        (flvs!- rhs (+ 3 idx)    (+ (flvr lhsp (+ 3 iij)) (flvr rhs (+ 3 idx))))
        (flvs!- rhs (+ 4 idx)    (+ (flvr lhsm (+ 3 iij)) (flvr rhs (+ 4 idx)))))

      (for ([ii (in-range (- nii2 1) -1 -1)])
        (let* ([iij (* ii jsize4)]
               [idx1 (+ (* (+ i 1) isize1) (* j jsize1) (* k ksize1))]
               [idx2 (+ (* (+ i 2) isize1) (* j jsize1) (* k ksize1))])
          (for ([m (in-range 3)])
            (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
              (flvs!- rhs idx (- (* (flvr lhs (+ 3 iij)) (flvr rhs (+ m idx1)))
                                 (* (flvr lhs (+ 4 iij)) (flvr rhs (+ m idx2)))))))

          (flvs!- rhs (+ 3 idx)    (* (flvr lhsp (+ 3 iij)) (flvr rhs (+ 3 idx1)))
                                   (* (flvr lhsp (+ 4 iij)) (flvr rhs (+ 3 idx2))))
          (flvs!- rhs (+ 4 idx)    (* (flvr lhsm (+ 3 iij)) (flvr rhs (+ 4 idx1)))
                                   (* (flvr lhsm (+ 4 iij)) (flvr rhs (+ 4 idx2))))))
      (inverter))))

; (- nz2 1)  = grid_points[2]-3
; nz2        = grid_points[2]-2
; (+ nz2 1)  = grid_points[2]-1

(define x_solve)
  
)

(define (checkSum arr)
  (for*/fold ([csum 0.0]) ([k (in-range (add1 nz2))]
         [j (in-range (add1 ny2))]
         [i (in-range (add1 nx2))]
         [m (in-range 5)])
    (let* ([offset (+ m (* i isize1) (* j jsize1) (* k ksize1))]
           [arro (flvr arr offset)])
      (+ csum (/ (sqr arro) (* (+ 2 nz2) (+ 2 ny2) (+ 2 nx2) 5))))))

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

