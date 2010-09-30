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
  (let ([args (parse-cmd-line-args argv "LU Decomposition")]) 
    (run-benchmark args)))

(define make-fxvector make-vector)

(define (run-benchmark args) 
  (define maxlevel 11)
  (let ([bmname "LU"]
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

          [dxi (/ 1.0  (- nx0 1))]
          [deta (/ 1.0  (- ny0 1))]
          [dzeta (/ 1.0  (- nz0 1))]
          [tx1 (/ 1.0  (* dxi dxi))]
          [tx2 (/ 1.0  (* 2.0 dxi))]
          [tx3 (/ 1.0  dxi)]
          [ty1 (/ 1.0  (* deta deta))]
          [ty2 (/ 1.0  (* 2.0 deta))]
          [ty3 (/ 1.0  deta)]
          [tz1 (/ 1.0  (* dzeta dzeta))]
          [tz2 (/ 1.0  (* 2.0 dzeta ))]
          [tz3 (/ 1.0  dzeta)]
;;;//---------------------------------------------------------------------
;;;//   diffusion coefficients
;;;//---------------------------------------------------------------------
          [dx1  0.75]
          [dx2  dx1]
          [dx3  dx1]
          [dx4  dx1]
          [dx5  dx1]
          [dy1  0.75]
          [dy2  dy1]
          [dy3  dy1]
          [dy4  dy1]
          [dy5  dy1]
          [dz1  1.00]
          [dz2  dz1]
          [dz3  dz1]
          [dz4  dz1]
          [dz5  dz1]
;;;//---------------------------------------------------------------------
;;;//   fourth difference dissipation
;;;//---------------------------------------------------------------------
          [dssp (/ (max dx1 dy1 dz1) 4.0)]
 

      ;(get-input-pars)
      (domain)

      ;(set-coefficients)
      (set-boundary-variables)
      (set-initial-values
      (compute-forcing-term)
      (timer-start 1)

      (SSOR-serial)
    
      (compute-error)

      (compute-surface-integral)

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

(define (get-mflops total-time niter nx ny nz)
  (if (not (= total-time 0.0))
    (let* ([n3 (* nx ny nz)]
           [t  (/ (+ nx ny nz) 3.0)])
      (/ (* (+ (* 1984.77 n3)
               (* -10923.3 (sqr t))
               (* 27770.9 t)
               -144010.0)
            niter)
          (* total-time 1000000.0)))
      0.0))


(define-syntax-rule (IDX m i j k) (+ m (* i isize1) (* j jsize1) (* k ksize1)))

(bltw m i j k - M1)
(bltw m i j k + M0)

(define-syntax-rule (bltw m i j k M1 M0 XM ij1 ijk1 BODY ...)
  (define-syntax-rule (M1 a ...) (- (fr v mijk) a ...))
  (define-syntax-rule (M0 a ...) (begin a ...))

  (for ([j (in-range (sub1 jst) jend)]
        [i (in-range (sub1 ist) iend)])
    (let ([ij1  (+ (* i isize1) (* j jsize1))]
          [ijk  (+ (* i isize1) (* j jsize1) (* k ksize1))]
    (for ([m (in-range 5)])
      (let ([mij1  (+ m ij1)]
            [mjk4  (+ m (* i jsize4) (* j ksize4))]
            [mijk  (+ m ijk1)]
            [i-jk (+ (* (IM i) isize1) (* j jsize1) (* k ksize1))]
            [ij-k (+ (* i isize1) (* (JM j) jsize1) (* k ksize1))]
            [ijk- (+ (* i isize1) (* j jsize1) (* (KM k) ksize1))])
        (f! tv mij1 (MODIFIER
            (* omega
              (for/fold ([T 0.0]) ([n (in-range 5)])
                (+ T
                   (* (fr ldz (+ (* n isize4) mjk4)) (fr v (+ n ijk-)))
                   (* (fr ldy (+ (* n isize4) mjk4)) (fr v (+ n ij-j)))
                   (* (fr ldx (+ (* n isize4) mjk4)) (fr v (+ n i-jk))))))))))
      
      (let ([jk4 (+ (* i jsize4) (* j ksize4))])
        (for ([m (in-range 5)])
          (define-syntax-rule (doit I1 I2) (f! tmat (+ m I1) (fr d (+ m (* I2 isize4) jk4)))) 
          (doit 0  0) 
          (doit 5  1)
          (doit 10 2)
          (doit 15 3)
          (doit 20 4)))

      (let* ([C 0] 
             [tmp1 (/ 1.0 (fr tmat (+ C (* C 5))))])
        (define-syntax-rule (doit A)
          (let ([tmp (* tmp1 (fr tmat (+ A (* C 5))))])
            (define-syntax-rule (DD V B) (f!- V (+ A B) (* tmp (fr V (+ C B)))))
            (define-syntax-rule (D B) (DD tmat B))
            (D 5) (D 10) (D 15) (D 20)
            (DD tv ijx)))
        (doit 1) (doit 2) (doit 3) (doit 4))

      (let* ([C 1] 
             [tmp1 (/ 1.0 (fr tmat (+ C (* C 5))))])
        (define-syntax-rule (doit A)
          (let ([tmp (* tmp1 (fr tmat (+ A (* C 5))))])
              (define-syntax-rule (DD V B) (f!- V (+ A B) (* tmp (fr V (+ C B)))))
            (define-syntax-rule (D B) (DD tmat B))
            (D 10) (D 15) (D 20)
            (DD tv ijx)))
        (doit 2) (doit 3) (doit 4))
      (let* ([C 2] 
             [tmp1 (/ 1.0 (fr tmat (+ C (* C 5))))])
        (define-syntax-rule (doit A)
          (let ([tmp (* tmp1 (fr tmat (+ A (* C 5))))])
            (define-syntax-rule (DD V B) (f!- V (+ A B) (* tmp (fr V (+ C B)))))
            (define-syntax-rule (D B) (DD tmat B))
            (D 15) (D 20)
            (DD tv ijx)))
        (doit 3) (doit 4))

      (let* ([C 3] 
             [tmp1 (/ 1.0 (fr tmat (+ C (* C 5))))])
        (define-syntax-rule (doit A)
          (let ([tmp (* tmp1 (fr tmat (+ A (* C 5))))])
            (define-syntax-rule (DD V B) (f!- V (+ A B) (* tmp (fr V (+ C B)))))
            (define-syntax-rule (D B) (DD tmat B))
            (D 20)
            (DD tv ijx)))
        (doit 4))

      BODY ...))


(define (blts ldmx ldmy ldmz nx ny nz k omega v tv ldz ldy ldx d ist iend jst jend nx0 ny0)
  (bltw m i j k - MO M1 M1 ij1 ijk1
    (f! v (+ 4 ijk1) (/ (fr tv (+ 4 ij1) (fr tmat (+ 4 20)))))
    (f!- tv (+ 3 ij1) (* (fr tmat (+ 3 20)) (fr v (+ 4 ijk1))))

    (f! v (+ 3 ijk1) (/ (fr tv (+ 3 ij1) (fr tmat (+ 3 15)))))
    (f!- tv (+ 2 ij1) (* (fr tmat (+ 2 15)) (fr v (+ 3 ijk1)))
                      (* (fr tmat (+ 2 20)) (fr v (+ 4 ijk1))))

    (f! v (+ 2 ijk1) (/ (fr tv (+ 2 ij1) (fr tmat (+ 2 10)))))
    (f!- tv (+ 1 ij1) (* (fr tmat (+ 1 10)) (fr v (+ 2 ijk1)))
                      (* (fr tmat (+ 1 15)) (fr v (+ 3 ijk1)))
                      (* (fr tmat (+ 1 20)) (fr v (+ 4 ijk1))))

    (f! v (+ 1 ijk1) (/ (fr tv (+ 1 ij1) (fr tmat (+ 1  5)))))
    (f!- tv (+ 0 ij1) (* (fr tmat (+ 0  5)) (fr v (+ 1 ijk1)))
                      (* (fr tmat (+ 0 15)) (fr v (+ 2 ijk1)))
                      (* (fr tmat (+ 0 15)) (fr v (+ 3 ijk1)))
                      (* (fr tmat (+ 0 20)) (fr v (+ 4 ijk1))))

    (f! v (+ 0 ijkx) (/ (fr tv (+ 0 ijx) (fr tmat (+ 0  0)))))
))

(define (buts ldmx ldmy ldmz nx ny nz k omega v tv ldz ldy ldx d ist iend jst jend nx0 ny0)
  (bltw m i j k + M0 ij1 ijk1
    (let ([tv4 (/ (fr tv (+ 4 ij1)) (fr tmat 24))])
      (f! tv (+ 4 ij1) tv4)
      (let ([tv3 (/ (- (fr tv (+ 3 ij1)) 
                           (* (fr tmat 23) t4)) (fr tmat 18))])
        (f! tv (+ 3 ij1) tv3)
        (let ([tv2 (/ (- (fr tv (+ 2 ij1)) 
                             (* (fr tmat 17) tv3)
                             (* (fr tmat 22) tv4)) (fr tmat 12))])
          (f! tv (+ 2 ij1) tv2)
          (let ([tv1 (/ (- (fr tv (+ 1 ij1)) 
                             (* (fr tmat 11) tv2)
                             (* (fr tmat 16) tv3)
                             (* (fr tmat 21) tv4)) (fr tmat 6))])
            (f! tv (+ 1 ij1) tv1)
            (let ([tv0 (/ (- (fr tv (+ 0 ij1)) 
                               (* (fr tmat  5) tv1)
                               (* (fr tmat 10) tv2)
                               (* (fr tmat 15) tv3)
                               (* (fr tmat 20) tv4)) (fr tmat 0))])
              (f! tv (+ 0 ij1) tv0)
              (f!- v (+ 0 ijk1) tv0)
              (f!- v (+ 1 ijk1) tv1)
              (f!- v (+ 2 ijk1) tv2)
              (f!- v (+ 3 ijk1) tv3)
              (f!- v (+ 4 ijk1) tv4))))))))

(define (domain nx ny nz)
  (when (or ( nx . < . 4 ) ( ny . < . 4 ) ( nz . < . 4 ))
    (printf "     SUBDOMAIN SIZE IS TOO SMALL - \n")
    (printf "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n")
    (printf "     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n")
    (printf "     TO 4 THEY ARE CURRENTLY ~a ~a ~a\n" nx ny nz)
    (exit 0))
     

  (when (or ( nx . > . isiz1 ) ( ny . > . isiz2 ) ( nz . > . isiz3 ))
    (printf "     SUBDOMAIN SIZE IS TOO LARGE - \n")
    (printf "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n")
    (printf "     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL\n")
    (printf "     TO ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY. THEY ARE\n")
    (printf "     CURRENTLY ~a ~a ~a\n" nx ny nz)
    (exit 0)))

(define-syntax-rule (fourth-order-dissipation 
  (let* ([zmost ZMOST])
    (let ([i1 (+ (* 1 ZS) zmost)]
          [i2 (+ (* 2 ZS) zmost)]
          [i3 (+ (* 3 ZS) zmost)]
          [i4 (+ (* 4 ZS) zmost)])
      (for ([m (in-range 5)])
        (let* ([m1 (+ m i1)]
               [m2 (+ m i2)]
               [m3 (+ m i3)]
               [m4 (+ m i4)]
               [didx1  m1])
               [didx2  m2])
          (f!- V didx1 (* dssp (+ (*  5.0 (fr V2 m1))
                                  (* -4.0 (fr V2 m2))
                                          (fr V2 m3))))
          (f!- V didx2 (* dssp (+ (* -4.0 (fl V2 m1))
                                  (*  6.0 (fl V2 m2))
                                  (* -4.0 (fl V2 m3))
                                          (fl V2 m4))))))
    (for ([ii (in-range 3 (- nx 3))])
        (let ([idx IDX])
           [i-2 (+ (* (- ii 2) ZS) zmost)]
           [i-1 (+ (* (- ii 1) ZS) zmost)]
           [ix  (+ (* ii       ZS) zmost)]
           [i+1 (+ (* (+ ii 1) ZS) zmost)]
           [i+2 (+ (* (+ ii 1) ZS) zmost)]
          (for ([m (in-range 5)])
            (let* ([m-2 (+ m i-2)]
                   [m-1 (+ m i-1)]
                   [mx  (+ mx ix)]
                   [m+1 (+ m i+1)]
                   [m+2 (+ m i+s)]
                   [didx  mx])
              (f!- V didx (* dssp (+         (flvr V2 m-2)
                                     (* -4.0 (flvr V2 m-1))
                                     (*  6.0 (flvr V2 mx))
                                     (* -4.0 (flvr V2 m+1))
                                             (flvr V2 m+2))))))))
    (let ([i5 (+ (* (- NZ 5) ZS) zmost)]
          [i4 (+ (* (- NZ 4) ZS) zmost)]
          [i3 (+ (* (- NZ 3) ZS) zmost)]
          [i2 (+ (* (- NZ 2) ZS) zmost)])
      (for ([m (in-range 5)])
        (let* ([m5 (+ m i5)]
               [m4 (+ m i4)]
               [m3 (+ m i3)]
               [m2 (+ m i2)]
               [didx1  m3])
               [didx2  m2])
          (f!- V didx1 (* dssp (+         (fl V2 m5)
                                  (* -4.0 (fl V2 m4))
                                  (*  6.0 (fl V2 m3))
                                  (* -4.0 (fl V2 m2)))))
          (f!- V didx2 (* dssp (+         (fr V2 m4)
                                  (* -4.0 (fr V2 m3))
                                  (*  5.0 (fr V2 m2)))))))))
(define (compute-forcing-term)
  (zero-matrix frct)

  (for ([k (in-range nz)])
    (let ([zeta (* k dnzm1)])
      (for ([j (in-range ny)])
        (let ([eta (* j dnym1)])
          (for ([i (in-range nx)])
            (let ([xi (* i dnxm1)])
              (for ([m (in-range 5)]) 
                (f! rsd (+ m (* i isize1) (* j jsize1) (* k ksize1))
                   (+ (fr ce m)
                      (* xi (+ (fr ce (+ m 5))
                               (* xi (+ (fr ce (+ m (* 4 5)))
                                        (* xi (+ (fr ce (+ m (* 7 5)))
                                                 (* xi (fr ce (+ m (* 10 5))))))))))
                      (* eta (+ (fr ce (+ m (* 2 5)))
                               (* eta (+ (fr ce (+ m (* 5 5)))
                                        (* eta (+ (fr ce (+ m (* 8 5)))
                                                 (* eta (fr ce (+ m (* 11 5))))))))))
                      (* zeta (+ (fr ce (+ m (* 3 5)))
                               (* zeta (+ (fr ce (+ m (* 6 5)))
                                        (* zeta (+ (fr ce (+ m (* 9 5)))
                                                 (* zeta (fr ce (+ m (* 12 5)))))))))))))))))))

(define-syntax-rule (NOOP a ...) (void))
(define-syntax-rule (RHSU a ...) (begin a ...))
(define-syntax-rule (SEL1 a b c) a)
(define-syntax-rule (SEL2 a b c) b)
(define-syntax-rule (SEL3 a b c) c)
(define-syntax-rule (MULT43 a ...) (* D43 a ...))
(define-syntax-rule (NOOPF a ...) (begin a ...))

(flux-flow i j k k j i 1 jst 0 (sub1 ist) (sub1 nz) jend nx iend rsd   rsd ijk1 frct MULT43 NOOPF NOOPF NOOP SEL1 tx3
(flux-flow i j k k i j 1 ist 0 (sub1 jst) (sub1 nz) iend jend ny rsd   rsd ijk1 frct NOOPF MULT43 NOOPF NOOP SEL2 ty3
(flux-flow i j k j i k jst ist 0 1 jend iend nz (sub1 nz)        rsd   rsd ijk1 frct NOOPF NOOPF MULT43 NOOP SEL3 tz3
(flux-flow i j k k j i 1 jst 0 (sub1 ist) (sub1 nz) jend nx iend   u rho_i ijk3  rsd MULT43 NOOPF NOOPF RHSU SEL1 tx3
(flux-flow i j k k i j 1 ist 0 (sub2 jst) (sub1 nz) iend ny jend   u rho_i ijk3  rsd NOOPF MULT43 NOOPF RHSU SEL2 ty3
(flux-flow i j k j i k jst ist 0 1 jend iend nz (sub1 nz)          u rho_i ijk3  rsd NOOPF NOOPF MULT43 RHSU SEL3 tz3

(define-syntax-rule (flux-flow i j k kk jj ii skk sjj sii sii2 nkk njj nii nii2 RSD/U RHO RHOI FRCT FM1 FM2 FM3 RHSU SELECTOR)
  t_3 d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1)
  (for ([kk (in-range skk nkk)])
    (for ([jj (in-range sjj njj)])
      (for ([ii (in-range sii nii)])
        (let* ([ijk1 (+ (* i isize1) (* j jsize1) (* k ksize1))]
        (let* ([ijk3 (+ i (* j jsize3) (* k ksize3))]
               [ii2 (* ii isize2)]
               [rhotmp (/ 1.0 (fr RSD/u ijk1))]
               [u21 (* (fr RSD/U (+ A ijk1)) rhotmp)]
               [q (* 0.50 (* (+ (sqr (fr RSD/U (+ 1 ijkx)))
                                (sqr (fr RSD/U (+ 2 ijkx)))
                                (sqr (fr RSD/U (+ 3 ijkx))))
                              rhotmp))])
          (RHSU
            (f! rho_i ijk3  rhotmp)
            (f! qs ijk3 q))

          (f! flux (+ 0 iix) (* (fr RSD/U (+ 1 ijk1))))
          (f! flux (+ 1 iix) (+ (* (fr RSD/U (+ 1 ijk1)) u21)
                                (* c2 (- (fr RSD/U (+ 4 ijkx)) q))))
          (f! flux (+ 2 iix) (* (fr RSD/U (+ 2 ijkx)) u21))
          (f! flux (+ 3 iix) (* (fr RSD/U (+ 3 ijkx)) u21))
          (f! flux (+ 4 iix) (* u21 (- (* c1 (fr rsd (+ 4 ijkx))) 
                                       (* c2 q))))))
      (for ([ii (in-range sii2 nii2)])
        (let ([ii+2 (* (+ ii 1) isize2)]
              [ii-2 (* (- ii 1) isize2)])
          (for ([m (in-range 5)])
            (f!- FRCT mijk1 (* t_2 (- (fr flux (+ m ii+2))
                                      (fr flux (+ m ii-2)))))))

      (for ([ii (in-range sii2 nii)])
        (let* ([ijk1 (+ (* i isize1) (* j jsize1) (* k ksize1))]
               [i-jk1 (+ (* (- i 1)  isize1) (* j jsize1) (* k ksize1))]
               [iix (* ii isize2)]
               [tmp (/ 1.0 (fr RHO (+ 0 RHOI)))]
               [u2li (* tmp (fr RSD/U (+ 1 ijk1)))]
               [u3li (* tmp (fr RSD/U (+ 2 ijk1)))]
               [u4li (* tmp (fr RSD/U (+ 3 ijk1)))]
               [u5li (* tmp (fr RSD/U (+ 4 ijk1)))]
               [tmp2 (/ 1.0 (fr RHO (+ 0 RHOI-)))]
               [u2lim1 (* tmp2 (fr RSD/U (+ 1 i-jkx)))]
               [u3lim1 (* tmp2 (fr RSD/U (+ 2 i-jkx)))]
               [u4lim1 (* tmp2 (fr RSD/U (+ 3 i-jkx)))]
               [u5lim1 (* tmp2 (fr RSD/U (+ 4 i-jkx)))])
          (f! flux (+ 1 iix) (FM1 (* t_3 (- u2li u2lim1))))
          (f! flux (+ 2 iix) (FM2 (* t_3 (- u3li u3lim1))))
          (f! flux (+ 3 iix) (FM3 (* t_3 (- u4li u4lim1))))
          (f! flux (+ 4 iix) (+ (* 0.50 (- 1.0 c1c5)
                                   t_3 (- (+ (sqr u2li) (sqr u3li) (sqr u4li))
                                          (+ (sqr u2lim1) (sqr u3lim1) (sqr u4lim1))))
                                (/ 1.0 6.0) 
                                (* t_3 (- (sqr (SELECTOR u2li u3li u4li)) 
                                          (sqr (SELECTOR u2lim1 u3lim1 u4lim1)))
                                (* c1c5 t_3 (- u5li u5lim1))))))

      (for ([ii (in-range sii2 nii2)])
        (let* ([ijk1 (+ (* i isize1) (* j jsize1) (* k ksize1))]
               [iix (* ii isize2)]
               [iix+ (* (add1 ii) isize2)])
            (define-syntax-rule (citA C C1 C2 C3) (* C (+ (- C1 (* 2.0 C2)) C3)))
            (define-syntax-rule (cit3S C V I1 I2 I3) (citA C (flvr V I1) (flvr V I2) (flvr V I3)))
            (define-syntax-rule (mid d__t_1 AA)
              (let ([idxA (+ ijk1 AA)]
                    [idxz+A (+ idxz+ AA)]
                    [idxz-A (+ idxz- AA)])
                (flvs!+ FRCT idxA
                  (cit3S d__t_1 rsd idxz+A idxA idxz-A)
                  (* t_3 c3 c4 (- (fr flux (+ AA ii+x))
                                  (fr flux (+ AA iix)))))))

            (flvs!+ FRCT ijkx (cit3S d_1t_1 RSD/U idxz+ idx idxz-))
            (mid d_2t_1 1)
            (mid d_3t_1 2)
            (mid d_4t_1 3)
            (mid d_5t_1 4)))
      
      (fourth-order-dissipation FRCT RSD/U ZMOST ZS NZ))))

(define (erhs)
  (compute-forcing-term)
#|
  (for ([k (in-range nz)])
    (for ([j (in-range ny)])
      (for ([i (in-range nx)])
        (for ([m (in-range 5)])
          (f! rsd midx (vr frct midx)))
      (let ([tmp (/ 1.0 (fr u idx))])
        (f! rho_i idx tmp)
        (f! qs ijk3 (* 0.50 tmp (+ (sqr (fr u (+ 1 idx)))
                                   (sqr (fr u (+ 2 idx)))
                                   (sqr (fr u (+ 3 idx))))))))))
|#

  (flux-flow i j k k j i 1 jst 0 (sub1 ist) (sub1 nz) jend nx iend rsd   rsd ijk1 frct MULT43 NOOPF NOOPF NOOP SEL1 tx3
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)
  (flux-flow i j k k i j 1 ist 0 (sub1 jst) (sub1 nz) iend jend ny rsd   rsd ijk1 frct NOOPF MULT43 NOOPF NOOP SEL2 ty3
    dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)
  (flux-flow i j k j i k jst ist 0 1 jend iend nz (sub1 nz)        rsd   rsd ijk1 frct NOOPF NOOPF MULT43 NOOP SEL3 tz3
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1))

(define (rhs)
  (for ([k (in-range nz)])
    (for ([j (in-range ny)])
      (for ([i (in-range nx)])
        (for ([m (in-range 5)])
          (f! rsd midx (vr frct midx))))))
|#
      (let ([tmp (/ 1.0 (fr u idx))])
        (f! rho_i idx tmp)
        (f! qs ijk3 (* 0.50 tmp (+ (sqr (fr u (+ 1 idx)))
                                   (sqr (fr u (+ 2 idx)))
                                   (sqr (fr u (+ 3 idx))))))))))
|#

  (flux-flow i j k k j i 1 jst 0 (sub1 ist) (sub1 nz) jend nx iend   u rho_i ijk3  rsd MULT43 NOOPF NOOPF RHSU SEL1 tx3
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)
  (flux-flow i j k k i j 1 ist 0 (sub2 jst) (sub1 nz) iend ny jend   u rho_i ijk3  rsd NOOPF MULT43 NOOPF RHSU SEL2 ty3
    dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)
  (flux-flow i j k j i k jst ist 0 1 jend iend nz (sub1 nz)          u rho_i ijk3  rsd NOOPF NOOPF MULT43 RHSU SEL3 tz3
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1))

(define-syntax-rule (flux-flow i j k kk jj ii skk sjj sii sii2 nkk njj nii nii2 RSD/U RHO RHOI FRCT FM1 FM2 FM3 RHSU SELECTOR)
  t_3 d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1)
  (flux-flow k j i)
  (flux-flow k i j)
  (flux-flow j i k))

(define (error nx0 ny0 nz0)
  (define n03 (* (- nx0 2) (- ny0 2) (- nz0 2)))
  (for ([m (in-range 5)]) (f! errnm m 0.0))

  (for* ([k (in-range 1 (- nz 1))]
         [j (in-range (- jst 1) jend)]
         [i (in-range (- ist 1) iend)])
    (exact (add1 i) (add1 j) (add1 k) u000ijk)
    (for ([m (in-range 5)])
      (let* ([ijkx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
        (f!+ errnm m (sqr (- (fr u000ijk m) (fr u mijkx)))))))

  (for ([m (in-range 5)]) 
    (f! errnm m (/ (sqrt (fr errnm m)) n03))))
  
(define-syntax-rule (DIAG V BB0 BB1 BB2 CCC1 CCC2 CCC3 t_1 t_2 d__ 
  v2 v3 v4 V342 V343 V344 V01 V02 V03 ijkz3-1 NEG NNEG)
  (let* (; CONSTATNS
         [dtt_1 (* dt t_1)]
         [dtt_2 (* dt t_2)]
         [dtt_1d__ (* dtt_1 d__)]
         [dtt_1d_1 (* dtt_1 d_1)]

         ; LOOP DEPENDENT
         [tmp1 (fr rho_i ijkz3-1)]
         [tmp2 (sqr tmp1)]
         [tmp3 (* tmp1 tmp2)]
         [uBB0 (fr u (+ BB0 ijkz3-1))]
         [uBB1 (fr u (+ BB1 ijkz3-1))]
         [uBB2 (fr u (+ BB2 ijkz3-1))]
         [U1 (fr u (+ 1 ijkz3-1))]
         [U2 (fr u (+ 2 ijkz3-1))]
         [U3 (fr u (+ 3 ijkz3-1))]
         [U4 (fr u (+ 4 ijkz3-1))]
         [QS (fr qs jkz3-1x3)]
         [dtt_1tmp1 (* dtt_1 tmp1)]
         [dtt_2tmp1 (* dtt_2 tmp1)]
         [dtt_1tmp2 (* dtt_1 tmp2)]
         [dtt_2tmp2 (* dtt_2 tmp2)]
         [dtt_2tmp1c2 (* dtt_2tmp1 c2)]
         [dtt_1tmp2c34 (* dtt_1tmp2 c34)]
         [dtt_1tmp1c34 (* dtt_1tmp1 c34)]
         [dtt_2tmp1uBB0 (* dtt_2tmp1 uBB0)]
         [dtt_2tmp2uBB0 (* dtt_2tmp2 uBB0)])

    (let ([v2 0.0]
          [v3 0.0]
          [v4 (NEG dtt_2)])
      (DIAG0 V 0 (- dtt_1d_1) V2 V3 V4 0.0))

    (define-syntax-rule (RX1/2 BB UU v2 v3 v4 V2 V3 V4)
      (let* ([v1 (NEG (* UU (+ dtt_2tmp2uBB0 
                               dtt_1tmp2c34)))]
             [v2 (NEG (+ dtt_2tmp1uBB0
                       dtt_1tmp1c34
                       dtt_1d__))]
             [v3 0.0]
             [v4 (NEG (* dtt_2tmp1 UU))])
        (DIAG0 V BB v1 V2 V3 V4 0.0)))
    (define-syntax-rule (RX1 v2 v3 v4 V2 V3 V4)
      (RX1/2 BB1 uBB1 v2 v3 v4 V2 V3 V4))
    (define-syntax-rule (RX2 v2 v3 v4 V2 V3 V4)
      (RX1/2 BB2 uBB2 v3 v2 v4 V2 V3 V4))

    (define-syntax-rule (RX3 v2 v3 v4 V2 V3 V4)
      (let* ([v1 (NEG (* dtt_2 (- (* c2 QS tmp1) 
                             (sqr (* uBB0 tmp1)))))]
             [v2 (NNEG (* dtt_2tmp1c2 uBB1))]
             [v3 (NNEG (* dtt_2tmp1c2 uBB2))]
             [v4 (+ (NEG (* dtt_2tmp1uBB0 (- 2.0 c2)))
                    (-   (* dtt_1tmp1c34 r43 ))
                    (- dtt_1d__))]
             [v5 (NEG (* dtt_2 c2))])
        (DIAG0 V BB0 v1 V2 V3 V4 v5)))
              
    (define-syntax-rule (RX4 v2 v3 v4 V2 V3 V4)
      (let* ([v1 (+ (NEG (* dtt_2tmp2uBB0 (- (* c2 QS) (* c1 U4))))
                    (-   (* dtt_1 (- (+ (* CCC1   tmp3 (sqr U1))
                                   (* CCC2   tmp3 (sqr U2))
                                   (* CCC3   tmp3 (sqr U3))
                                   (* c1345 tmp2      U4))))))]
             [dtt_2tmp2c2 (NEG (* dtt_2tmp2 c2))]
             [dtt_1tmp2C  (-  (* dtt_1tmp2 (- c34 c1345)))]
             [BIGC (+ dtt_2tmp2c2 dtt_1tmp2C)]
             [v2 (* BIGC uBB1)]
             [v3 (+ BIGC uBB2)]
             [v4 (+ (NEG (* dtt_2 (- (* c1 U4 tmp1)
                                (* c2 (+ (* (sqr uBB0) tmp2)
                                         (* QS  tmp1))))))
                    (-   (* dtt_1tmp2uBB0 (- (* r43 c34) c1345))))]

             [v5 (+ (NEG (* dtt_2tmp1uBB0 c2))
                    (-   (* dtt_1tmp1 r1345))
                    (-   (* dtt_1 d_5)))]
        (DIAG0 V 4 v1 V2 V3 V4 v5))))
             
    (RX1 V v2 v3 v4  v342 v343 v344)
    (RX2 V v2 v3 v4  v342 v343 v344)
    (RX3 V v2 v3 v4  v342 v343 v344)
    (RX4 V v2 v3 v4  v342 v343 v344)
))

(define-syntax-rule (NEG a ...) (- a ...))
(define-syntax-rule (NNEG a ...) (begin a ...))

(define (jacld k)
  (for [j (in-range (- jst 1) jend)]
    (let ([jk4 (+ (* j jsize4) (* k ksize4))])
      (for ([i (in-range (- ist 1) iend)])

    (let ([C1 (- c34 c1345 )]
          [C2 (- (* r43 c34) c1345)])
      (define-syntax-rule (C1C2__ C1 C2 C3) (+ (* tx1 C1) (* ty1 C2) (* tz1 C3)))

    ;; CONSTANTS
    (let ([r43 (/ 4.0 3.0)]
          [c1345 (* c1 c3 c4 c5)]
          [c34 (* c3 c4)]
          [t_1s (+ tx1 ty1 tz1)]
          [tds1 (+ tx1dx1 ty1dy1 tz1dz1)]
          [tds2 (+ tx1dx2 ty1dy2 tz1dz2)]
          [tds3 (+ tx1dx3 ty1dy3 tz1dz3)]
          [tds4 (+ tx1dx4 ty1dy4 tz1dz4)]
          [tds5 (+ tx1dx5 ty1dy5 tz1dz5)]
          [C1C2_1 (C1C2__ C2 C1 C1)]
          [C1C2_2 (C1C2__ C1 C2 C1)]
          [C1C2_3 (C1C2__ C1 C1 C2)]

          ;; LOOP DEPENDENT VARIABLES
          [ijk1 (+ (* i isize1) (* j jsize1) (* k ksize1))]
          [ijk3 (+    i         (* j jsize3) (* k ksize3))]
          [i0jk4  jk4]
          [i1jk4  (+ (* 1 isize4) (* j jsize4) (* k ksize4))]
          [i2jk4  (+ (* 2 isize4) (* j jsize4) (* k ksize4))]
          [i3jk4  (+ (* 3 isize4) (* j jsize4) (* k ksize4))]
          [i4jk4  (+ (* 4 isize4) (* j jsize4) (* k ksize4))])

      (define-syntax-rule (DIAG0 V A V1 V2 V3 V4 V5)
        (f! V (+ A i0jk4) V1)
        (f! V (+ A i1jk4) V2)
        (f! V (+ A i2jk4) V3)
        (f! V (+ A i3jk4) V4)
        (f! V (+ A i4jk4) V5))
     
;;;//---------------------------------------------------------------------
;;;//   form the block daigonal
;;;//---------------------------------------------------------------------
    (let* ([tmp1 (fr rho_i ijkx3)]
           [tmp2 (sqr tmp1)]
           [tmp3 (* tmp1 tmp2)])
 
      (define-syntax-rule (RX123 A v2 v3 v4 V2 V3 V4 t_143 t_1_d__)
        (let ([v1 (- (* dt 2.0 t_1r43 c34 tmp2 (fr u (+ A ijk1))))]
              [v2 (+ 1.0 
                     (* dt 2.0 t_143 c34 tmp1)
                     (* dt 2.0 t_1_d__))]
              [v3 0.0]
              [v4 0.0])
          (DIAG0 d A v1 V2 V3 V4 0.0)))
  
      (define-syntax-rule (RX4M A C1C2I) (* dt2tmp2 (fr u (+ A ijk1)) C1C2I ))
      (define-syntax-rule (RX4L A C1C2I) (* (sqr (fr u (+ A ijk1))) C1C2I))

      (DIAG0 d 0 (+ 1.0 (* dt 2.0 tds1)) 0.0 0.0 0.0 0.0)

      (RX123 1 v2 v3 v4 v2 v3 v4 (+ (* tx1 r43) ty1 tz1) tds2)
      (RX123 2 v2 v3 v4 v3 v2 v4 (+ (* ty1 r43) tx1 tz1) tds3) 
      (RX123 3 v2 v3 v4 v3 v4 v2 (+ (* tz1 r43) tx1 ty1) tds4)

      (DIAG0 d 4 (* -2.0 dt (+ (* tmp3
                                  (RX4L 1 C1C2_1)
                                  (RX4L 2 C1C2_2)
                                  (RX4L 3 C1C2_3))
                               (* t_1s c1345 tmp1 (fr u (+ 4 ijk1)))))

        (RX4M 1 C1C2_1)
        (RX4M 2 C1C2_2)
        (RX4M 3 C1C2_3)
        (+ 1.0 
           (* dt2tmp1 t_1s c1345) 
           (* dt2 tds5))))

;;;//---------------------------------------------------------------------
;;;//   form the first,second,third block sub-diagonal
;;;//---------------------------------------------------------------------
  
  (DIAG a 3 1 2 C1 C1 C2 tz1 tz2 dz2 
    v2 v3 v4 
    v2 v3 v4
    (+ i (* j jsize3) (* (- k 1) ksize3))
    NEG NNEG)
  (DIAG b 2 1 3 C1 C2 C1 ty1 ty2 dy3 
    v2 v3 v4 
    v2 v4 v3
    (+ i (* (- j 1) jsize3) (* k ksize3))
    NEG NNEG)
  (DIAG c 1 2 3 C2 C1 C1 tx1 tx2 dx2 
    v2 v3 v4 
    v4 v2 v3 
    (+ (- i 1) (* j jsize3) (* k ksize3))
    NEG NNEG)
))))))


(define (jacu k)
  (DIAG a 1 2 3 C2 C1 C1 tx1 tx2 dx2 
    v2 v3 v4 
    v4 v2 v3 
    (+ (- i 1) (* j jsize3) (* k ksize3))
    NNEG NEG)
  (DIAG b 2 1 3 C1 C2 C1 ty1 ty2 dy3 
    v2 v3 v4 
    v2 v4 v3
    (+ i (* (- j 1) jsize3) (* k ksize3))
    NNEG NEG)
  (DIAG c 3 1 2 C1 C1 C2 tz1 tz2 dz2 
    v2 v3 v4 
    v2 v3 v4
    (+ i (* j jsize3) (* (- k 1) ksize3))
    NNEG NEG)
)

(define (l2norm ldx ldy ldz nx0 ny0 nzo ust iend jst jend V sum)
  (for ([m (in-range 5)]) (f! sum m 0.0))

  (for* ([k (in-range 1 (- nz0 1))]
         [j (in-range (sub1 jst) jend)]
         [i (in-range (sub1 ist) iend)])
    (f! sum m (for/fold ([x (fr sum m)]) ([m (in-range 5)])
      (sqr (fr v (+ m (* i isize1) (* j jsize1) (* k ksize1)))))))
  
  (let ([ns (* (- nx0 2) (- ny0 2) (- nz 2))])
    (for ([m (in-range 5)]) 
      (f! sum m (sqrt (fl/ (fr sum m)  ns))))))

(define (pintgr)
  (define phi1 (make-flvector (sqr ( + isiz3 2)) 0.0))
  (define phi2 (make-flvector (sqr ( + isiz3 2)) 0.0))

  (for* ([j (in-range (sub1 jbeg) jfin)]
         [i (in-range (sub1 ibeg) ifin)])
    (let ([ij (+ i (* j isize5))]
          [ij_ (+ (* i isize1) (* j isize1))])
      (define_syntax-rule (kit PHI KVAL)
        (let ([k (sub1 KVAL)]
              [ijk (+ ij_ (* k ksuze1))])
          (f! PHI  ij (- (* c2 (fr u (+ 4 ijk)))
                         (* 0.5 (/ (+ (sqr (fr u (+ 1 ijk)))
                                      (sqr (fr u (+ 2 ijk)))
                                      (sqr (fr u (+ 3 ijk))))
                                   (fr u (+ 0 ijk))))))))
      (kit ph1 (sub1 ki1))
      (kit ph2 (sub1 ki2))))

  (define frc1 (* dxi deta 
    (for*/fold ([frc1 0.0])  
          ([j (in-range (sub1 jbeg) jfin)]
           [i (in-range (sub1 ibeg) ifin)])
      (define-syntax-rule (stencil V)
        (+ (fr V (+ i (* j isize5)))
           (fr V (+ i 1 (* j isize5)))
           (fr V (+ i (* (+ j 1) isize5)))
           (fr V (+ i 1 (* (+ j 1) isize5)))))
      (+ frc1 (stencil phi1) (stencil phi2)))))



  (define-syntax-rule (FRCIT i j k TEST1 TEST2 I IBEG IFIN IJK1 IJK2 DDD)

    (for* ([i (in-range (+ isiz2 2))]
           [k (in-range (+ ksiz3 2))])
      (let ([ik (+ i (* k isize5))])
        (f! phi1 ik 0.0)
        (f! phi2 ik 0.0)))


    (define-syntax-rule (LOOP i j k TEST_ I_ PHI_ IJK_)
      (if TEST_
        (for* ([k (in-range (sub1 ki1) ki2)]
               [I (in-range (sub1 IBEG) IFIN)])
            (let (ijk IJK_)
            (f! PHI_  ij (- (* c2 (fr u (+ 4 ijk)))
                           (* 0.5 (/ (+ (sqr (fr u (+ 1 ijk)))
                                        (sqr (fr u (+ 2 ijk)))
                                        (sqr (fr u (+ 3 ijk))))
                                     (fr u (+ 0 ijk))))))))))

    (LOOP i j k TEST1 phi1 IJK1)
    (LOOP i j k TEST2 phi2 IJK2)

    (* DDD dzeta (for*/fold ([X 0.0]) 
          ([k (in-range (sub1 ki1) (sub1 ki2))]
           [i (in-range (sub1 IBEG) IFIN)])
      (define-syntax-rule (stencil V)
        (+ (fr V (+ i (* k isize5)))
           (fr V (+ i 1 (* k isize5)))
           (fr V (+ i (* (+ k 1) isize5)))
           (fr V (+ i 1 (* (+ k 1) isize5)))))
      (+ X (stencil phi1) (stencil phi2)))))

  (define frc2 (FRCIT i J k (= jbeg ji1) (= jfin ji2) i ibeg ifin
      (* i isize1) (* (- jbeg 1) jsize1) (* k ksize1)
      (* i isize1) (* (- jfin 1) jsize1) (* k ksize1) dxi))

  (define fr3 (FRCIT i j k (= ibeg ii1) (= ifin ii2) jbeg jfin
      (* (- ibeg 1) isize1) (* j jsize1) (* k ksize1)
      (* (- ifin 1) isize1) (* j jsize1) (* k ksize1) deta))

  (* 0.25 (+ frc1 frc2 frc3)))

(define (set-boundary-variables)
  (define-syntax-rule (INITIT i j k II JJ NII NJJ
    E1 E2 IDX1 IDX2)
    (for* ([II (in-range NII)]
           [JJ (in-range NJJ)])
      E1
      E2
      (for ([m (in-range 5)])
        (f! u IDX1 (vr temp1 m))
        (f! u IDX2 (vr temp2 m)))))

  (INITIT i j k j i ny nx
    (exact (+ i 1) (+ j 1) 1 temp1)
    (exact (+ i 1) (+ j 1) nz temp2)
    (+ m (* i isize1) (* j jsize1) 0)
    (+ m (* i isize1) (* j jsize1) (* (- nz 1) ksize1)))
  (INITIT i j k k i nz nx
    (exact (+ i 1) 1 (+ k 1) temp1)
    (exact (+ i 1) ny (+ k 1) temp2)
    (+ m (* i isize1) 0 (* k ksize1))
    (+ m (* i isize1) (* (- ny 1) jsize1) (* k ksize1)))
  (INITIT i j k k k nz ny
    (exact 1 (+ j 1) (+ k 1) temp1)
    (exact nx (+ j 1) (+ k 1) temp2)
    (+ m 0 (* j jsize1) (* i isize1))
    (+ m (* (- nx 1) isize1) (* j jsize1) (* i isize1))))

(define (set-initial-values)
  (define ie_1jk (make-flvector 5 0.0))
  (define ie_i1k (make-flvector 5 0.0))
  (define ie_ij1 (make-flvector 5 0.0))
  (define ie_nx0jk (make-flvector 5 0.0))
  (define ie_iny0k (make-flvector 5 0.0))
  (define ie_ijnz0 (make-flvector 5 0.0))
  (let ([Pface (make-flvector (* 5 3 2) 0.0)])
    (for ([k (in-range 1 (- nz 1))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range 1 (- ny 1))])
          (let ([eta (* j dnym1)])
            (for ([i (in-range 1 (- nx 1))])
              (let ([xi (* i dnxm1)])
                  (exact 1 (+ j 1) (+ k 1) ie_1jk)
                  (exact nx0 (+ j 1) (+ k 1) ie_nx0jk)
                  (exact (+ i 1) 1 (+ k 1) ie_i1k)
                  (exact (+ i 1) ny0 (+ k 1) ie_iny0k)
                  (exact (+ i 1) (+ j 1) 1 ie_ij1)
                  (exact (+ i 1) (+ j 1) nz0 ie_ijnz0)
                (for ([m (in-range 5)])
                  (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
                        [pxi (+ (- 1.0 xi) * (fr ue1jk m) (* xi (fr ue_nx0jk m)))]
                        [peta (+ (- 1.0 eta) * (fr uei1k m) (* eta (fr ue_iny0k m)))]
                        [pzeta (+ (- 1.0 zeta) * (fr ueij1 m) (* zeta (fr ue_ijnz0 m)))])
                    (f! u idx (+ pxi peta pzeta
                                 (- (* pxi peta))
                                 (- (* peta pzeta))
                                 (- (* pzeta pxi))
                                 (* pxi pet pzeta)))))))))))))

(define (sssor)
  (for ([j (in-range isiz2)]
        [i (in-range isiz1)]
        [n (in-range 5)]
        [m (in-range 5)])
    (let ([idx (+ m (* n isize4) (*i jsize4) (* j ksize4))])
      (f! a idx 0.0)
      (f! b idx 0.0)
      (f! c idx 0.0)
      (f! d idx 0.0)))

  (rhs)

  (lnorm isiz1 isiz2 isiz3 nx0 ny0 nz0 ist iend jst jend)

  (timer-start 1)

  (for ([istep (in-range 1 (add itmax))])
    (when (or (zero? (modulo istep 20)) (= istep itmax) (= istep 1))
      (printf " Time step ~a\n" istep))

    (for ([k (in-range 1 (sub1 nz))]
          [j (in-range (sub1 jst) jend)]
          [j (in-range (sub1 jst) jend)]
          [m (in-range 5)])
      (f!* rsd (+ m (* i isize1) (*j jsize1) (* k ksize1)) dt))

    (for ([k (in-range 1 (sub1 nz))]) 
      (jacld k)
      (blts isiz1 isiz2 isiz3 nx ny nz k omega rsd tv a b c d ist iend jst jend nx0 ny0))

    (for ([k (in-range (- nz 2) 0 -1)]) 
      (jacu k)
      (buts isiz1 isiz2 isiz3 nx ny nz k omega rsd tv d a b c ist iend jst jend nx0 ny0))

    (for ([k (in-range 1 (sub1 nz))]
          [j (in-range (sub1 jst) jend)]
          [j (in-range (sub1 jst) jend)]
          [m (in-range 5)])
      (let ([idx (+ m (* i isize1) (*j jsize1) (* k ksize1))])
        (f!+ u idx (* tmp (fr rsd idx)))))

    (if (zero? (modulo istep inorm))
      (lnorm isiz1 isiz2 isiz3 nx0 ny0 nz0 ist iend jst jedn rsd rsdnm))

    (if (and ((fr rsdnm 0) . < . (fr tolrsd 0))
             ((fr rsdnm 1) . < . (fr tolrsd 1))
             ((fr rsdnm 2) . < . (fr tolrsd 2))
             ((fr rsdnm 3) . < . (fr tolrsd 3))
             ((fr rsdnm 4) . < . (fr tolrsd 4)))
      (timer-stop 1))))


(define (verify xcr xcell xci)
  (define-values (xcrref xceref dtref xciref)
    (case class
      [(#\S)
        (values
          (flvector
      1.6196343210976702E-2
      2.1976745164821318E-3
      1.5179927653399185E-3
      1.5029584435994323E-3
      3.4264073155896461E-2)
          (flvector
      6.4223319957960924E-4
      8.4144342047347926E-5
      5.8588269616485186E-5
      5.8474222595157350E-5
      1.3103347914111294E-3)

      0.5 
      7.8418928865937083)]
      [(#\W)
        (values
          (flvector
      0.1236511638192E+02
      0.1317228477799E+01
      0.2550120713095E+01
      0.2326187750252E+01
      0.2826799444189E+02)
          (flvector
      0.4867877144216
      0.5064652880982E-1
      0.9281818101960E-1
      0.8570126542733E-1
      0.1084277417792E+01)
      1.5E-3
      0.1161399311023E+02)]
      [(#\A)
        (values
          (flvector
      7.7902107606689367E+02
      6.3402765259692870E+01
      1.9499249727292479E+02
      1.7845301160418537E+02
      1.8384760349464247E+03)
          (flvector
      2.9964085685471943E+01
      2.8194576365003349
      7.3473412698774742
      6.7139225687777051
      7.0715315688392578E+01)
      2.0
      2.6030925604886277E+01)]
      [(#\B)
        (values
          (flvector
      1.03766980323537846E+04
      8.92212458801008552E+02
      2.56238814582660871E+03
      2.19194343857831427E+03
      1.78078057261061185E+04)
          (flvector
      2.15986399716949279E+02
      1.55789559239863600E+01
      5.41318863077207766E+01
      4.82262643154045421E+01
      4.55902910043250358E+02)
      2.0
      6.66404553572181300E+01)]
      [else
        (values
          (make-flvector 5 1.0)
          (make-flvector 5 1.0)
          0.0001
          0.0001)]))

  (for ([m (in-range 5)])
    (f! xcrdif m (abs (/ (- (fr xcr m) (fr xcrref m)) (fr xcrref m))))
    (f! xcedif m (abs (/ (- (fr xce m) (fr xceref m)) (fr xceref m)))))

  (define xcidif (abs (/ (- xci xciref) xciref)))
  (define epsilon 1.0E-8)

  (if (not (equal? class #\U))
    (begin
      (printf " Verification being performed for class ~a\n" class)
      (printf " Accuracy setting for epsilon = ~a\n" epsilon)
      (if ((abs (- dt dtref)) . <= . epsilon)
        #t
        (printf " DT= ~a does not match the reference value of ~a\n" dt dtref)))
    (printf " Unknown class\n"))

  (if (not (equal? class #\U))
    (printf " Comparison of RMS-norms of residual\n")
    (printf " RMS-norms of residual"))

  (for ([m (in-range (flvector-length xcr))])
    (printf "~a. ~a ~a ~a\n" m (flvr xcr m) (flvr xcrref m) (flvr xcrdif m)))

  (if (not (equal class #\U))
    (printf " Comparison of RMS-norms of solution error\n")
    (printf " RMS-norms of solution error"))

  (for ([m (in-range (flvector-length xce))])
    (printf "~a. ~a ~a ~a\n" m (flvr xce m) (flvr xceref m) (flvr xcedif m)))

  (if (not (equal class #\U))
    (printf " Comparison of surface integral")
    (printf " Surface integral"))

  (for ([m (in-range (flvector-length xci))])
    (printf "~a. ~a ~a ~a\n" m (flvr xci m) (flvr xciref m) (flvr xcidif m))))

    

  






