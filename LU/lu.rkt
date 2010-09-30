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

(define-syntax-rule (omega-worker m i j k isize4 ldz jkx ijk-)
  (* omega (+ (* (fr ldz (+ (* 0 isize4) jkx)) (fr v (+ 0 ijk-)))
              (* (fr ldz (+ (* 1 isize4) jkx)) (fr v (+ 1 ijk-)))
              (* (fr ldz (+ (* 2 isize4) jkx)) (fr v (+ 2 ijk-)))
              (* (fr ldz (+ (* 3 isize4) jkx)) (fr v (+ 3 ijk-)))
              (* (fr ldz (+ (* 4 isize4) jkx)) (fr v (+ 4 ijk-))))))

(define-syntax-rule (bltw m i j k v BODY1 KK)
  (for ([j (in-range (sub1 jst) jend)]
        [i (in-range (sub1 ist) iend)]
        [m (in-range 5)])
    (let ([ijx  (+ m (* i isize1) (* j jsize1))]
          [jkx  (+ m (* i jsize4) (* j ksize4))]
          [ijk  (+ (* i isize1) (* j jsize1) (* k ksize1))]
          [ijk- (+ (* i isize1) (* j jsize1) (* (- k 1) ksize1))])
      (f! tv ijx (+ (fr v (+ m ijk))
                    (omega-worker m i j k isize4 ldz jkx ijk-)))))

  (for ([j (in-range (sub1 jst) jend)])
    (for ([i (in-range (sub1 ist) iend)])
      (for ([m (in-range 5)])
        (let ([ijx (+ m (* i isize1) (* j jsize1))]
              [idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
              [jkx (+ m (* i jsize4) (* j ksize4))]
              [i-x (+ (* (- i 1)  isize1) (* j jsize1) (* k ksize1))]
              [j-x (+ (* i isize1) (* (- j 1)  jsize1) (* k ksize1))])
          (f!- tv ijx (- (* omega (+ 
                                     (* (fr ldy (+ (* 0 isize4) jkx)) (fr v (+ 0 j-x)))
                                     (* (fr ldy (+ (* 1 isize4) jkx)) (fr v (+ 1 j-x)))
                                     (* (fr ldy (+ (* 2 isize4) jkx)) (fr v (+ 2 j-x)))
                                     (* (fr ldy (+ (* 3 isize4) jkx)) (fr v (+ 3 j-x)))
                                     (* (fr ldy (+ (* 4 isize4) jkx)) (fr v (+ 4 j-x)))
                                     (* (fr ldx (+ (* 0 isize4) jkx)) (fr v (+ 0 i-x)))
                                     (* (fr ldx (+ (* 1 isize4) jkx)) (fr v (+ 2 i-x)))
                                     (* (fr ldx (+ (* 2 isize4) jkx)) (fr v (+ 3 i-x)))
                                     (* (fr ldx (+ (* 3 isize4) jkx)) (fr v (+ 4 i-x)))
                                     (* (fr ldx (+ (* 4 isize4) jkx)) (fr v (+ 5 i-x)))))))))
      
      (let ([ijx (+ (* i isize1) (* j jsize1))]
            [jkx (+ (* i jsize4) (* j ksize4))])
        (for ([m (in-range 5)])
          (define-syntax-rule (doit I1 I2) (f! tmat (+ m I1) (fr d (+ m (* I2 isize4) jkx)))) 
          (doit 0  0) 
          (doit 5  1)
          (doit 10 2)
          (doit 15 3)
          (doit 20 4))

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

      (f! v (+ 4 ijkx) (/ (fr tv (+ 4 ijx) (fr tmat (+ 4 20)))))
      (f!- tv (+ 3 ijx) (* (fr tmat (+ 3 20)) (fr v (+ 4 ijkx))))

      (f! v (+ 3 ijkx) (/ (fr tv (+ 3 ijx) (fr tmat (+ 3 15)))))
      (f!- tv (+ 2 ijx) (* (fr tmat (+ 2 15)) (fr v (+ 3 ijkx)))
                        (* (fr tmat (+ 2 20)) (fr v (+ 4 ijkx))))

      (f! v (+ 2 ijkx) (/ (fr tv (+ 2 ijx) (fr tmat (+ 2 10)))))
      (f!- tv (+ 1 ijx) (* (fr tmat (+ 1 10)) (fr v (+ 2 ijkx)))
                        (* (fr tmat (+ 1 15)) (fr v (+ 3 ijkx)))
                        (* (fr tmat (+ 1 20)) (fr v (+ 4 ijkx))))

      (f! v (+ 1 ijkx) (/ (fr tv (+ 1 ijx) (fr tmat (+ 1  5)))))
      (f!- tv (+ 0 ijx) (* (fr tmat (+ 0  5)) (fr v (+ 1 ijkx)))
                        (* (fr tmat (+ 0 15)) (fr v (+ 2 ijkx)))
                        (* (fr tmat (+ 0 15)) (fr v (+ 3 ijkx)))
                        (* (fr tmat (+ 0 20)) (fr v (+ 4 ijkx))))

      (f! v (+ 0 ijkx) (/ (fr tv (+ 0 ijx) (fr tmat (+ 0  0)))))))))

(define (blts ldmx ldmy ldmz nx ny nz k omega v tv ldz ldy ldx d ist iend jst jend nx0 ny0)
  (bltsw i j k m v
    (- (fr v (IDX m i j k)) (bltw1 m i j k imega isize4 ldz jkx ijk)) (- k 1)
(define (buts ldmx ldmy ldmz nx ny nz k omega v tv ldz ldy ldx d ist iend jst jend nx0 ny0)
  (bltsw i j k m v (bltw1 m i j k imega isize4 udz jkx ijk)) (+ k 1)

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
(define-syntax-rule (flux-flow)
  (for ([k (in-range 1 (sub nz))])
    (for ([j (in-range (sub1 jst) jend)])
      (for ([i (in-range nx)])
        (let* ([ijkx (+ (* i isize1) (* j jsize1) (* k ksize1))]
               [iix (* ii isize2)]
               [u21 (/ (fr rsd (+ 1 ijkx)) (fr rsd (+ 0 ijkx)))]
               [q (* 0.50 (/ (+ (sqr (fr rsd (+ 1 ijkx)))
                                 (sqr (fr rsd (+ 2 ijkx)))
                                 (sqr (fr rsd (+ 3 ijkx))))
                              (fr rsd (+ 0 ijkx))))])
          (f! flux (+ 0 iix) (* (fr rsd (+ 1 ijkx))))
          (f! flux (+ 1 iix) (+ (* (fr rsd (+ 1 ijkx)) u21)
                                    (* (- (fr rsd (+ 4 ijkx)) q) c2)))
          (f! flux (+ 2 iix) (* (fr rsd (+ 2 ijkx)) u21))
          (f! flux (+ 3 iix) (* (fr rsd (+ 3 ijkx)) u21))
          (f! flux (+ 4 iix) (* (- (* (fr rsd (+ 4 ijkx)) c1) (* 2 q)) uw1))))
      (for ([i (in-range (sub1 ist) iend)])
        (for ([m (in-range 5)])
          (f!- frct (* tx2 (- (fr flux (+ m (* ip1 isize2)))
                              (fr flux (+ m (* im1 isize2))))))))

      (for ([i (in-range (sub1 ist) nx)])
        (let* ([ijkx (+ (* i isize1) (* j jsize1) (* k ksize1))]
               [ijkx (+ (* (- i 1)  isize1) (* j jsize1) (* k ksize1))]
               [iix (* ii isize2)]
               [tmp (/ 1.0 (fr rsd (+ 0 ijkx)))]
               [u2li (* tmp (fr rsd (+ 1 ijkx)))]
               [u3li (* tmp (fr rsd (+ 2 ijkx)))]
               [u4li (* tmp (fr rsd (+ 3 ijkx)))]
               [u5li (* tmp (fr rsd (+ 4 ijkx)))]
               [tmp2 (/ 1.0 (fr rsd (+ 0 ijkx)))]
               [u2lim1 (* tmp2 (fr rsd (+ 1 i-jkx)))]
               [u3lim1 (* tmp2 (fr rsd (+ 2 i-jkx)))]
               [u4lim1 (* tmp2 (fr rsd (+ 3 i-jkx)))]
               [u5lim1 (* tmp2 (fr rsd (+ 4 i-jkx)))])
          (f! flux (+ 1 iix) (* (/ 4.0 3.0) tx3 (- u2li u2lim1)))
          (f! flux (+ 2 iix) (* tx3 (- u3li u3lim1)))
          (f! flux (+ 3 iix) (* tx3 (- u4li u4lim1)))
          (f! flux (+ 4 iix) (+ (* 0.50 (- 1.0 (* c1 c5))
                                   tx3 (- (+ (sqr u2li) (sqr u3li) (sqr u4li))
                                          (+ (sqr u2lim1) (sqr u3lim1) (sqr u4lim1))))
                                (* (/ 1.0 6.0) tx3 (- (sqr u2li) (sqr u2lim1)))
                                (* c1 c5 tx3 (- u5li u5lim1))))))

      (for ([i (in-range (sub1 ist) nx)])
        (let* ([ijkx (+ (* i isize1) (* j jsize1) (* k ksize1))])
            (define-syntax-rule (citA C C1 C2 C3) (* C (+ (- C1 (* 2.0 C2)) C3)))
            (define-syntax-rule (cit3S C V I1 I2 I3) (citA C (flvr V I1) (flvr V I2) (flvr V I3)))
            (define-syntax-rule (mid d__t_1 AA)
              (let ([idxA (+ idx AA)]
                    [idxz+A (+ idxz+ AA)]
                    [idxz-A (+ idxz- AA)])
                (flvs!+ fct (+ ijkx AA)
                  (cit3S d__t_1 rsd idxz+A idxA idxz-A)
                  (* t_3 c3 c4 (- (fr flux (+ AA (iip1 isize2)))
                                  (fr flux (+ AA (ii isize2))))))))

            (flvs!+ frct ijkx (cit3S d_1t_1 rsd idxz+ idx idxz-))
            (mid d_2t_1 1)
            (mid d_3t_1 2)
            (mid d_4t_1 3)
            (mid d_5t_1 4)))
      
      (fourth-order-dissipation frct rsd ZMOST ZS NZ))))

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









