#lang racket/base

(provide main)
  
(require "bm-args.rkt") 
(require "bm-results.rkt") 
(require "rand-generator.rkt")
(require "timer.rkt")
(require "parallel-utils.rkt")
(require racket/match)
(require racket/math)
(require "macros.rkt")
(require "debug.rkt")
(require (for-syntax "macros.rkt"))
(require (for-syntax racket/base))
(require (for-meta 2 racket/base))

(require (only-in racket/flonum make-flvector make-shared-flvector flvector-length flvector))
#|
(require (rename-in racket/flonum
                    [flvector-ref fr] 
                    [flvector-set! f!]))

|#

#|
|#
(require (rename-in racket/unsafe/ops
                    [unsafe-flvector-ref fr] 
                    [unsafe-flvector-set! f!]
))
 
(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values  12  0.5    50)]
    [(#\W) (values  33  0.0015 300)]
    [(#\A) (values  64  2      250)]
    [(#\B) (values  102 2      250)]
    [(#\C) (values  250 2      250)]
    [else (raise "Unknown class")]))

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


(define (run-benchmark args) 
  (define maxlevel 11)
  (let ([bmname "LU"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(problem_size dt_default niter_default) (get-class-size CLASS)])
    (let* (
          [isiz1 problem_size]
          [isiz2 problem_size]
          [isiz3 problem_size]
          [dt dt_default]
          [niter niter_default]
          [isize1 5]
          [jsize1 (* 5 (add1 isiz1))]
          [ksize1 (* 5 (add1 isiz1) (add1 isiz2))]
          [isize2 5]
          [jsize3 (add1 isiz1)]
          [ksize3 (* (add1 isiz1) (add1 isiz2))]
          [isize4 5]
          [jsize4 (* 5 5)]
          [ksize4 (* 5 5 (add1 isiz1))]
          [s1 (* ksize1 isiz3)]
          [s2 (* ksize3 isiz3)]
          [s3 (* ksize4 isiz2)]
          [u (make-shared-flvector s1 0.0)]
          [rsd (make-shared-flvector s1 0.0)]
          [frct (make-shared-flvector s1 0.0)]

          [qs (make-shared-flvector s2 0.0)]
          [rho_i (make-shared-flvector s2 0.0)]

          [a (make-shared-flvector s3 0.0)]
          [b (make-shared-flvector s3 0.0)]
          [c (make-shared-flvector s3 0.0)]
          [d (make-shared-flvector s3 0.0)]
          [tv (make-shared-flvector (* 5 (+ isiz1 1) isiz2) 0.0)]

          [flux (make-shared-flvector (* 5 isiz1)  0.0)]
          [rsdnm (make-shared-flvector 5 0.0)]
          [tolrsd (make-shared-flvector 5 0.0)]
          [errnm (make-shared-flvector 5 0.0)]
        
          [nx isiz1]
          [ny isiz2]
          [nz isiz3]

          [c1 1.4]
          [c2 0.4]
          [c3 0.1]
          [c4 1.0]
          [c5 1.4]
          [c1c5 (* c1 c5)]
          [c3c4 (* c3 c4)]
          [r43 (/ 4.0 3.0)]
          [c1345 (* c1 c3 c4 c5)]
          [c34 (* c3 c4)]
          [omega 1.2]

          [dnxm1 (/ 1.0  (- nx 1))]
          [dnym1 (/ 1.0  (- ny 1))]
          [dnzm1 (/ 1.0  (- nz 1))]
          [dxi (/ 1.0  (- nx 1))]
          [deta (/ 1.0  (- ny 1))]
          [dzeta (/ 1.0  (- nz 1))]
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
          [dx1tx1 (* dx1 tx1)]
          [dx2tx1 (* dx2 tx1)]
          [dx3tx1 (* dx3 tx1)]
          [dx4tx1 (* dx4 tx1)]
          [dx5tx1 (* dx5 tx1)]
          [dy1ty1 (* dy1 ty1)]
          [dy2ty1 (* dy2 ty1)]
          [dy3ty1 (* dy3 ty1)]
          [dy4ty1 (* dy4 ty1)]
          [dy5ty1 (* dy5 ty1)]
          [dz1tz1 (* dz1 tz1)]
          [dz2tz1 (* dz2 tz1)]
          [dz3tz1 (* dz3 tz1)]
          [dz4tz1 (* dz4 tz1)]
          [dz5tz1 (* dz5 tz1)]
;;;//---------------------------------------------------------------------
;;;//   fourth difference dissipation
;;;//---------------------------------------------------------------------
          [dssp (/ (max dx1 dy1 dz1) 4.0)])
 

      ;(get-input-pars)
      (domain nx ny nz nx ny nz)
      
      (print-banner "LU" args) 
      (printf "~a: Iterations=~a dt=~a\n" bmname niter dt)

      ;(set-coefficients)
      (set-boundary-variables u nx ny nz isize1 jsize1 ksize1)
      (set-initial-values u nx ny nz dnzm1 dnym1 dnxm1 isize1 jsize1 ksize1)

      (compute-forcing-term frct rsd nz ny nx dnzm1 dnym1 dnxm1 isize1 jsize1 ksize1)
      (compute-forcing-term-exact flux rsd frct nz ny nx
        c1 c2 c1c5 c3c4 r43 dssp
        isize1 jsize1 ksize1 isize2 jsize3 ksize3
        tx2 tx3 ty2 ty3 tz2 tz3
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
        dnzm1 dnym1 dnxm1)
  

      (CGspawn (if serial 0 num-threads) 
        sssor a b c d u rsd rsdnm tv nx ny nz niter niter dt omega isize1 jsize1 ksize1 isize4 jsize4 ksize4
        rho_i flux qs frct
        r43 c1345 c1 c2 c1c5 c3c4 dssp
        isize2 jsize3 ksize3
        tx1 tx2 tx3
        ty1 ty2 ty3
        tz1 tz2 tz3
        dx1 dy1 dz1
        dx2 dy3 dz2
        dx5 dy5 dz5
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

      (compute-error errnm u nx ny nz isize1 jsize1 ksize1)

      (let* ([verified (verify CLASS rsdnm errnm (compute-surface-integral u nx ny nz isiz1 isiz2 isiz3 isize1 jsize1 ksize1 c2 dxi deta dzeta) dt
)])

        (print-verification-status CLASS verified bmname)
        (let* ([time (/ (read-timer 1) 1000)]
               [results (new-BMResults bmname CLASS nx ny nz niter time 
                                       (get-mflops time niter nx ny nz)
                                       "floating point" 
                                       (if verified 1 0)
                                       serial 
                                       num-threads 
                                       -1)]) 
            (print-results results)))))))

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

; ist = 2
; jst = 2
; iend = nx -1 
; jend = ny -1 
; ii1 = 2
; ji1 = 2
; ki1 = 3
; ii2 = nx -1
; ji2 = ny -2
; ki2 = nx -1
(define-syntax-rule (istR nx) (in-range 1 (sub1 nx))) 
(define-syntax-rule (jstR ny) (in-range 1 (sub1 ny))) 
(define-syntax-rule (kstR nz) (in-range 1 (sub1 nz))) 

(define-syntax-case (gen-b_ts DNAME KS)
  (with-syntax ([XM (PICK2 #'KS fx- fx+)]
                [iRR (PICK2 #'KS (in-range 1 (sub1 nx)) (in-range (- nx 2) 0 -1))]
                [jRR (PICK2 #'KS (p-range cg (in-range 1 (sub1 ny))) (p-range cg (in-range (fx- ny 2) 0 -1)))]
                [f!or- (PICK2 #'KS f! f!-)]
                [(MODIFIER ...) (PICK2 #'KS (fl- (fr v mijk)) (values) )])

  #'(define (DNAME cg tmatX v tv ldz ldy ldx d nx ny nz k omega
      isize1 jsize1 ksize1 isize4 jsize4 ksize4)
  (define tmat (make-flvector (fx* 5 5) 0.0))
    (for ([j jRR])
      (for  ([i iRR])
      (let* ([ij1  (fx+ (fx* i isize1) (fx* j jsize1))]
             [ijk1 (fx+ ij1 (fx* k ksize1))]
             [jk4 (fx+ (fx* i jsize4) (fx* j ksize4))])
      (for ([m (in-range 5)])
        (let* ([mij1 (fx+ m ij1)]
               [mjk4 (fx+ m jk4)]
               [mijk (fx+ m ijk1)]
               [i-jk (XM ijk1 isize1)]
               [ij-k (XM ijk1 jsize1)]
               [ijk- (XM ijk1 ksize1)])
          (f! tv mij1 (MODIFIER ...
              (fl* omega
                (for/fold ([T 0.0]) ([n (in-range 5)]
                                     [n5 (in-range 0 21 5)])
                  (let ([nmjk4 (fx+ (fx* n isize4) mjk4)])
                    (f! tmat (fx+ m n5) (fr d nmjk4))
;                    (printf "~a ~a ~a ~a ~a ~a\n" (CGid cg) k j i m n) (flush-output)
                    (fl+ T
                       (fl* (fr ldz nmjk4) (fr v (fx+ n ijk-)))
                       (fl* (fr ldy nmjk4) (fr v (fx+ n ij-k)))
                       (fl* (fr ldx nmjk4) (fr v (fx+ n i-jk)))))))))))

      (for ([C  (in-range 4)]
            [C5 (in-range 0 16 5)])
        (let ([tmp1 (/ 1.0 (fr tmat (fx+ C C5)))])
        (for ([A  (in-range (fx+ C 1) 5)])
          (let ([tmp (fl* tmp1 (fr tmat (fx+ A C5)))])
            (for ([N  (in-range (fx+ C 1) 5)]
                  [N5  (in-range (fx+ 5 C5) 21 5)])
              (f!- tmat (fx+ A N5) (fl* tmp (fr tmat (fx+ C N5)))))
            (f!- tv (fx+ A ij1) (fl* tmp (fr tv (fx+ C ij1))))))))
    
      (let* ([tv4 (fl/ (fr tv (fx+ 4 ij1)) (fr tmat 24))]
             [tv3 (fl/ (fl- (fr tv (fx+ 3 ij1)) 
                        (fl* (fr tmat 23) tv4)) (fr tmat 18))]
             [tv2 (fl/ (fl- (fr tv (fx+ 2 ij1)) 
                        (fl* (fr tmat 17) tv3)
                        (fl* (fr tmat 22) tv4)) (fr tmat 12))]
             [tv1 (fl/ (fl- (fr tv (fx+ 1 ij1)) 
                        (fl* (fr tmat 11) tv2)
                        (fl* (fr tmat 16) tv3)
                        (fl* (fr tmat 21) tv4)) (fr tmat 6))]
             [tv0 (fl/ (fl- (fr tv (fx+ 0 ij1)) 
                        (fl* (fr tmat  5) tv1)
                        (fl* (fr tmat 10) tv2)
                        (fl* (fr tmat 15) tv3)
                        (fl* (fr tmat 20) tv4)) (fr tmat 0))])
        (f!or- v (fx+ 0 ijk1) tv0)
        (f!or- v (fx+ 1 ijk1) tv1)
        (f!or- v (fx+ 2 ijk1) tv2)
        (f!or- v (fx+ 3 ijk1) tv3)
        (f!or- v (fx+ 4 ijk1) tv4))))))))

(gen-b_ts blts 1)
(gen-b_ts buts 2)

(define (domain nx ny nz isiz1 isiz2 isiz3)
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

(define-syntax-rule (fourth-order-dissipation V V2 ZMOST ZS NZ dssp)
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
               [didx1  m1]
               [didx2  m2])
          (f!- V didx1 (* dssp (+ (*  5.0 (fr V2 m1))
                                  (* -4.0 (fr V2 m2))
                                          (fr V2 m3))))
          (f!- V didx2 (* dssp (+ (* -4.0 (fr V2 m1))
                                  (*  6.0 (fr V2 m2))
                                  (* -4.0 (fr V2 m3))
                                          (fr V2 m4)))))))
    (for ([ii (in-range 3 (- NZ 3))])
      (let ([i-2 (+ (* (- ii 2) ZS) zmost)]
            [i-1 (+ (* (- ii 1) ZS) zmost)]
            [ix  (+ (* ii       ZS) zmost)]
            [i+1 (+ (* (+ ii 1) ZS) zmost)]
            [i+2 (+ (* (+ ii 2) ZS) zmost)])
        (for ([m (in-range 5)])
          (let* ([m-2 (+ m i-2)]
                 [m-1 (+ m i-1)]
                 [mx  (+ m ix)]
                 [m+1 (+ m i+1)]
                 [m+2 (+ m i+2)]
                 [didx  mx])
            (f!- V didx (* dssp (+         (fr V2 m-2)
                                   (* -4.0 (fr V2 m-1))
                                   (*  6.0 (fr V2 mx))
                                   (* -4.0 (fr V2 m+1))
                                           (fr V2 m+2))))))))
    (let ([i5 (+ (* (- NZ 5) ZS) zmost)]
          [i4 (+ (* (- NZ 4) ZS) zmost)]
          [i3 (+ (* (- NZ 3) ZS) zmost)]
          [i2 (+ (* (- NZ 2) ZS) zmost)])
      (for ([m (in-range 5)])
        (let* ([m5 (+ m i5)]
               [m4 (+ m i4)]
               [m3 (+ m i3)]
               [m2 (+ m i2)]
               [didx1  m3]
               [didx2  m2])
          (f!- V didx1 (* dssp (+         (fr V2 m5)
                                  (* -4.0 (fr V2 m4))
                                  (*  6.0 (fr V2 m3))
                                  (* -4.0 (fr V2 m2)))))
          (f!- V didx2 (* dssp (+         (fr V2 m4)
                                  (* -4.0 (fr V2 m3))
                                  (*  5.0 (fr V2 m2))))))))))
(define (compute-forcing-term frct rsd nz ny nx dnzm1 dnym1 dnxm1 isize1 jsize1 ksize1)

  (for ([m (in-range (flvector-length frct))]) (f! frct m 0.0))

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
                                                 (* zeta (fr ce (+ m (* 12 5))))))))))))))))))))

(define-syntax-case (flux-flow cg UD A flux RSD/U QS RHO FRCT nz ny nx
  c1 c2 c1c5 c3c4 r43 dssp
  isize1 jsize1 ksize1 isize2 jsize3 ksize3
  t_2 t_3 d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1)
(with-syntax-values ([(kk jj ii) (PICK3 #'A (k j i) (k i j) (j i k))]
                     [(nkk njj nii nii2) (PICK3 #'A ((fx-- nz) (fx-- ny) nx (fx-- nx))
                                                    ((fx-- nz) (fx-- nx) ny (fx-- ny))
                                                    ((fx-- ny) (fx-- nx) nz (fx-- nz)))]
                     [(ZMOST ZS NZ) (PICK3 #'A ((fx+ (fx* j jsize1) (fx* k ksize1)) isize1 nx)
                                               ((fx+ (fx* i isize1) (fx* k ksize1)) jsize1 ny)
                                               ((fx+ (fx* i isize1) (fx* j jsize1)) ksize1 nz))]
                     [(ZMOST3 ZS3)  (PICK3 #'A ((fx+ (fx* j jsize3) (fx* k ksize3)) 1)
                                               ((fx+ i (fx* k ksize3)) jsize3)
                                               ((fx+ i (fx* j jsize3)) ksize3))]
                     [((FM1 ...) (FM2 ...) (FM3 ...)) (REPDET* #'A (values) (values) (values) (fl* r43) (fl* r43) (fl* r43))])
  (with-syntax-values
                     ([(RHOI RHOI-) (PICK2 #'UD (ijk1z ijk1z-) ((fx+ (fx* ii ZS3) ZMOST3)
                                                                (fx+ (fx* (fx- ii 1) ZS3) ZMOST3)))])
  #'(let ()
  (define-syntax-rule (RHSU a (... ...)) (PICK2M UD (void) (begin a (... ...))))
  (define-syntax-rule (oneo a (... ...)) (PICK2M UD (/ 1.0 a (... ...)) (begin a (... ...))))
  (define-syntax-case (RETDET R2 N1 N2 N3 V1 V2 V3 BODY (... ...))
    (with-syntax ([LETBINDINGS (case (syntax->datum #'R2)
                            [(1) #'([N1 V1] [N2 V2] [N3 V3])]
                            [(2) #'([N1 V2] [N2 V1] [N3 V3])]
                            [(3) #'([N1 V3] [N2 V2] [N3 V1])]
                            [else (raise (format "invalid RETDET value ~a" #'R2))])])
    #'(let LETBINDINGS
      BODY (... ...)
    )))

  (CGfor cg ([kk (in-range 1 nkk)])
    (for ([jj (in-range 1 njj)])
      (let ([zmost ZMOST])

      (for ([ii (in-range 0 nii)])
        (let* ([ijk1z (fx+ zmost (fx* ii ZS))]
               [ijk3 (fx+ i (fx* j jsize3) (fx* k ksize3))]
               [iix (fx* ii isize2)]
               [rhotmp (fl/ 1.0 (fr RSD/U ijk1z))]
               [u21 (fl* (fr RSD/U (fx+ A ijk1z)) rhotmp)]
               [q (fl* 0.50 (fl* (fl+ (sqr (fr RSD/U (fx+ 1 ijk1z)))
                                (sqr (fr RSD/U (fx+ 2 ijk1z)))
                                (sqr (fr RSD/U (fx+ 3 ijk1z))))
                              rhotmp))])

          (RETDET A A1 A2 A3 1 2 3
          (RETDET A N2 N3 N4
            (fl+ (fl* (fr RSD/U (fx+ A1 ijk1z)) u21)
               (fl* c2 (fl- (fr RSD/U (fx+ 4 ijk1z)) q)))
            (fl* (fr RSD/U (fx+ A2 ijk1z)) u21)
            (fl* (fr RSD/U (fx+ A3 ijk1z)) u21)

            (f! flux (fx+ 0 iix)    (fr RSD/U (fx+ A ijk1z)))
            (f! flux (fx+ 1 iix) N2)
            (f! flux (fx+ 2 iix) N3)
            (f! flux (fx+ 3 iix) N4)
            (f! flux (fx+ 4 iix) (fl* u21 (fl- (fl* c1 (fr RSD/U (fx+ 4 ijk1z))) 
                                         (fl* c2 q))))))))
      (for ([ii (in-range 1 nii2)])
        (let ([ii+2 (fx* (fx+ ii 1) isize2)]
              [ii-2 (fx* (fx- ii 1) isize2)])
          (for ([m (in-range 5)])
            (let ([mijk1 (fx+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1))])
              (f!- FRCT mijk1 (fl* t_2 (fl- (fr flux (fx+ m ii+2))
                                        (fr flux (fx+ m ii-2)))))))))

      (for ([ii (in-range 1 nii)])
        (let* ([ijk1z (fx+ zmost (fx* ii ZS))]
               [ijk1z- (fx+ zmost (fx* (fx- ii 1) ZS))]
               [iix (fx* ii isize2)]
               [tmp (oneo (fr RHO RHOI))]
               [u2li (fl* tmp (fr RSD/U (fx+ 1 ijk1z)))]
               [u3li (fl* tmp (fr RSD/U (fx+ 2 ijk1z)))]
               [u4li (fl* tmp (fr RSD/U (fx+ 3 ijk1z)))]
               [u5li (fl* tmp (fr RSD/U (fx+ 4 ijk1z)))]
               [tmp2 (oneo (fr RHO RHOI-))]
               [u2lim1 (fl* tmp2 (fr RSD/U (fx+ 1 ijk1z-)))]
               [u3lim1 (fl* tmp2 (fr RSD/U (fx+ 2 ijk1z-)))]
               [u4lim1 (fl* tmp2 (fr RSD/U (fx+ 3 ijk1z-)))]
               [u5lim1 (fl* tmp2 (fr RSD/U (fx+ 4 ijk1z-)))])
          (f! flux (fx+ 1 iix) (FM1 ... (fl* t_3 (fl- u2li u2lim1))))
          (f! flux (fx+ 2 iix) (FM2 ... (fl* t_3 (fl- u3li u3lim1))))
          (f! flux (fx+ 3 iix) (FM3 ... (fl* t_3 (fl- u4li u4lim1))))
          (f! flux (fx+ 4 iix) (fl+ (fl* 0.50 (fl- 1.0 c1c5)
                                   t_3 (fl- (fl+ (sqr u2li) (sqr u3li) (sqr u4li))
                                          (fl+ (sqr u2lim1) (sqr u3lim1) (sqr u4lim1))))
                                (fl* (fl/ 1.0 6.0) t_3 (fl- (sqr (PICK3M A u2li u3li u4li)) 
                                          (sqr (PICK3M A u2lim1 u3lim1 u4lim1))))
                                (fl* c1c5 t_3 (fl- u5li u5lim1))))))

      (for ([ii (in-range 1 nii2)])
        (let* ([ijk1z (fx+ zmost (fx* ii ZS))]
               [iix (fx* ii isize2)]
               [iix+ (fx* (add1 ii) isize2)]
               [ijk1z+ (fx+ zmost (fx* (fx+ ii 1) ZS))]
               [ijk1z- (fx+ zmost (fx* (fx- ii 1) ZS))])

            (define-syntax-rule (citA C C1 C2 C3) (fl* C (fl+ (fl- C1 (fl* 2.0 C2)) C3)))
            (define-syntax-rule (cit3S C V I1 I2 I3) (citA C (fr V I1) (fr V I2) (fr V I3)))
            (define-syntax-rule (mid d__t_1 AA)
              (let ([idxA (fx+ ijk1z AA)]
                    [idxz+A (fx+ ijk1z+ AA)]
                    [idxz-A (fx+ ijk1z- AA)])
                (f!+ FRCT idxA
                  (cit3S d__t_1 RSD/U idxz+A idxA idxz-A)
                  (fl* t_3 c3c4 (fl- (fr flux (fx+ AA iix+))
                                     (fr flux (fx+ AA iix)))))))

            (f!+ FRCT ijk1z (cit3S d_1t_1 RSD/U ijk1z+ ijk1z ijk1z-))
            (mid d_2t_1 1)
            (mid d_3t_1 2)
            (mid d_4t_1 3)
            (mid d_5t_1 4)))
      
      (fourth-order-dissipation FRCT RSD/U zmost ZS NZ dssp))))))))

(define (compute-forcing-term-exact flux rsd frct nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    tx2 tx3 ty2 ty3 tz2 tz3
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
    dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
    dnzm1 dnym1 dnxm1)

  (define cgs (CGSingle))
  (flux-flow cgs 1 1 flux rsd   qs rsd frct nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    tx2 tx3 dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)
  (flux-flow cgs 1 2 flux rsd   qs rsd frct nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    ty2 ty3 dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)
  (flux-flow cgs 1 3 flux rsd   qs rsd frct nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    tz2 tz3 dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1))

;(define (act flux rsd        frct nz ny nx
;             flux rsd qs rsd frct nz ny nx
(define (rhs cg flux u   qs rho_i rsd frct nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    tx2 tx3 ty2 ty3 tz2 tz3
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
    dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)
  (for ([k (in-range nz)])
    (for ([j (in-range ny)])
      (for ([i (in-range nx)])
        (let ([idx (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))])
        (for ([m (in-range 5)])
          (let ([midx (fx+ m idx)])
          (f! rsd midx (fl- (fr frct midx)))))
        (let* ([ijk1z idx]
               [ijk3 (fx+ i (fx* j jsize3) (fx* k ksize3))]
               [rhotmp (fl/ 1.0 (fr u ijk1z))]
               [q (fl* 0.50 (fl* (fl+ (sqr (fr u (fx+ 1 ijk1z)))
                                (sqr (fr u (fx+ 2 ijk1z)))
                                (sqr (fr u (fx+ 3 ijk1z))))
                              rhotmp))])
          (f! qs ijk3 q)
          (f! rho_i ijk3 rhotmp))))))

  (CG-B cg)  

  (flux-flow cg 2 1 flux u qs rho_i rsd nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    tx2 tx3 dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)

  (CG-B cg)  

  (flux-flow cg 2 2 flux u qs rho_i rsd nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    ty2 ty3 dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)

  (CG-B cg)  

  (flux-flow cg 2 3 flux u qs rho_i rsd nz ny nx
    c1 c2 c1c5 c3c4 r43 dssp
    isize1 jsize1 ksize1 isize2 jsize3 ksize3
    tz2 tz3 dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1))

(define (exactKBT i j k u0 nx ny nz)
 (let ([xi   (/ (- i 1) (- nx 1))]
       [eta  (/ (- j 1) (- ny 1))]
       [zeta (/ (- k 1) (- nz 1))])
   (for ([m (in-range 5)])
     (f! u0 m 
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
                                     (* zeta (fr ce (+ m (* 12 5)))))))))))))))

(define (compute-error errnm u nx ny nz isize1 jsize1 ksize1)
  (define u000ijk (make-flvector 5 0.0))
  (define ns-2 (/ 1.0 (* (- nx 2) (- ny 2) (- nz 2))))
  (for ([m (in-range 5)]) (f! errnm m 0.0))

  (for* ([k (kstR nz)]
         [j (jstR ny)]
         [i (istR nx)])
    (exactKBT (add1 i) (add1 j) (add1 k) u000ijk nx ny nz)
    (for ([m (in-range 5)])
      (let* ([mijk1 (+ m (* i isize1) (* j jsize1) (* k ksize1))])
        (f!+ errnm m (sqr (- (fr u000ijk m) (fr u mijk1)))))))

  (for ([m (in-range 5)]) 
    (f! errnm m (sqrt (* (fr errnm m) ns-2)))))

(define-syntax-rule (let-with-syntax ([a b ...] ...) c ...)
  (let-syntax ([a (make-set!-transformer (lambda (stx) b ...))] ...)
        (begin c ...)))

(define-syntax-case (MIXER R2 MACRO V A V0 V1 V2 V3 V4)
  (let ([r1 (syntax->datum (local-expand #'R2 'expression #f))])
  (with-syntax ([X1 (case r1 [(1) #'V1] [(2) #'V2] [(3) #'V3][else (raise (format "invalid MIXER X1 value ~a" #'R2))])]
                [X2 (case r1 [(1) #'V2] [(2) #'V1] [(3) #'V2][else (raise (format "invalid MIXER X1 value ~a" #'R2))])]
                [X3 (case r1 [(1) #'V3] [(2) #'V3] [(3) #'V1][else (raise (format "invalid MIXER X1 value ~a" #'R2))])])
  #'(MACRO V A V0 X1 X2 X3 V4))))

(define-syntax-case (MIXERN R2 C2)
  (let ([r1 (syntax->datum #'R2)]
        [c1 (syntax->datum #'C2)])
    (case r1 [(1) (case c1 [(1) #'1] [(2) #'2] [(3) #'3][else (raise (format "invalidA RETDET value ~a" #'R2))])]
             [(2) (case c1 [(1) #'2] [(2) #'1] [(3) #'3][else (raise (format "invalidA RETDET value ~a" #'R2))])]
             [(3) (case c1 [(1) #'3] [(2) #'2] [(3) #'1][else (raise (format "invalidA RETDET value ~a" #'R2))])][else (raise (format "invalidA RETDET value ~a" #'R2))])))

(define-syntax-case (MIXERN3 R2 C2 V1 V2 V3)
  (let ([r1 (syntax->datum #'R2)]
        [c1 (syntax->datum #'C2)])
    (case r1 [(1) (case c1 [(1) #'V1] [(2) #'V2] [(3) #'V3][else (raise (format "invalidA RETDET value ~a" #'R2))])]
             [(2) (case c1 [(1) #'V2] [(2) #'V1] [(3) #'V3][else (raise (format "invalidA RETDET value ~a" #'R2))])]
             [(3) (case c1 [(1) #'V3] [(2) #'V2] [(3) #'V1][else (raise (format "invalidA RETDET value ~a" #'R2))])][else (raise (format "invalidA RETDET value ~a" #'R2))])))

(define-syntax-case (RETDET R2 N1 N2 N3 V1 V2 V3 BODY ...)
  (with-syntax ([LETBINDINGS (case (syntax->datum #'R2)
                          [(1) #'([N1 V1] [N2 V2] [N3 V3])]
                          [(2) #'([N1 V2] [N2 V1] [N3 V3])]
                          [(3) #'([N1 V3] [N2 V2] [N3 V1])]
                          [else (raise (format "invalidA RETDET value ~a" #'R2))])])
  #'(let LETBINDINGS
    BODY ...
  )))
(define-syntax-case (REPDET R N1 N2 N3 V1 V2 V3 VV1 VV2 VV3 BODY (... ...))
  (with-syntax ([LETBINDINGS (case (syntax->datum #'R)
                          [(1) #'([N1 VV1] [N2 V2] [N3 V3])]
                          [(2) #'([N1 V1] [N2 VV2] [N3 V3])]
                          [(3) #'([N1 V1] [N2 V2] [N3 VV3])]
                          [else (raise "invalid REPDET value")])])
  #'(let LETBINDINGS
    BODY (... ...)
  )))


(define-syntax-case (define-jac DNAME)
  #'(define (DNAME cg a b c d u k rho_i qs nx ny c1 c2 c34 r43 c1345
              isize1 jsize1 ksize1 jsize3 ksize3 isize4 jsize4 ksize4
              dt
              tx1 ty1 tz1
              tx2 ty2 tz2
              dx1 dy1 dz1
              dx2 dy3 dz2
              dx5 dy5 dz5
              dx1tx1 dy1ty1 dz1tz1
              dx2tx1 dy2ty1 dz2tz1
              dx3tx1 dy3ty1 dz3tz1
              dx4tx1 dy4ty1 dz4tz1
              dx5tx1 dy5ty1 dz5tz1)

    (let ([C1 (fl- c34 c1345 )]
          [C2 (fl- (fl* r43 c34) c1345)])
      (define-syntax-rule (C1C2__ C1 C2 C3) (fl+ (fl* tx1 C1) (fl* ty1 C2) (fl* tz1 C3)))

      ;; CONSTANTS
      (let ([t_1s (fl+ tx1 ty1 tz1)]
            [tds1 (fl+ dx1tx1 dy1ty1 dz1tz1)]
            [tds2 (fl+ dx2tx1 dy2ty1 dz2tz1)]
            [tds3 (fl+ dx3tx1 dy3ty1 dz3tz1)]
            [tds4 (fl+ dx4tx1 dy4ty1 dz4tz1)]
            [tds5 (fl+ dx5tx1 dy5ty1 dz5tz1)]
            [C1C2_1 (C1C2__ C2 C1 C1)]
            [C1C2_2 (C1C2__ C1 C2 C1)]
            [C1C2_3 (C1C2__ C1 C1 C2)])

      (for ([j (p-range cg (in-range 1 (sub1 ny)))])
        (for ([i (istR nx)])
          (let* (
            ;; LOOP DEPENDENT VARIABLES
            [ijk1 (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
            [ijk3 (fx+    i         (fx* j jsize3) (fx* k ksize3))]
            [jk4 (fx+ (fx* i jsize4) (fx* j ksize4))]
            [i0jk4  jk4]
            [i1jk4  (fx+ (fx* 1 isize4) jk4)]
            [i2jk4  (fx+ (fx* 2 isize4) jk4)]
            [i3jk4  (fx+ (fx* 3 isize4) jk4)]
            [i4jk4  (fx+ (fx* 4 isize4) jk4)])

            (define-syntax-case (DIAG0 V A V1 V2 V3 V4 V5)
              #'(begin
                (f! V (fx+ A i0jk4) V1)
                (f! V (fx+ A i1jk4) V2)
                (f! V (fx+ A i2jk4) V3)
                (f! V (fx+ A i3jk4) V4)
                (f! V (fx+ A i4jk4) V5)))

            (define-syntax-case (ROTASN V A R1 V1 V2 V3 V4 V5)
              #'(MIXER R1 DIAG0 V A V1 V2 V3 V4 V5))

       
  ;;;//---------------------------------------------------------------------
  ;;;//   form the block daigonal
  ;;;//---------------------------------------------------------------------
            (let* ([tmp1 (fr rho_i ijk3)]
                   [tmp2 (sqr tmp1)]
                   [tmp3 (fl* tmp1 tmp2)]
                   [dt2 (fl* dt 2.0)]
                   [dt2tmp2 (fl* dt2 tmp2)]
                   [dt2tmp1 (* dt2 tmp1)])
         
              (define-syntax-rule (RX123 A t_1r43 t_1_d__)
                (ROTASN d A A 
                  (fl- (fl* dt 2.0 t_1r43 c34 tmp2 (fr u (fx+ A ijk1))))
                  (fl+ 1.0 
                     (fl* dt 2.0 t_1r43 c34 tmp1)
                     (fl* dt 2.0 t_1_d__))
                  0.0
                  0.0
                  0.0))
          
              (define-syntax-rule (RX4M A C1C2I) (fl* dt2tmp2 (fr u (fx+ A ijk1)) C1C2I ))
              (define-syntax-rule (RX4L A C1C2I) (fl* (sqr (fr u (fx+ A ijk1))) C1C2I))

              (DIAG0 d 0 (fl+ 1.0 (fl* dt 2.0 tds1)) 0.0 0.0 0.0 0.0)

              (RX123 1 (fl+ (fl* tx1 r43) ty1 tz1) tds2)
              (RX123 2 (fl+ (fl* ty1 r43) tx1 tz1) tds3) 
              (RX123 3 (fl+ (fl* tz1 r43) tx1 ty1) tds4)

              (DIAG0 d 4 
                (fl* -2.0 dt (fl+ (fl* tmp3
                                 (fl+ (RX4L 1 C1C2_1)
                                    (RX4L 2 C1C2_2)
                                    (RX4L 3 C1C2_3)))
                              (fl* t_1s c1345 tmp2 (fr u (fx+ 4 ijk1)))))
                (RX4M 1 C1C2_1)
                (RX4M 2 C1C2_2)
                (RX4M 3 C1C2_3)
                (fl+ 1.0 
                   (fl* dt2tmp1 t_1s c1345) 
                   (fl* dt2 tds5)))))))

;;;//---------------------------------------------------------------------
;;;//   form the first,second,third block sub-diagonal
;;;//---------------------------------------------------------------------
(define-syntax-case (DIAG UD R)
  (let ()
  (define-syntax-case (PICK2 R2 V1 V2) 
    #'(begin
;      (printf "PICK2 ~a ~a ~a\n" (syntax->datum R2) (syntax->datum V1) (syntax->datum V2))
      (case (syntax->datum R2) [(1) V1] [(2) V2] 
      [else (printf "PICK2 error ~a~n" (syntax->datum R2))])))
  (define-syntax-case (PICK3 R2 V1 V2 V3) 
    #'(begin 
      (case (syntax->datum R2) [(1) V1] [(2) V2] [(3) V3]
      [else (raise (format "~a" (syntax->datum R2)))])))
  (with-syntax ([+- (PICK2 #'UD #'fx- #'fx+)])
  (with-syntax ([V (PICK2 #'UD (PICK3 #'R #'c #'b #'a) (PICK3 #'R #'a #'b #'c))]
                [t_1 (PICK3 #'R #'tx1 #'ty1 #'tz1)]
                [t_2 (PICK3 #'R #'tx2 #'ty2 #'tz2)] 
                [d__ (PICK3 #'R #'dx2 #'dy3 #'dz2)] 
                [d_1 (PICK3 #'R #'dx1 #'dy1 #'dz1)] 
                [d_5 (PICK3 #'R #'dx5 #'dy5 #'dz5)] 
                [NEG (PICK2 #'UD #'- #'begin)]
                [NNEG (PICK2 #'UD #'begin #'- )]
                [ijkz3-1E (PICK3 #'R 
                           #'(fx+ (+- i 1) (fx* j jsize3) (fx* k ksize3))
                           #'(fx+ i (fx* (+- j 1) jsize3) (fx* k ksize3))
                           #'(fx+ i (fx* j jsize3) (fx* (+- k 1) ksize3)))]
                [ijkz1-1E (PICK3 #'R 
                           #'(fx+ (fx* (+- i 1) isize1) (fx* j jsize1) (fx* k ksize1))
                           #'(fx+ (fx* i isize1) (fx* (+- j 1) jsize1) (fx* k ksize1))
                           #'(fx+ (fx* i isize1) (fx* j jsize1) (fx* (+- k 1) ksize1)))]
                [CCC1 (PICK3 #'R #'C2 #'C1 #'C1)]
                [CCC2 (PICK3 #'R #'C1 #'C2 #'C1)]
                [CCC3 (PICK3 #'R #'C1 #'C1 #'C2)]
                [BB0 (PICK3 #'R #'1 #'2 #'3)]
                [BB1 (PICK3 #'R #'2 #'1 #'2)]
                [BB2 (PICK3 #'R #'3 #'3 #'1)])
  

    #'(let* (; CONSTATNS
           [dtt_1 (fl* dt t_1)]
           [dtt_2 (fl* dt t_2)]
           [dtt_1d__ (fl* dtt_1 d__)]
           [dtt_1d_1 (fl* dtt_1 d_1)])

      (for ([j (jstR ny)])
        (for ([i (istR nx)])
          (let* (
            ;; LOOP DEPENDENT VARIABLES
            [ijk1 (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
            [ijk3 (fx+    i           (fx* j jsize3) (fx* k ksize3))]
            [ijkz3-1 ijkz3-1E]
            [ijkz1-1 ijkz1-1E]
            [jk4 (fx+ (fx* i jsize4) (fx* j ksize4))]
            [i0jk4  jk4]
            [i1jk4  (fx+ (fx* 1 isize4) jk4)]
            [i2jk4  (fx+ (fx* 2 isize4) jk4)]
            [i3jk4  (fx+ (fx* 3 isize4) jk4)]
            [i4jk4  (fx+ (fx* 4 isize4) jk4)]

            [tmp1 (fr rho_i ijkz3-1)]
            [tmp2 (sqr tmp1)]
            [tmp3 (fl* tmp1 tmp2)]
            [uBB0 (fr u (fx+ BB0 ijkz1-1))]
            [uBB1 (fr u (fx+ BB1 ijkz1-1))]
            [uBB2 (fr u (fx+ BB2 ijkz1-1))]
            [U1 (fr u (fx+ 1 ijkz1-1))]
            [U2 (fr u (fx+ 2 ijkz1-1))]
            [U3 (fr u (fx+ 3 ijkz1-1))]
            [U4 (fr u (fx+ 4 ijkz1-1))]
            [QS (fr qs ijkz3-1)]
            [dtt_1tmp1 (fl* dtt_1 tmp1)]
            [dtt_2tmp1 (fl* dtt_2 tmp1)]
            [dtt_1tmp2 (fl* dtt_1 tmp2)]
            [dtt_2tmp2 (fl* dtt_2 tmp2)]
            [dtt_2tmp1c2 (fl* dtt_2tmp1 c2)]
            [dtt_1tmp2c34 (fl* dtt_1tmp2 c34)]
            [dtt_1tmp1c34 (fl* dtt_1tmp1 c34)]
            [dtt_2tmp1uBB0 (fl* dtt_2tmp1 uBB0)]
            [dtt_2tmp2uBB0 (fl* dtt_2tmp2 uBB0)]
            [dtt_1tmp2uBB0 (fl* dtt_1tmp2 uBB0)])

            (define-syntax-case (DIAG0 V A V1 V2 V3 V4 V5)
              #'(begin
                (f! V (fx+ A i0jk4) V1)
                (f! V (fx+ A i1jk4) V2)
                (f! V (fx+ A i2jk4) V3)
                (f! V (fx+ A i3jk4) V4)
                (f! V (fx+ A i4jk4) V5)))

            (define-syntax-case (ROTASN V A R1 V1 V2 V3 V4 V5)
              #'(MIXER R1 DIAG0 V A V1 V2 V3 V4 V5))

            (define-syntax-case (RX2/3 R1 R2)
              (with-syntax ([UU #'(MIXERN3 R1 R2 U1 U2 U3)]
                            [A #'(MIXERN3 R1 R2 1 2 3)])
                (with-syntax ([V0 #'(fl+ (fl* UU (fl+ (NNEG dtt_2tmp2uBB0) dtt_1tmp2c34)))]
                            [V1 #'(fl- (NEG dtt_2tmp1uBB0) dtt_1tmp1c34 dtt_1d__)]
                            [V2 #'(NEG (fl* dtt_2tmp1 UU))]
                            [V3 #'0.0]
                            [V4 #'0.0])
                  (let ([rR2 (syntax->datum #'R2)])
                  (case (syntax->datum #'R1)
                    [(1) (case rR2 [(2) #'(DIAG0 V A V0 V2 V1 V3 V4)] [(3) #'(DIAG0 V A V0 V2 V3 V1 V4)])]
                    [(2) (case rR2 [(2) #'(DIAG0 V A V0 V1 V2 V3 V4)] [(3) #'(DIAG0 V A V0 V3 V2 V1 V4)])]
                    [(3) (case rR2 [(3) #'(DIAG0 V A V0 V1 V3 V2 V4)] [(2) #'(DIAG0 V A V0 V3 V1 V2 V4)])])))))

            (ROTASN V 0 R (fl- dtt_1d_1) (NEG dtt_2) 0.0 0.0 0.0)
            (ROTASN V (MIXERN R 1) (MIXERN R 1)
              (fl+ (NEG (fl* dtt_2 (fl- (fl* c2 QS tmp1) 
                               (sqr (fl* uBB0 tmp1)))))
                 (fl* dtt_1tmp2c34 r43 uBB0))

              (fl+ (NEG (fl* dtt_2tmp1uBB0 (fl- 2.0 c2)))
                 (fl-   (fl* dtt_1tmp1c34 r43 ))
                 (fl- dtt_1d__))
              (NNEG (fl* dtt_2tmp1c2 uBB1))
              (NNEG (fl* dtt_2tmp1c2 uBB2))
              (NEG  (fl* dtt_2 c2)))
            (RX2/3 R 2)
            (RX2/3 R 3)


            (let ([BIGC (fl+ (NNEG (fl* dtt_2tmp2 uBB0 c2))
                           (fl- (fl* dtt_1tmp2 (- c34 c1345))))])
            (ROTASN V 4 R
              (fl+ (NEG (fl* dtt_2tmp2uBB0 (fl- (fl* 2.0 c2 QS) (fl* c1 U4))))
                   (fl-   (fl* dtt_1 (fl- (fl+ (fl* CCC1  tmp3 (sqr U1))
                                               (fl* CCC2  tmp3 (sqr U2))
                                               (fl* CCC3  tmp3 (sqr U3))
                                               (fl* c1345 tmp2      U4))))))
              (fl+ (NEG (fl* dtt_2 (fl- (fl* c1 U4 tmp1)
                             (fl* c2 (fl+ (fl* (sqr uBB0) tmp2)
                                      (fl* QS  tmp1))))))
                 (fl-   (fl* dtt_1tmp2uBB0 (fl- (fl* r43 c34) c1345))))
              (fl* BIGC uBB1)
              (fl* BIGC uBB2)
              (fl+ (NEG (fl* dtt_2tmp1uBB0 c1))
                 (fl-   (fl* dtt_1tmp1 c1345))
                 (fl-   (fl* dtt_1 d_5)))))
        ))))))))

(let-syntax ([DIAGBODY (lambda (stx)
  (if (equal? (syntax->datum #'DNAME) 'jacld)
    #'(begin 
      (DIAG 1 3)
      (DIAG 1 2)
      (DIAG 1 1)
)
    #'(begin
      (DIAG 2 1)
      (DIAG 2 2)
      (DIAG 2 3))))])
  (DIAGBODY)) 
))))

(define-jac jacld)
(define-jac jacu)

(define (l2norm nx ny nz V sum isize1 jsize1 ksize1)
  (for ([m (in-range 5)]) (f! sum m 0.0))

  (for ([k (kstR nz)])
    (let ([ks1 (fx* k ksize1)])
    (for ([j (jstR ny)])
    (let ([jks1 (fx+ ks1 (fx* j jsize1))])
    (for ([i (istR nx)])
    (let ([ijks1 (fx+ jks1 (fx* i isize1))])
    (for ([m (in-range 5)])
    (f!+ sum m (sqr (fr V (fx+ m ijks1)))))))))))
  
  (let ([ns (exact->inexact (fx* (fx- nx 2) (fx- ny 2) (fx- nz 2)))])
    (for ([m (in-range 5)]) 
      (f! sum m (sqrt (fl/ (fr sum m) ns))))))

(define-syntax-rule (let-syntax-rule ([(N a ...) b ...] ...) B ...)
  (let-syntax ([N (lambda (stx) (syntax-case stx () [(N a ...) (syntax b ...)]))] ...)
    B ...))

(define (compute-surface-integral u nx ny nz isiz1 isiz2 isiz3 isize1 jsize1 ksize1 c2 dxi deta dzeta)
; ist = 2
; jst = 2
; iend = nx -1 
; jend = ny -1 
; ii1 = 2
; ji1 = 2
; ki1 = 3
; ii2 = nx -1
; ji2 = ny -2
; ki2 = nz -1
; ibeg = ii1;
; ifin = ii2;
; jbeg = ji1;
; jfin = ji2;
; ifin1 = ifin - 1;
; jfin1 = jfin - 1;

  (define isize5 (+ ny 2))
  (define phi1 (make-flvector (sqr ( + isiz3 2)) 0.0))
  (define phi2 (make-flvector (sqr ( + isiz3 2)) 0.0))

  (define-syntax-case (FRCIT A)
    (with-syntax-values ([(I IBEG IFIN K KBEG KFIN JBEG JFIN ijkidx1 ijkidx2 DDD)
      (PICK3 #'A ( i 1 (- nx 1) j 1 (- ny 2) 3 (- nz 1) 
                    (+ (* i isize1) (* j jsize1) (* 2 ksize1))
                    (+ (* i isize1) (* j jsize1) (* (- nz 2) ksize1)) (* dxi deta))
                 ( i 1 (- nx 1) k 2 (- nz 1) 2 (- ny 2)
                    (+ (* i isize1) (* 1 jsize1) (* k ksize1))
                    (+ (* i isize1) (* (- ny 3) jsize1) (* k ksize1)) (* dxi dzeta))
                 ( j 1 (- ny 2) k 2 (- nz 1) 2 (- nx 1)
                    (+ (* 1 isize1) (* j jsize1) (* k ksize1))
                    (+ (* (- nx 2) isize1) (* j jsize1) (* k ksize1)) (* deta dzeta)))])
    #'(let ()   
      (for* ([i (in-range (+ isiz2 2))]
             [k (in-range (+ isiz3 2))])
        (let ([ik (+ i (* k isize5))])
          (f! phi1 ik 0.0)
          (f! phi2 ik 0.0)))


      (for* ([K (in-range KBEG KFIN)]
             [I (in-range IBEG IFIN)])
        (let ([ij (+ I (* K isize5))])
          (define-syntax-rule (PHI_CALC PHI_ IJK_)
            (let ([ijk IJK_])
              (f! PHI_  ij (* c2 (- (fr u (+ 4 ijk))
                             (* 0.5 (/ (+ (sqr (fr u (+ 1 ijk)))
                                          (sqr (fr u (+ 2 ijk)))
                                          (sqr (fr u (+ 3 ijk))))
                                       (fr u (+ 0 ijk)))))))))
          (PHI_CALC phi1 ijkidx1)
          (PHI_CALC phi2 ijkidx2)))

      (* DDD (for*/fold ([X 0.0]) 
            ([K (in-range KBEG (sub1 KFIN))]
             [I (in-range IBEG (sub1 IFIN))])
        (define-syntax-rule (stencil V)
          (+ (fr V (+ I (* K isize5)))
             (fr V (+ I 1 (* K isize5)))
             (fr V (+ I (* (+ K 1) isize5)))
             (fr V (+ I 1 (* (+ K 1) isize5)))))
        (+ X (stencil phi1) (stencil phi2)))))))

      (* 0.25 (+ (FRCIT 1) (FRCIT 2) (FRCIT 3))))
   

(define (set-boundary-variables u nx ny nz isize1 jsize1 ksize1)
  (define temp1 (make-flvector 5 0.0))
  (define temp2 (make-flvector 5 0.0))
  (define-syntax-case (INITIT R)
    (let ()
      (define-syntax-rule (with-syntax-values ([(a ...) b] ...) body ...)
        (syntax-case (list b ...) ()
          [( (a ...) ...) (let () body ...)]))

      (define-syntax-case (REPDET* R2 V1 V2 V3 VV1 VV2 VV3)
        #'(case (syntax->datum R2)
          [(1) #'(VV1 V2 V3)]
          [(2) #'(V1 VV2 V3)]
          [(3) #'(V1 V2 VV3)]
          [else (raise (format "invalid REPDEP value ~a" (syntax->datum #'R2)))]))

      (define-syntax-case (PICK3 R2 V1 V2 V3) 
        #'(begin 
          (case (syntax->datum R2) [(1) #'V1] [(2) #'V2] [(3) #'V3]
          [else (raise (format "~a" (syntax->datum R2)))])))

      (with-syntax-values ([(II JJ NII NJJ ISIZE1 JSIZE1)
                              (PICK3 #'R (k j nz ny ksize1 jsize1)
                                         (k i nz nx ksize1 isize1)
                                         (j i ny nx jsize1 isize1))]
                           [(E1 E2 E3) (REPDET* #'R (+ i 1) (+ j 1) (+ k 1) 1 1 1)]
                           [(E4 E5 E6) (REPDET* #'R (+ i 1) (+ j 1) (+ k 1) nx ny nz)]
                           [(KS1) (list (PICK3 #'R (* (- nx 1) isize1) (* (- ny 1) jsize1) (* (- nz 1) ksize1)))])
        #'(for* ([II (in-range NII)]
                 [JJ (in-range NJJ)])
            (let ([IS1 (* II ISIZE1)]
                  [JS1 (* JJ JSIZE1)]
                  [KKS1 KS1])
            (exactKBT E1 E2 E3 temp1 nx ny nz)
            (exactKBT E4 E5 E6 temp2 nx ny nz)
          (for ([m (in-range 5)])
            (let ([idx1 (+ m IS1 JS1)]
                  [idx2 (+ m IS1 JS1 KKS1)]) 
              (f! u idx1 (fr temp1 m))
              (f! u idx2 (fr temp2 m)))))))))

  (INITIT 3)
  (INITIT 2)
  (INITIT 1))

(define (set-initial-values u nx ny nz dnzm1 dnym1 dnxm1 isize1 jsize1 ksize1)
  (define ue_1jk (make-flvector 5 0.0))
  (define ue_i1k (make-flvector 5 0.0))
  (define ue_ij1 (make-flvector 5 0.0))
  (define ue_nx0jk (make-flvector 5 0.0))
  (define ue_iny0k (make-flvector 5 0.0))
  (define ue_ijnz0 (make-flvector 5 0.0))
  (let ([Pface (make-flvector (* 5 3 2) 0.0)])
    (for ([k (in-range 1 (- nz 1))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range 1 (- ny 1))])
          (let ([eta (* j dnym1)])
            (for ([i (in-range 1 (- nx 1))])
              (let ([xi (* i dnxm1)])
                  (exactKBT 1 (+ j 1) (+ k 1) ue_1jk nx ny nz)
                  (exactKBT nx (+ j 1) (+ k 1) ue_nx0jk nx ny nz)
                  (exactKBT (+ i 1) 1 (+ k 1) ue_i1k nx ny nz)
                  (exactKBT (+ i 1) ny (+ k 1) ue_iny0k nx ny nz)
                  (exactKBT (+ i 1) (+ j 1) 1 ue_ij1 nx ny nz)
                  (exactKBT (+ i 1) (+ j 1) nz ue_ijnz0 nx ny nz)
                (for ([m (in-range 5)])
                  (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
                        [pxi   (+ (* (- 1.0 xi) (fr ue_1jk m)) (* xi (fr ue_nx0jk m)))]
                        [peta  (+ (* (- 1.0 eta) (fr ue_i1k m)) (* eta (fr ue_iny0k m)))]
                        [pzeta (+ (* (- 1.0 zeta) (fr ue_ij1 m)) (* zeta (fr ue_ijnz0 m)))])
                    (f! u idx (+ pxi peta pzeta
                                 (- (* pxi peta))
                                 (- (* peta pzeta))
                                 (- (* pzeta pxi))
                                 (* pxi peta pzeta)))))))))))))

(define (sssor cg a b c d u rsd rsdnm tv_ nx ny nz itmax inorm dt omega isize1 jsize1 ksize1 isize4 jsize4 ksize4
          rho_i flux_ qs frct
          r43 c1345 c1 c2 c1c5 c3c4 dssp
          isize2 jsize3 ksize3
          tx1 tx2 tx3
          ty1 ty2 ty3
          tz1 tz2 tz3
          dx1 dy1 dz1
          dx2 dy3 dz2
          dx5 dy5 dz5
          dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
          dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
          dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

  (define delunm (make-flvector 5 0.0))
  (define tmat (make-flvector (fx* 5 5) 0.0))
  (define tolrsd (make-flvector 5 0.00000001))
  (define flux (make-flvector (fx* 5 nx)  0.0))
  (define tv (make-flvector (fx* 5 (fx++ nx) nx)   0.0))


  (CG-n0-only cg
  (for ([j (in-range ny)]
        [i (in-range nx)]
        [n (in-range 5)]
        [m (in-range 5)])
    (let ([idx (fx+ m (fx* n isize4) (fx* i jsize4) (fx* j ksize4))])
      (f! a idx 0.0)
      (f! b idx 0.0)
      (f! c idx 0.0)
      (f! d idx 0.0))))

    (rhs cg flux u qs rho_i rsd frct nz ny nx
      c1 c2 c1c5 c3c4 r43 dssp
      isize1 jsize1 ksize1 isize2 jsize3 ksize3
      tx2 tx3 ty2 ty3 tz2 tz3
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)


  (CG-n0-only cg

    (l2norm nx ny nz rsd rsdnm isize1 jsize1 ksize1)

    (timer-start 1))
  
  (for ([istep (in-range 1 (fx++ itmax))])
    (CG-n0-only cg
      (when (or (zero? (modulo istep 20)) (fx= istep itmax) (fx= istep 1))
        (printf " Time step ~a\n" istep)))

    (for* ([k (p-range cg (in-range 1 (fx-- nz)))]
           [j (jstR ny)]
           [i (istR nx)]
           [m (in-range 5)])
      (f!* rsd (fx+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1)) dt))

    (CG-B cg)

    (for ([k (in-range 1 (fx-- nz))]) 
      (jacld cg a b c d u k rho_i qs nx ny c1 c2 c3c4 r43 c1345
        isize1 jsize1 ksize1 jsize3 ksize3 isize4 jsize4 ksize4
        dt
        tx1 ty1 tz1
        tx2 ty2 tz2
        dx1 dy1 dz1
        dx2 dy3 dz2
        dx5 dy5 dz5
        dx1tx1 dy1ty1 dz1tz1
        dx2tx1 dy2ty1 dz2tz1
        dx3tx1 dy3ty1 dz3tz1
        dx4tx1 dy4ty1 dz4tz1
        dx5tx1 dy5ty1 dz5tz1)
;      (CGpipeline cg
;        (blts cg tmat rsd tv a b c d nx ny nz k omega isize1 jsize1 ksize1 isize4 jsize4 ksize4))
      (CG-n0-only cg
        (blts (CGSingle) tmat rsd tv a b c d nx ny nz k omega isize1 jsize1 ksize1 isize4 jsize4 ksize4))
      ;(CG-n0-only cg (void))
      ;(CG-n0-only cg (pflv rsd "RSD"))
)


;    (CG-n0-only cg (pflv rsd "RSD"))
;    (exit 0)
    (CG-B cg)


    (for ([k (in-range (fx- nz 2) 0 -1)]) 
      (jacu cg a b c d u k rho_i qs nx ny c1 c2 c3c4 r43 c1345
        isize1 jsize1 ksize1 jsize3 ksize3 isize4 jsize4 ksize4
        dt
        tx1 ty1 tz1
        tx2 ty2 tz2
        dx1 dy1 dz1
        dx2 dy3 dz2
        dx5 dy5 dz5
        dx1tx1 dy1ty1 dz1tz1
        dx2tx1 dy2ty1 dz2tz1
        dx3tx1 dy3ty1 dz3tz1
        dx4tx1 dy4ty1 dz4tz1
        dx5tx1 dy5ty1 dz5tz1)
;      (CGpipeline cg
;        (buts cg tmat rsd tv c b a d nx ny nz k omega isize1 jsize1 ksize1 isize4 jsize4 ksize4))
      (CG-n0-only cg
        (buts (CGSingle) tmat rsd tv c b a d nx ny nz k omega isize1 jsize1 ksize1 isize4 jsize4 ksize4))
      )

    (CG-n0-only cg
    (let ([tmp (fl/ 1.0 (* omega (- 2.0 omega)))]) 
      (for* ([k (kstR nz)]
             [j (jstR ny)]
             [i (istR nx)]
             [m (in-range 5)])
        (let ([idx (fx+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1))])
          (f!+ u idx (fl* tmp (fr rsd idx)))))))

    (CG-n0-only cg
      (when (zero? (modulo istep inorm))
        (l2norm nx ny nz rsd delunm isize1 jsize1 ksize1)))

    (rhs cg flux u qs rho_i rsd frct nz ny nx
      c1 c2 c1c5 c3c4 r43 dssp
      isize1 jsize1 ksize1 isize2 jsize3 ksize3
      tx2 tx3 ty2 ty3 tz2 tz3
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

    (CG-n0-only cg
      (when (or (zero? (modulo istep 20)) (fx= istep itmax))
        (l2norm nx ny nz rsd rsdnm isize1 jsize1 ksize1))

      (when (and ((fr rsdnm 0) . < . (fr tolrsd 0))
               ((fr rsdnm 1) . < . (fr tolrsd 1))
               ((fr rsdnm 2) . < . (fr tolrsd 2))
               ((fr rsdnm 3) . < . (fr tolrsd 3))
               ((fr rsdnm 4) . < . (fr tolrsd 4)))
        (timer-stop 1))))
  (timer-stop 1))


(define (verify class xcr xce xci dt)
  (define-values (xcrdif xcedif) (values (make-flvector 5 0.0) (make-flvector 5 0.0)))
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
    (printf "~a. ~a ~a ~a\n" m (fr xcr m) (fr xcrref m) (fr xcrdif m)))

  (if (not (equal? class #\U))
    (printf " Comparison of RMS-norms of solution error\n")
    (printf " RMS-norms of solution error"))

  (for ([m (in-range (flvector-length xce))])
    (printf "~a. ~a ~a ~a\n" m (fr xce m) (fr xceref m) (fr xcedif m)))

  (if (not (equal? class #\U))
    (printf " Comparison of surface integral\n")
    (printf " Surface integral"))

  (printf "~a ~a ~a\n" xci xciref xcidif)
  
  (and ((abs (- dt dtref)) . <= . epsilon)
       (<epsilon-vmap xcrdif epsilon)
       (<epsilon-vmap xcedif epsilon)
       (xcidif . < . epsilon)))

#|
 vim set makeprg=/home/tewk/Tools/bin/rktl\ \-tm\ %\ SERIAL\ CLASS=S
 vim set errorformat=%f:%l:%m
|#
