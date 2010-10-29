#lang racket/base

(provide main)
  
(require "../bm-args.rkt") 
(require "../bm-results.rkt") 
(require "../rand-generator.rkt")
(require "../timer.rkt")
(require "../parallel-utils.rkt")
(require "../debug.rkt")
(require (rename-in "../macros.rkt" [f!+ flvs!+] [f!- flvs!-] [f!* flvs!*] [f!/ flvs!/]))
(require (for-syntax "../macros.rkt"))
(require racket/match)
(require racket/math)
(require (for-syntax scheme/base))

#|
(require (only-in scheme/flonum make-flvector 
                                make-shared-flvector 
                                shared-flvector
                                flvector-length
                                flvector
                                flmax
                                flvector-set!
                                flvector-ref))
(define vr vector-ref)
(define vs! vector-set!)
(define flvs! flvector-set!)
(define flvr flvector-ref)
|#


(require (only-in scheme/flonum make-flvector 
                                make-shared-flvector 
                                shared-flvector
                                flvector-length
                                flvector
                                ))

(require (rename-in scheme/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-flvector-ref flvr] 
                    [unsafe-flvector-set! flvs!]
))

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values  12 0.015 100)]
    [(#\W) (values  36 0.0015 400)]
    [(#\A) (values  64 0.0015 400)]
    [(#\B) (values  102 0.001 400)]
    [(#\C) (values  162 0.00067 400)]
    [else (error "Unknown class")]))

(define (get-verify-values class)
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
        0.00067)]
    [else (values
        (make-flvector 5 1.0)
        (make-flvector 5 1.0)
        0.00001)]))


 (define ce (shared-flvector 
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

(define (run-benchmark args) 
  (define maxlevel 11)
  (let ([bmname "SP"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(problem_size dt_default niter_default) (get-class-size CLASS)])
    (let* (
          [niter niter_default]
          [dt dt_default]
          [IMAX problem_size]
          [JMAX problem_size]
          [KMAX problem_size]
          [nx problem_size]
          [ny problem_size]
          [nz problem_size]
          [nx2 (- problem_size 2)]
          [ny2 (- problem_size 2)]
          [nz2 (- problem_size 2)]
          [dnxm1 (/ 1.0 (- problem_size 1))]
          [dnym1 (/ 1.0 (- problem_size 1))]
          [dnzm1 (/ 1.0 (- problem_size 1))]
          [isize1 5]
          [jsize1 (* 5 (add1 IMAX))]
          [ksize1 (* 5 (add1 IMAX) (add1 JMAX))]
          [jsize2 (add1 IMAX)]
          [ksize2 (* (add1 IMAX) (add1 JMAX))]
          [jsize3 problem_size]
          [jsize4 5]
          [s1 (* ksize1 KMAX)]
          [s2 (* ksize2 KMAX)]
          [s3 (* 5 (add1 problem_size))]
          [u       (make-shared-flvector s1 0.0)]
          [rhs     (make-shared-flvector s1 0.0)]
          [forcing (make-shared-flvector s1 0.0)]
          [us      (make-shared-flvector s2 0.0)]
          [vs      (make-shared-flvector s2 0.0)]
          [ws      (make-shared-flvector s2 0.0)]
          [qs      (make-shared-flvector s2 0.0)]
          [rho_i   (make-shared-flvector s2 0.0)]
          [speed   (make-shared-flvector s2 0.0)]
          [square  (make-shared-flvector s2 0.0)]
          [lhs     (make-shared-flvector s3 0.0)]
          [lhsp    (make-shared-flvector s3 0.0)]
          [lhsm    (make-shared-flvector s3 0.0)]
          [cv      (make-shared-flvector problem_size 0.0)]
          [rhon    (make-shared-flvector problem_size 0.0)]
          [rhos    (make-shared-flvector problem_size 0.0)]
          [rhoq    (make-shared-flvector problem_size 0.0)]
          [cuf     (make-shared-flvector problem_size 0.0)]
          [q       (make-shared-flvector problem_size 0.0)]
          [ue      (make-shared-flvector (* 5 problem_size) 0.0)]
          [buf     (make-shared-flvector (* 5 problem_size) 0.0)]
          [tx1     (/ 1.0 (sqr dnxm1))]
          [ty1     (/ 1.0 (sqr dnym1))]
          [tz1     (/ 1.0 (sqr dnzm1))]
          [tx2     (/ 1.0 (* 2.0 dnxm1))]
          [ty2     (/ 1.0 (* 2.0 dnym1))]
          [tz2     (/ 1.0 (* 2.0 dnzm1))]
          [tx3     (/ 1.0 dnxm1)]
          [ty3     (/ 1.0 dnym1)]
          [tz3     (/ 1.0 dnzm1)]
          [bt (sqrt 0.5)]
          [c1 1.4]
          [c2 0.4]
          [c1c2 (* 1.4 0.4)]
          [c1c5 (* 1.4 1.4)]
          [c3c4 (* 0.1 1.0)]
          [c1345 (* 1.4 1.4 0.1 1.0)]
          [c3c4tx3 (* c3c4 tx3)]
          [c3c4ty3 (* c3c4 ty3)]
          [c3c4tz3 (* c3c4 tz3)]
          [con43   (/ 4.0 3.0)]
          [conz1   (- 1.0 c1c5)]
          [con16   (/ 1.0 6.0)]
          [dssp    0.25]
          [xxcon1  (* c3c4tx3 con43 tx3)]
          [xxcon2  (* c3c4tx3 tx3)]
          [xxcon3  (* c3c4tx3 conz1 tx3)]
          [xxcon4  (* c3c4tx3 con16 tx3)]
          [xxcon5  (* c3c4tx3 c1c5 tx3)]
          [yycon1  (* c3c4ty3 con43 ty3)]
          [yycon2  (* c3c4ty3 ty3)]
          [yycon3  (* c3c4ty3 conz1 ty3)]
          [yycon4  (* c3c4ty3 con16 ty3)]
          [yycon5  (* c3c4ty3 c1c5 ty3)]
          [zzcon1  (* c3c4tz3 con43 tz3)]
          [zzcon2  (* c3c4tz3 tz3)]
          [zzcon3  (* c3c4tz3 conz1 tz3)]
          [zzcon4  (* c3c4tz3 con16 tz3)]
          [zzcon5  (* c3c4tz3 c1c5 tz3)]
          [dx1tx1  (* 0.75 tx1)]
          [dx2tx1  (* 0.75 tx1)]
          [dx3tx1  (* 0.75 tx1)]
          [dx4tx1  (* 0.75 tx1)]
          [dx5tx1  (* 0.75 tx1)]
          [dy1ty1  (* 0.75 ty1)]
          [dy2ty1  (* 0.75 ty1)]
          [dy3ty1  (* 0.75 ty1)]
          [dy4ty1  (* 0.75 ty1)]
          [dy5ty1  (* 0.75 ty1)]
          [dz1tz1  (* 1.0 tz1)]
          [dz2tz1  (* 1.0 tz1)]
          [dz3tz1  (* 1.0 tz1)]
          [dz4tz1  (* 1.0 tz1)]
          [dz5tz1  (* 1.0 tz1)]
          [dx1 0.75]
          [dx2 0.75]
          [dx3 0.75]
          [dx4 0.75]
          [dx5 0.75]
          [dy1 0.75]
          [dy2 0.75]
          [dy3 0.75]
          [dy4 0.75]
          [dy5 0.75]
          [dz1 1.00]
          [dz2 1.00]
          [dz3 1.00]
          [dz4 1.00]
          [dz5 1.00]
          [dttx1 (* dt tx1)]
          [dttx2 (* dt tx2)]
          [dtty1 (* dt ty1)]
          [dtty2 (* dt ty2)]
          [dttz1 (* dt tz1)]
          [dttz2 (* dt tz2)]
          [c2dttx1 (* 2.0 dttx1)]
          [c2dtty1 (* 2.0 dtty1)]
          [c2dttz1 (* 2.0 dttz1)]
          [dxmax (flmax dx3 dx4)]
          [dymax (flmax dy2 dy4)]
          [dzmax (flmax dz2 dz3)]
          [dtdssp (* dt dssp)]
          [comz1 dtdssp]
          [comz4 (* 4.0 dtdssp)]
          [comz5 (* 5.0 dtdssp)]
          [comz6 (* 6.0 dtdssp)]
)
  (define (compute_rhs_thunk)
    (compute_rhs (CGSingle) isize1 jsize1 ksize1 jsize2 ksize2 u us vs ws rho_i square qs speed
        c1c2 rhs forcing nx2 ny2 nz2 c1 c2 dssp
        tx2 ty2 tz2 con43 dt
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        xxcon2 xxcon3 xxcon4 xxcon5
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        yycon2 yycon3 yycon4 yycon5
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
        zzcon2 zzcon3 zzcon4 zzcon5))

   (print-banner bmname args) 

;;;//---------------------------------------------------------------------
;;;//      Read input file (if it exists), else take
;;;//      defaults from parameters
;;;//---------------------------------------------------------------------
;      (get-input-pars)
      (printf "No input file inputsp.data,Using compiled defaults\n")
      (printf "Size: ~a X ~a X ~a\n" nx ny nz)
      (printf "Iterations: ~a dt: ~a\n" niter dt)
      (initialize u nx ny nz isize1 jsize1 ksize1 dnxm1 dnym1 dnzm1)
      (exact_rhs nx2 ny2 nz2 isize1 jsize1 ksize1 jsize3 forcing dnxm1 dnym1 dnzm1 ue buf cuf q
        rhs u c1 c2 0.25
        tx2 ty2 tz2
        xxcon1 xxcon2 xxcon3 xxcon4 xxcon5
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        yycon1 yycon2 yycon3 yycon4 yycon5
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        zzcon1 zzcon2 zzcon3 zzcon4 zzcon5
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

      (CGspawn (if serial 0 num-threads) sp-body 
        u us vs ws cv rhon rhoq rhos rho_i square qs speed rhs forcing lhs lhsp lhsm
        nx ny nz 
        nx2 ny2 nz2 
        c1 c2 dssp c1c2 c1c5 con43 c3c4
        tx2 ty2 tz2 dt bt niter
        dx2 dx5 dx1 dttx1 dttx2 c2dttx1
        dy3 dy5 dy1 dtty1 dtty2 c2dtty1
        dz4 dz5 dz1 dttz1 dttz2 c2dttz1
        dxmax dymax dzmax 
        dnxm1 dnym1 dnzm1
        comz1 comz4 comz5 comz6
        isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        xxcon2 xxcon3 xxcon4 xxcon5
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        yycon2 yycon3 yycon4 yycon5
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
        zzcon2 zzcon3 zzcon4 zzcon5)

      (let* ([verified (verify CLASS niter dt compute_rhs_thunk
        nx2 ny2 nz2 isize1 jsize1 ksize1 u rhs dnzm1 dnym1 dnxm1)])
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

(define (sp-body cg u us vs ws cv_ rhon rhoq rhos rho_i square qs speed rhs forcing lhs_ lhsp_ lhsm_
      nx ny nz 
      nx2 ny2 nz2 c1 c2 dssp c1c2 c1c5 con43 c3c4
      tx2 ty2 tz2 dt bt niter
      dx2 dx5 dx1 dttx1 dttx2 c2dttx1
      dy3 dy5 dy1 dtty1 dtty2 c2dtty1
      dz4 dz5 dz1 dttz1 dttz2 c2dttz1
      dxmax dymax dzmax 
      dnxm1 dnym1 dnzm1
      comz1 comz4 comz5 comz6
      isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
      xxcon2 xxcon3 xxcon4 xxcon5
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
      yycon2 yycon3 yycon4 yycon5
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
      zzcon2 zzcon3 zzcon4 zzcon5)

  (define s3 (* 5 (add1 nx)))
  (define lhs  (make-flvector s3 0.0))
  (define lhsp (make-flvector s3 0.0))
  (define lhsm (make-flvector s3 0.0))
  (define rho (make-flvector nx 0.0))
  (define cv (make-flvector nx 0.0))

  (define (adi step)
    (compute_rhs cg isize1 jsize1 ksize1 jsize2 ksize2 u us vs ws rho_i square qs speed
      c1c2 rhs forcing nx2 ny2 nz2 c1 c2 dssp
      tx2 ty2 tz2 con43 dt
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
      xxcon2 xxcon3 xxcon4 xxcon5
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
      yycon2 yycon3 yycon4 yycon5
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
      zzcon2 zzcon3 zzcon4 zzcon5)

    (CG-B cg)

    (txinvr cg nz2 ny2 nx2 isize1 jsize1 ksize1 jsize2 ksize2 rho_i us vs ws qs speed rhs c2 bt)
    
    (CG-B cg)

    (x_solve cg nz2 ny2 nx2 us rho
      dx2 dx5 dx1 dttx1 dttx2 c2dttx1
      isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
      rho_i cv con43 c1c5 c3c4 dxmax bt rhs lhs speed lhsp lhsm
      comz1 comz4 comz5 comz6)

    (CG-B cg)

    (y_solve cg nz2 ny2 nx2 vs rho
      dy3 dy5 dy1 dtty1 dtty2 c2dtty1
      isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
      rho_i cv con43 c1c5 c3c4 dymax bt rhs lhs speed lhsp lhsm
      comz1 comz4 comz5 comz6)

    (CG-B cg)

    (z_solve cg nz2 ny2 nx2 ws rho
      dz4 dz5 dz1 dttz1 dttz2 c2dttz1
      isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
      rho_i cv con43 c1c5 c3c4 dzmax bt rhs lhs speed lhsp lhsm
      comz1 comz4 comz5 comz6
      u us vs qs)

    (CG-B cg)

    (add cg nz2 ny2 nx2 isize1 jsize1 ksize1 u rhs))

;;;//---------------------------------------------------------------------
;;;//      do one time step to touch all code, and reinitialize
;;;//---------------------------------------------------------------------
  (adi "I")

  (CG-n0-only cg
    (initialize u nx ny nz isize1 jsize1 ksize1 dnxm1 dnym1 dnzm1))

  (CG-n0-only cg
    (timer-start 1))

  (for ([step (in-range 1 (add1 niter))])
    (CG-n0-only cg
      (when (or (zero? (modulo step 20)) (= step 1) (= step niter))
        (printf "Time step ~a\n" step)))

    (adi step))
  


  (CG-n0-only cg
    (timer-stop 1)))



(define (get-mflops total-time niter nx ny nz)
  (if (not (= total-time 0.0))
    (let* ([n3 (* nx ny nz)]
           [t  (/ (+ nx ny nz) 3.0)])
      (/ (* (+ (* 881.174 n3)
               (* -4683.91 (* t t))
               (* 11484.5 t)
               -19272.4)
            niter)
          (* total-time 1000000.0)))
      0.0))

(define (exact_solution xi eta zeta dtemp offset)
  (for ([m (in-range 5)])
    (flvs! dtemp (+ m offset) 
      (+ (flvr ce m)
         (* xi (+ (flvr ce (+ m 5))
                  (* xi (+ (flvr ce (+ m (* 4 5)))
                           (* xi (+ (flvr ce (+ m (* 7 5)))
                                    (* xi (flvr ce (+ m (* 10 5))))))))))
         (* eta (+ (flvr ce (+ m (* 2 5)))
                  (* eta (+ (flvr ce (+ m (* 5 5)))
                           (* eta (+ (flvr ce (+ m (* 8 5)))
                                    (* eta (flvr ce (+ m (* 11 5))))))))))
         (* zeta (+ (flvr ce (+ m (* 3 5)))
                  (* zeta (+ (flvr ce (+ m (* 6 5)))
                           (* zeta (+ (flvr ce (+ m (* 9 5)))
                                    (* zeta (flvr ce (+ m (* 12 5))))))))))))))
  
(define (initialize u nx ny nz isize1 jsize1 ksize1 dnxm1 dnym1 dnzm1)
  (for* ([k (in-range nz)]
         [j (in-range ny)]
         [i (in-range nx)])
    (let ([idx (+ (* i isize1) (* j jsize1) (* k ksize1))])
      (flvs! u (+ 0 idx) 1.0)
      (flvs! u (+ 1 idx) 0.0)
      (flvs! u (+ 2 idx) 0.0)
      (flvs! u (+ 3 idx) 0.0)
      (flvs! u (+ 4 idx) 1.0)))

  (let ([Pface (make-flvector (* 5 3 2) 0.0)])
    (for ([k (in-range nz)])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range ny)])
          (let ([eta (* j dnym1)])
            (for ([i (in-range nx)])
              (let ([xi (* i dnxm1)])
                (for ([ix (in-range 2)])
                  (exact_solution ix eta zeta Pface (+ 0 (* 0 5) (* ix 15))))
                (for ([ix (in-range 2)])
                  (exact_solution xi ix zeta Pface (+ 0 (* 1 5) (* ix 15))))
                (for ([ix (in-range 2)])
                  (exact_solution xi eta ix Pface (+ 0 (* 2 5) (* ix 15))))

                (for ([m (in-range 5)])
                  (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
                        [pxi   (+ (* xi   (flvr Pface (+ m (* 0 5) (* 1 15))))
                                 (* (- 1.0   xi) (flvr Pface (+ m (* 0 5) (* 0 15)))))]
                        [peta  (+ (* eta  (flvr Pface (+ m (* 1 5) (* 1 15))))
                                 (* (- 1.0  eta) (flvr Pface (+ m (* 1 5) (* 0 15)))))]
                        [pzeta (+ (* zeta (flvr Pface (+ m (* 2 5) (* 1 15))))
                                 (* (- 1.0 zeta) (flvr Pface (+ m (* 2 5) (* 0 15)))))])
                    (flvs! u idx  (+ (- (+ pxi peta pzeta) (* pxi peta) (* pxi pzeta) (* peta pzeta))
                                     (* pxi peta pzeta))))))))))))
  (let ([temp (make-flvector 5 0.0)]
        [temp2 (make-flvector 5 0.0)]
        [i2 (sub1 nx)]
        [j2 (sub1 ny)]
        [k2 (sub1 nz)])
    (for ([k (in-range nz)])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range ny)])
          (let ([eta (* j dnym1)])
            (exact_solution 0.0 eta zeta temp 0)
            (exact_solution 1.0 eta zeta temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* j jsize1) (* k ksize1))]
                     [idx2 (+ (* i2 isize1) idx)])
                (flvs! u idx (flvr temp m))
                (flvs! u idx2 (flvr temp2 m))))))
            
        (for ([i (in-range nx)])
          (let ([xi (* i dnxm1)])
            (exact_solution xi 0.0 zeta temp 0)
            (exact_solution xi 1.0 zeta temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* i isize1) (* k ksize1))]
                     [idx2 (+ (* j2 jsize1) idx)])
                (flvs! u idx (flvr temp m))
                (flvs! u idx2 (flvr temp2 m))))))))
  
    (for ([j (in-range ny)])
      (let ([eta (* j dnym1)])
        (for ([i (in-range nx)])
          (let ([xi (* i dnxm1)])
            (exact_solution xi eta 0.0 temp 0)
            (exact_solution xi eta 1.0 temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* i isize1) (* j jsize1))]
                     [idx2 (+ (* k2 ksize1) idx)])
                (flvs! u idx (flvr temp m))
                (flvs! u idx2 (flvr temp2 m))))))))))

(define (lhsinit size jsize4 lhs lhsp lhsm)
  (for* ([i (in-range 0 (add1 size) size)]
         [n (in-range 5)])
    (let ([nij (+ n (* i jsize4))])
      (flvs! lhs nij 0.0)
      (flvs! lhsp nij 0.0)
      (flvs! lhsm nij 0.0)))
  (for ([i (in-range 0 (add1 size) size)])
    (let ([nij (+ 2 (* i jsize4))])
      (flvs! lhs nij 1.0)
      (flvs! lhsp nij 1.0)
      (flvs! lhsm nij 1.0))))

(define (add cg nz2 ny2 nx2 isize1 jsize1 ksize1 u rhs)
  (CGfor cg ([k (in-range 1 (add1 nz2))])
    (let ([ks1 (* k ksize1)])
      (for ([j (in-range 1 (add1 ny2))])
        (let ([js1 (+ (* j jsize1) ks1)])
          (for ([i (in-range 1 (add1 nx2))])
            (let ([is1 (+ (* i isize1) js1)])
              (for ([m (in-range 5)])
                (let ([ms1 (+ m is1)])
                (flvs!+ u ms1 (flvr rhs ms1)))))))))))

(define (error-norm rms nx2 ny2 nz2 isize1 jsize1 ksize1 u dnzm1 dnym1 dnxm1)
  (for ([m (in-range 5)]) (flvs! rms m 0.0))

  (let ([u-exact (make-flvector 5 0.0)])
    (for ([k (in-range (+ nz2 2))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range (+ ny2 2))])
          (let ([eta (* j dnym1)])
            (for ([i (in-range (+ nx2 2))])
              (let ([xi (* i dnxm1)])
                (exact_solution xi eta zeta u-exact 0)
                (for ([m (in-range 5)])
                  (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
                         [add (- (flvr u idx) (flvr u-exact m))])
                    (flvs!+ rms m (* add add)))))))))))

  
  (for ([m (in-range 5)])
    (flvs! rms m (sqrt (/ (flvr rms m) nx2 ny2 nz2)))))

(define (rhs-norm rms nz2 ny2 nx2 isize1 jsize1 ksize1 rhs)
  (for ([m (in-range 5)]) (flvs! rms m 0.0))

  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let* ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))]
           [add (flvr rhs idx)])
      (flvs!+ rms m (* add add))))
  
  (for ([m (in-range 5)])
    (flvs! rms m (sqrt (/ (flvr rms m) nx2 ny2 nz2)))))

(define-syntax-rule (fourth-order-dissipation ii nII2 V V2 i j k m midx idx
  MIDX MIDX+ MIDX+2 MIDX- MIDX-2 
  DIDX IDX dssp)

  (begin
  (let* ([ii 1]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx+  MIDX+]
             [midx+2 MIDX+2]
             [didx   DIDX])
        (flvs!+ V didx (- (fl* dssp (+ (fl*  5.0 (flvr V2 midx))
                                           (fl* -4.0 (flvr V2 midx+))
                                                     (flvr V2 midx+2))))))))
  (let* ([ii 2]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx+  MIDX+]
             [midx+2 MIDX+2]
             [midx-  MIDX-]
             [didx   DIDX])
        (flvs!+ V didx (- (fl* dssp (+ (fl* -4.0 (flvr V2 midx-)) 
                                           (fl*  6.0 (flvr V2 midx))
                                           (fl* -4.0 (flvr V2 midx+))
                                                     (flvr V2 midx+2))))))))
  (for ([ii (in-range 3 (sub1 nII2))])
    (let ([idx IDX])
      (for ([m (in-range 5)])
        (let* ([midx   MIDX]
               [midx+  MIDX+]
               [midx+2 MIDX+2]
               [midx-  MIDX-]
               [midx-2 MIDX-2]
               [didx   DIDX])
          (flvs!+ V didx (- (fl* dssp (+           (flvr V2 midx-2)
                                             (fl* -4.0 (flvr V2 midx-))
                                             (fl*  6.0 (flvr V2 midx))
                                             (fl* -4.0 (flvr V2 midx+))
                                                       (flvr V2 midx+2)))))))))
  (let* ([ii (sub1 nII2)]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx+  MIDX+]
             [midx-  MIDX-]
             [midx-2 MIDX-2]
             [didx   DIDX])
        (flvs!+ V didx (- (fl* dssp (+           (flvr V2 midx-2)
                                           (fl* -4.0 (flvr V2 midx-))
                                           (fl*  6.0 (flvr V2 midx))
                                           (fl* -4.0 (flvr V2 midx+)))))))))
  (let* ([ii nII2]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx-  MIDX-]
             [midx-2 MIDX-2]
             [didx   DIDX])
        (flvs!+ V didx (- (fl* dssp (+           (flvr V2 midx-2)
                                           (fl* -4.0 (flvr V2 midx-))
                                           (fl*  5.0 (flvr V2 midx)))))))))))

(define (exact_rhs nx2 ny2 nz2 isize1 jsize1 ksize1 jsize3 forcing dnxm1 dnym1 dnzm1 ue buf cuf q
  rhs u c1 c2 dssp
  tx2 ty2 tz2
  xxcon1 xxcon2 xxcon3 xxcon4 xxcon5
  dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
  yycon1 yycon2 yycon3 yycon4 yycon5
  dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
  zzcon1 zzcon2 zzcon3 zzcon4 zzcon5
  dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
      (flvs!+ forcing idx 0.0)))

(define-syntax-rule (body2 ii jj kk NII2 NJJ2 NKK2 t_2 f1jsize3 jsize3 ue buf cuf q c1 c2 i j k hh dtemp
  t_2m1 t_2m2 t_2m3 
  _con_0 _con_1 _con_2 __con3 __con4 __con5
  d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1)
  (begin
    (for ([kk (in-range 1 (add1 NKK2))])
        (for ([jj (in-range 1 (add1 NJJ2))])
          (for ([ii (in-range 0 (+ NII2 2))])
              (let ([xi (* i dnxm1)]
                    [eta (* j dnym1)]
                    [zeta (* k dnzm1)])
                (exact_solution xi eta zeta dtemp 0)
                (let ([dtpp (/ 1.0 (flvr dtemp 0))])
                  (for ([m (in-range 5)])
                    (let ([idx (+ ii (* m jsize3))]
                          [dtv (flvr dtemp m)])
                    (flvs! ue idx dtv)
                    (flvs! buf idx (* dtpp dtv)))))
                (let* ([ij (+ ii jsize3)]
                       [i2j (+ ii (* 2 jsize3))]
                       [i3j (+ ii (* 3 jsize3))]
                       [bufij (flvr buf ij)]
                       [bufi2j (flvr buf i2j)]
                       [bufi3j (flvr buf i3j)]
                       [ueij (flvr ue ij)]
                       [uei2j (flvr ue i2j)]
                       [uei3j (flvr ue i3j)]
                       [bufij2 (sqr bufij)]
                       [bufi2j2 (sqr bufi2j)]
                       [bufi3j2 (sqr bufi3j)])
                (flvs! cuf ii (sqr (flvr buf (+ ii (* hh jsize3)))))
                (flvs! buf ii (+ bufij2 bufi2j2 bufi3j2))
                (flvs! q ii (* 0.5 (+ (* bufij ueij) (* bufi2j uei2j) (* bufi3j uei3j)))))))

          (for ([ii (in-range 1 (add1 NII2))])
            (let* ([ip1 (add1 ii)]
                   [im1 (sub1 ii)]
                   [didx (+ (* i isize1) (* j jsize1) (* k ksize1))]
                   [idx2 (+ ii f1jsize3)]
                   [idx2+ (+ ip1 f1jsize3)]
                   [idx2- (+ im1 f1jsize3)]
                   [A4J (+ ii (* 4 jsize3))]
                   [A4J+ (+ A4J 1)]
                   [A4J- (- A4J 1)])

              (define-syntax-rule (citA C C1 C2 C3) (* C (+ (- C1 (* 2.0 C2)) C3)))
              (define-syntax-rule (citS3 C V I1 I2 I3) (citA C (flvr V I1) (flvr V I2) (flvr V I3)))
              (define-syntax-rule (citS C V) (citA C (flvr V ip1) (flvr V ii) (flvr V im1)))

              (define-syntax-rule (t_2it l r) (- (* t_2 (-  l r))))
              (define-syntax-rule (t_2it3 l r o) (- (* t_2 (+ (-  l r) o))))
              (define-syntax-rule (t_2itlr UI ZSI) (* (flvr ue UI) (flvr buf ZSI)))
              (define-syntax-rule (t_2ito) (- (* c2 (- (flvr ue A4J+) (flvr q ip1)))
                                              (* c2 (- (flvr ue A4J-) (flvr q im1)))))

              (define-syntax-rule (mid d__t_1 __con2X A t_2itother)
                   (let* ([AJ (+ ii (* A jsize3))]
                          [AJ+ (+ AJ 1)]
                          [AJ- (- AJ 1)])
                (flvs!+ forcing (+ didx A)
                  (citS3 d__t_1 ue AJ+ AJ AJ-)
                  (citS3 __con2X buf AJ+ AJ AJ-)
                  (t_2it3 (t_2itlr AJ+ idx2+)
                         (t_2itlr AJ- idx2-)
                         (t_2itother (t_2ito))))))

              (flvs!+ forcing didx 
                   (t_2it (flvr ue idx2+)
                          (flvr ue idx2-))
                   (citS d_1t_1 ue))

              (mid d_2t_1 _con_0 1 t_2m1)
              (mid d_3t_1 _con_1 2 t_2m2)
              (mid d_4t_1 _con_2 3 t_2m3)

              (define-syntax-rule (t4clause I1 I2 I3)
                (* (flvr buf I1) (- (* c1 (flvr ue I2)) (* c2 (flvr q I3)))))
              (flvs!+ forcing (+ didx 4)
                   (t_2it (t4clause idx2+ A4J+ ip1) (t4clause idx2- A4J- im1))
                   (* 0.5 (citS __con3 buf))
                   (citS __con4 cuf)
                   (citS3 __con5 buf A4J+ A4J A4J-)
                   (citS3 d_5t_1 ue A4J+ A4J A4J-))))

;//---------------------------------------------------------------------
;//            Fourth-order dissipation
;//---------------------------------------------------------------------
    (fourth-order-dissipation ii NII2 forcing ue i j k m midx idx
              (+ ii (* m jsize3)) (+ midx 1) (+ midx 2) (- midx 1) (- midx 2) 
              (+ m idx) (+ (* i isize1) (* j jsize1) (* k ksize1)) dssp)
))))

  (define-syntax-rule (KZERO a ...) 0)
  (define-syntax-rule (KIDENT a ...) (begin a ...))

  (let ([dtemp (make-flvector 5 0.0)])
   (body2 i j k nx2 ny2 nz2 tx2 jsize3 jsize3 ue buf cuf q c1 c2 i j k 1 dtemp
      KIDENT KZERO KZERO
      xxcon1 xxcon2 xxcon2 xxcon3 xxcon4 xxcon5
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)

    (body2 j i k ny2 nx2 nz2 ty2 (* 2 jsize3) jsize3 ue buf cuf q c1 c2 i j k 2 dtemp 
      KZERO KIDENT KZERO
      yycon2 yycon1 yycon2 yycon3 yycon4 yycon5
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)

    (body2 k i j nz2 nx2 ny2 tz2 (* 3 jsize3) jsize3 ue buf cuf q c1 c2 i j k 3 dtemp
      KZERO KZERO KIDENT 
      zzcon2 zzcon2 zzcon1 zzcon3 zzcon4 zzcon5
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1))

    (for* ([k (in-range 1 (add1 nz2))]
           [j (in-range 1 (add1 ny2))]
           [i (in-range 1 (add1 nx2))]
           [m (in-range 5)])
      (flvs!* forcing (+ m (* i isize1) (* j jsize1) (* k ksize1)) -1.0))
)

(define (ninvr cg nz2 ny2 nx2 isize1 jsize1 ksize1 rhs bt)
  (CGfor cg ([k (in-range 1 (add1 nz2))])
    (for* ([j (in-range 1 (add1 ny2))]
           [i (in-range 1 (add1 nx2))])
    (let* ([idx (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
           [r1 (flvr rhs (fx+ 0 idx))]
           [r2 (flvr rhs (fx+ 1 idx))]
           [r3 (flvr rhs (fx+ 2 idx))]
           [r4 (flvr rhs (fx+ 3 idx))]
           [r5 (flvr rhs (fx+ 4 idx))]
           [t1 (fl* bt r3)]
           [t2 (fl* 0.5 (fl+ r4 r5))])
      (flvs! rhs (fx+ 0 idx) (fl- r2))
      (flvs! rhs (fx+ 1 idx) r1)
      (flvs! rhs (fx+ 2 idx) (fl* bt (fl- r4 r5)))
      (flvs! rhs (fx+ 3 idx) (fl- t2 t1))
      (flvs! rhs (fx+ 4 idx) (fl+ t1 t2))))))

(define (pinvr cg nz2 ny2 nx2 isize1 jsize1 ksize1 rhs bt)
  (CGfor cg ([k (in-range 1 (add1 nz2))])
    (for* ([j (in-range 1 (add1 ny2))]
           [i (in-range 1 (add1 nx2))])
    (let* ([idx (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
           [r1 (flvr rhs (fx+ 0 idx))]
           [r2 (flvr rhs (fx+ 1 idx))]
           [r3 (flvr rhs (fx+ 2 idx))]
           [r4 (flvr rhs (fx+ 3 idx))]
           [r5 (flvr rhs (fx+ 4 idx))]
           [t1 (fl* bt r1)]
           [t2 (fl* 0.5 (fl+ r4 r5))])
      (flvs! rhs (fx+ 0 idx) (fl* bt (fl- r4 r5)))
      (flvs! rhs (fx+ 1 idx) (fl- r3))
      (flvs! rhs (fx+ 2 idx) r2)
      (flvs! rhs (fx+ 3 idx) (fl- t2 t1))
      (flvs! rhs (fx+ 4 idx) (fl+ t1 t2))))))

(define (compute_rhs cg isize1 jsize1 ksize1 jsize2 ksize2 u us vs ws rho_i square qs speed
c1c2 rhs forcing nx2 ny2 nz2 c1 c2 dssp
    tx2 ty2 tz2 con43 dt
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
    xxcon2 xxcon3 xxcon4 xxcon5
    dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
    yycon2 yycon3 yycon4 yycon5
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
    zzcon2 zzcon3 zzcon4 zzcon5
)
  (CGfor cg ([k (in-range 0 (fx+ nz2 2))])
    (for* ([j (in-range (fx+ ny2 2))]
           [i (in-range (fx+ nx2 2))])
      (let* ([idx (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
             [idx2 (fx+ i (fx* j jsize2) (fx* k ksize2))]
             [rho_inv (/ 1.0 (flvr u idx))]
             [u1 (flvr u (fx+ idx 1))]
             [u2 (flvr u (fx+ idx 2))]
             [u3 (flvr u (fx+ idx 3))]
             [u4 (flvr u (fx+ idx 4))]
             [sq (fl* 0.5 (fl+ (sqr u1) (sqr u2) (sqr u3)) rho_inv)])
      (flvs! rho_i idx2 rho_inv)
      (flvs! us idx2 (fl* rho_inv u1))
      (flvs! vs idx2 (fl* rho_inv u2))
      (flvs! ws idx2 (fl* rho_inv u3))
      (flvs! square idx2 sq)
      (flvs! qs idx2 (fl* rho_inv sq))
      (flvs! speed idx2 (sqrt (fl* c1c2 rho_inv (fl- u4 sq))))

      (for* ([m (in-range 5)])
        (flvs! rhs (+ m idx) (flvr forcing (+ m idx)))))))

  (define-syntax-rule (KZERO a ...) 0)
  (define-syntax-rule (KIDENT a ...) (begin a ...))

  (define-syntax-case (DISSIP A)
    (with-syntax-values ([(nkk2 njj2 nii2) (PICK3 #'A (nz2 ny2 nx2) (nz2 nx2 ny2) (nx2 ny2 nz2))]
                         [(kk jj ii) (PICK3 #'A (k j i) (k i j) (i j k))]
                         [(kksize1 jjsize1 iisize1) (PICK3 #'A (ksize1 jsize1 isize1)
                                                               (ksize1 isize1 jsize1)
                                                               (isize1 jsize1 ksize1))]
                         [(iiiisize2 jjjjsize2 kkkksize2) (PICK3 #'A (i (fx* j jsize2) (fx* k ksize2))
                                                                     ((fx* j jsize2) i (fx* k ksize2))
                                                                     ((fx* k ksize2) (fx* j jsize2) i))]
                         [(d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1) (PICK3 #'A
                                                                 (dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1)
                                                                 (dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1)
                                                                 (dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1))]
                         [(_con_0 _con_1 _con_2 __con3 __con4 __con5) (PICK3 #'A
                                                                 ((fl* xxcon2 con43) xxcon2 xxcon2 xxcon3 xxcon4 xxcon5)
                                                                 (yycon2 (fl* yycon2 con43) yycon2 yycon3 yycon4 yycon5)
                                                                 (zzcon2 zzcon2 (fl* zzcon2 con43) zzcon3 zzcon4 zzcon5))]
                         [(t_2 t_2m1 t_2m2 t_2m3) (PICK3 #'A
                                                     (tx2 KIDENT KZERO KZERO)
                                                     (ty2 KZERO KIDENT KZERO)
                                                     (tz2 KZERO KZERO KIDENT))]
                         [(zs iisize2) (PICK3 #'A (us 1) (vs jsize2) (ws ksize2))])

    #'(CGfor cg ([kk (in-range 1 (add1 nkk2))])
      (for ([jj (in-range 1 (add1 njj2))])
        (let* ([JKIDX  (fx+ (fx* jj jjsize1) (fx* kk kksize1))]
               [JKIDX2 (fx+ jjjjsize2 kkkksize2)])
        (for ([ii (in-range 1 (add1 nii2))])
          (let* ([idx    (fx+ (fx* ii iisize1) JKIDX)]
                 [idx2   (fx+ iiiisize2 JKIDX2)]
                 [idx4   (fx+ idx 4)]
                 [idxz+  (fx+ idx iisize1)]
                 [idxz-  (fx- idx iisize1)]
                 [idxz+4 (fx+ idxz+ 4)]
                 [idxz-4 (fx+ idxz- 4)]
                 [idxz+2 (fx+ idxz+ iisize1)]
                 [idxz-2 (fx- idxz- iisize1)]
                 [idx2z+ (fx+ idx2 iisize2)]
                 [idx2z- (fx- idx2 iisize2)])
            (define-syntax-rule (citA C C1 C2 C3) (fl* C (fl+ (fl- C1 (fl* 2.0 C2)) C3)))
            (define-syntax-rule (cit3S C V I1 I2 I3) (citA C (flvr V I1) (flvr V I2) (flvr V I3)))
            (define-syntax-rule (cit2 C V)  (citA C (flvr V idx2z+) (flvr V idx2) (flvr V idx2z-)))
            (define-syntax-rule (cit3 C F V)(citA C (F (flvr V idx2z+)) (F (flvr V idx2)) (F (flvr V idx2z-))))
            (define-syntax-rule (t_2m)
              (fl* (fl+    (flvr u idxz+4) (fl-  (flvr square idx2z+))
                    (fl- (flvr u idxz-4))    (flvr square idx2z-)) c2))
            (define-syntax-rule (t_2it l r) (fl- (fl* t_2 (fl-  l r))))
;internal error: broken unbox depth
            (define-syntax-rule (t_2it3 l r o) (fl- (fl* t_2 (+ (-  l r) o))))
            (define-syntax-rule (t_2itlr UI ZSI) (fl* (flvr u UI) (flvr zs ZSI)))

            (define-syntax-rule (mid d__t_1 __con2X ZS AA t_2m_)
              (let ([idxA (fx+ idx AA)]
                    [idxz+A (fx+ idxz+ AA)]
                    [idxz-A (fx+ idxz- AA)])
                (flvs!+ rhs (fx+ idx AA)
                  (cit3S d__t_1 u idxz+A idxA idxz-A)
                  (cit2 __con2X ZS)
                  (t_2it3 (t_2itlr idxz+A idx2z+)
                          (t_2itlr idxz-A idx2z-)
                           (t_2m_ (t_2m))))))
        
            (flvs!+ rhs idx
              (cit3S d_1t_1 u idxz+ idx idxz-)
              (t_2it (flvr u (fx+ idxz+ A)) 
                     (flvr u (fx+ idxz- A))))

            (mid d_2t_1 _con_0 us 1 t_2m1)
            (mid d_3t_1 _con_1 vs 2 t_2m2)
            (mid d_4t_1 _con_2 ws 3 t_2m3)

            (define-syntax-rule (CONS5W UI RI) (fl* (flvr u UI) (flvr rho_i RI)))
            (define-syntax-rule (T_25W UI RI) (fl* (fl- (fl* c1 (flvr u UI))
                                                        (fl* c2 (flvr square RI)))
                                                   (flvr zs RI)))
            (flvs!+ rhs idx4
              (cit3S d_5t_1 u idxz+4 idx4 idxz-4)
              (cit2 __con3 qs)
              (cit3 __con4 sqr zs)
              (citA __con5 (CONS5W idxz+4 idx2z+)
                           (CONS5W idx4 idx2)
                           (CONS5W idxz-4 idx2z-))
              (t_2it (T_25W idxz+4 idx2z+)
                     (T_25W idxz-4 idx2z-)))))

        (fourth-order-dissipation ii nii2 rhs u i j k m midx idx 
        (fx+ m JKIDX (fx* ii iisize1))
        (fx+ m JKIDX (fx* ii iisize1) iisize1)
        (fx+ m JKIDX (fx* ii iisize1) iisize1 iisize1)
        (fx+ m JKIDX (fx- (fx* ii iisize1) iisize1))
        (fx+ m JKIDX (fx- (fx* ii iisize1) iisize1 iisize1))
        midx
        (fx+ JKIDX (fx* ii iisize1))
         dssp))))))

  (CG-B cg) 

  (DISSIP 1)

  (CG-B cg) 

  (DISSIP 2)

  (CG-B cg) 

  (DISSIP 3)

  (CG-B cg)

  (CGfor cg ([k (in-range 1 (add1 nz2))])
    (for* ([j (in-range 1 (add1 ny2))]
           [i (in-range 1 (add1 nx2))]
           [m (in-range 5)])
    (let* ([idx (fx+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1))])
      (flvs!* rhs idx dt)))))
  

(define (txinvr cg nz2 ny2 nx2 isize1 jsize1 ksize1 jsize2 ksize2 rho_i us vs ws qs speed rhs c2 bt)
  (CGfor cg ([k (in-range 1 (fx++ nz2))])
    (let ([ks1 (fx* k ksize1)]
          [ks2 (fx* k ksize2)])
    (for ([j (in-range 1 (fx++ ny2))])
    (let ([js1 (fx+ (fx* j jsize1) ks1)]
          [js2 (fx+ (fx* j jsize2) ks2)])
    (for ([i (in-range 1 (fx++ nx2))])
    (let* ([idx (fx+ (fx* i isize1) js1)]
           [idx2 (fx+ i js2)]

           [ru1 (flvr rho_i idx2)]
           [uu (flvr us idx2)]
           [vv (flvr vs idx2)]
           [ww (flvr ws idx2)]
           [ac (flvr speed idx2)]
           [ac2inv (fl/ 1.0 (sqr ac))]

           [r1 (flvr rhs (fx+ 0 idx))]
           [r2 (flvr rhs (fx+ 1 idx))]
           [r3 (flvr rhs (fx+ 2 idx))]
           [r4 (flvr rhs (fx+ 3 idx))]
           [r5 (flvr rhs (fx+ 4 idx))]
           [t1 (fl* c2 ac2inv (fl+ (fl- (fl* (flvr qs idx2) r1) (fl* uu r2) (fl* vv r3) (fl* ww r4)) r5))]
           [t2 (fl* bt ru1 (fl- (fl* uu r1) r2))]
           [t3 (fl* bt ru1 ac t1)])

      (flvs! rhs (fx+ 0 idx) (fl- r1 t1))
      (flvs! rhs (fx+ 1 idx) (fl- (fl* ru1 (fl- (fl* ww r1) r4))))
      (flvs! rhs (fx+ 2 idx)    (fl* ru1 (fl- (fl* vv r1) r3)))
      (flvs! rhs (fx+ 3 idx) (fl- t3 t2))
      (flvs! rhs (fx+ 4 idx) (fl+ t2 t3)))))))))

(define (tzetar cg u nz2 ny2 nx2 isize1 jsize1 ksize1 jsize2 ksize2 us vs ws qs speed rhs bt c2iv)
  (CGfor cg ([k (in-range 1 (fx++ nz2))])
    (for* ([j (in-range 1 (fx++ ny2))]
           [i (in-range 1 (fx++ nx2))])
    (let* ([idx (fx+ (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
           [idx2 (fx+ i (fx* j jsize2) (fx* k ksize2))]

           [xvel (flvr us idx2)]
           [yvel (flvr vs idx2)]
           [zvel (flvr ws idx2)]
           [ac (flvr speed idx2)]
           [ac2u (sqr ac)]

           [r1 (flvr rhs (fx+ 0 idx))]
           [r2 (flvr rhs (fx+ 1 idx))]
           [r3 (flvr rhs (fx+ 2 idx))]
           [r4 (flvr rhs (fx+ 3 idx))]
           [r5 (flvr rhs (fx+ 4 idx))]

           [uzik1 (flvr u idx)]
           [btuz  (fl* bt uzik1)]

           [t1 (fl* (fl/ btuz ac) (fl+ r4 r5))]
           [t2 (fl+ r3 t1)]
           [t3 (fl* btuz (fl- r4 r5))])
      (flvs! rhs (fx+ 0 idx) t2)
      (flvs! rhs (fx+ 1 idx) (fl+ (fl- (fl* uzik1 r2)) (fl* xvel t2)))
      (flvs! rhs (fx+ 2 idx) (fl+ (fl* uzik1 r1) (fl* yvel t2))) 
      (flvs! rhs (fx+ 3 idx) (fl+ (fl* zvel t2) t3))
      (flvs! rhs (fx+ 4 idx) (fl+ (fl* uzik1 (fl- (fl* yvel r1) (fl* xvel r2)))
                              (fl* (flvr qs idx2) t2)
                              (fl* c2iv ac2u t1)
                              (fl* zvel t3)))))))
(define (verify class no_time_steps dt compute_rhs_thunk
  nx2 ny2 nz2 isize1 jsize1 ksize1 u rhs dnzm1 dnym1 dnxm1)
  (define xcrdif (make-flvector 5 0.0))
  (define xcedif (make-flvector 5 0.0))
  (define xcr (make-flvector 5 0.0))
  (define xce (make-flvector 5 0.0))
  (define-values (xcrref xceref dtref) (get-verify-values class))

;;;//---------------------------------------------------------------------
;;;//   compute the error norm and the residual norm, and exit if not printing
;;;//---------------------------------------------------------------------
    (error-norm xce nx2 ny2 nz2 isize1 jsize1 ksize1 u dnzm1 dnym1 dnxm1)
    (compute_rhs_thunk)
    (rhs-norm xcr nz2 ny2 nx2 isize1 jsize1 ksize1 rhs)

    (for ([m (in-range 5)]) 
      (flvs!/ xcr m dt))

;;;//---------------------------------------------------------------------
;;;//    reference data for 12X12X12 grids after 100 time steps, with DT = 1.50d-02
;;;//---------------------------------------------------------------------
;;;//---------------------------------------------------------------------
;;;//    Reference values of RMS-norms of residual.
;;;//---------------------------------------------------------------------
;;;//---------------------------------------------------------------------
;;;//    Reference values of RMS-norms of solution error.
;;;//---------------------------------------------------------------------

    (for ([m (in-range 5)]) 
      (let ([xcrr (flvr xcrref m)]
            [xcer (flvr xceref m)])
        (flvs! xcrdif m (abs (/ (- (flvr xcr m) xcrr) xcrr)))
        (flvs! xcedif m (abs (/ (- (flvr xce m) xcer) xcer)))))

  (define epsilon 1.0E-8)
  (begin0
    (if (not (equal? class #\U))
      (let ([verified (and ((abs (- dt dtref)) . <= . epsilon)
                           (<epsilon-vmap xcrdif epsilon)
                           (<epsilon-vmap xcedif epsilon))])
        (printf " Verification being performed for class ~a\n" class)
        (printf " Accuracy setting for epsilon = ~a\n" epsilon)
        (unless verified
          (printf "DT does not match the reference value of ~a\n" dtref))
        verified)
      (begin
        (printf " Unknown CLASS")
        (printf " RMS-norms of residual")
        #f))

    (printf " Comparison of RMS-norms of residual\n")
    (for ([m (in-range (flvector-length xcr))])
      (printf "~a. ~a ~a ~a\n" m (flvr xcr m) (flvr xcrref m) (flvr xcrdif m)))
    (printf " Comparison of RMS-norms of solution error\n")
    (for ([m (in-range (flvector-length xce))])
      (printf "~a. ~a ~a ~a\n" m (flvr xce m) (flvr xceref m) (flvr xcedif m)))))

(define-syntax-case (__solve cg A kk jj ii nkk2 njj2 nii2 i j k inverter
_s rho_ rhs lhs speed lhsp lhsm
d__ d_5 d_1
dtt_1 dtt_2
c2dtt_
isize1 jsize1 ksize1 jsize2 ksize2 jsize4
rho_i cv con43 c1c5 c3c4 d_max bt
comz1 comz4 comz5 comz6
)
    (with-syntax-values (
;;;[(nkk2 njj2 nii2) (PICK3 #'A (nz2 ny2 nx2) (nz2 nx2 ny2) (nx2 ny2 nz2))]
;;;                         [(kk jj ii) (PICK3 #'A (k j i) (k i j) (i j k))]
                         [(kksize1 jjsize1 iisize1) (PICK3 #'A (ksize1 jsize1 isize1)
                                                               (ksize1 isize1 jsize1)
                                                               (isize1 jsize1 ksize1))]
                         [(iiiisize2 jjjjsize2 kkkksize2) (PICK3 #'A (i (fx* j jsize2) (fx* k ksize2))
                                                                     ((fx* j jsize2) i (fx* k ksize2))
                                                                     ((fx* k ksize2) (fx* j jsize2) i))]
                         [(zs iisize2) (PICK3 #'A (us 1) (vs jsize2) (ws ksize2))])

  #'(begin
  (CGfor cg ([kk (in-range 1 (add1 nkk2))])
    (for ([jj (in-range 1 (add1 njj2))])
      (let ([JKIDX (fx+ (fx* jj jjsize1) (fx* kk kksize1))]
            [JKIDX2 (fx+ jjjjsize2 kkkksize2)])
      (for ([ii (in-range (fx+ nii2 2))])
        (let* ([idx2 (fx+ i (fx* j jsize2) (fx* k ksize2))]
               [ru1 (fl* c3c4 (flvr rho_i idx2))])
          (flvs! cv ii (flvr _s idx2))
          (flvs! rho_ ii (flmax* (fl+ d__ (fl* con43 ru1))
                                 (fl+ d_5 (fl* c1c5  ru1))  
                                 (fl+ d_max ru1)  
                                 d_1))))
      (lhsinit (add1 nii2) jsize4 lhs lhsp lhsm)

      (for ([ii (in-range 1 (fx+ nii2 1))])
        (let ([iij (fx* ii jsize4)]
              [ii+ (add1 ii)]
              [ii- (sub1 ii)])
        (flvs! lhs (fx+ 0 iij) 0.0)
        (flvs! lhs (fx+ 1 iij) (fl- (- (fl* dtt_2 (flvr cv ii-))) (fl* dtt_1 (flvr rho_ ii-))))
        (flvs! lhs (fx+ 2 iij) (fl+ 1.0 (fl* c2dtt_ (flvr rho_ ii))))
        (flvs! lhs (fx+ 3 iij)    (fl- (fl* dtt_2 (flvr cv ii+)) (fl* dtt_1 (flvr rho_ ii+))))
        (flvs! lhs (fx+ 4 iij) 0.0)))
;;;//---------------------------------------------------------------------
;;;//      add fourth order dissipation                             
;;;//---------------------------------------------------------------------
      (let* ([ii 1]
            [iij (fx* ii jsize4)]
            [iij+ (fx* (add1 ii) jsize4)])
        (flvs!+ lhs (fx+ 2 iij) comz5)
        (flvs!- lhs (fx+ 3 iij) comz4)
        (flvs!+ lhs (fx+ 4 iij) comz1)

        (flvs!- lhs (fx+ 1 iij+) comz4)
        (flvs!+ lhs (fx+ 2 iij+) comz6)
        (flvs!- lhs (fx+ 3 iij+) comz4)
        (flvs!+ lhs (fx+ 4 iij+) comz1))
    
      (for ([ii (in-range 3 (sub1 nii2))])
        (let ([iij (fx* ii jsize4)])
          (flvs!+ lhs (fx+ 0 iij) comz1)
          (flvs!- lhs (fx+ 1 iij) comz4)
          (flvs!+ lhs (fx+ 2 iij) comz6)
          (flvs!- lhs (fx+ 3 iij) comz4)
          (flvs!+ lhs (fx+ 4 iij) comz1)))

      (let* ([ii (sub1 nii2)]
            [iij (fx* ii jsize4)]
            [iij+ (fx* (add1 ii) jsize4)])
        (flvs!+ lhs (fx+ 0 iij) comz1)
        (flvs!- lhs (fx+ 1 iij) comz4)
        (flvs!+ lhs (fx+ 2 iij) comz6)
        (flvs!- lhs (fx+ 3 iij) comz4)

        (flvs!+ lhs (fx+ 0 iij+) comz1)
        (flvs!- lhs (fx+ 1 iij+) comz4)
        (flvs!+ lhs (fx+ 2 iij+) comz5))
      
;;;//---------------------------------------------------------------------
;;;//      subsequently, fill the other factors (u+c), (u-c) by adding to 
;;;//      the first  
;;;//---------------------------------------------------------------------
      (for ([ii (in-range 1 (add1 nii2))])
        (let* ([iij (fx* ii jsize4)]
               [idx2z (fx+ JKIDX2 (fx* ii iisize2))]
               [d+ (fl* dtt_2 (flvr speed (fx+ idx2z iisize2)))]
               [d- (fl* dtt_2 (flvr speed (fx- idx2z iisize2)))])
          (flvs! lhsp (fx+ 0 iij)    (flvr lhs (fx+ 0 iij)))
          (flvs! lhsp (fx+ 1 iij) (fl- (flvr lhs (fx+ 1 iij)) d-))
          (flvs! lhsp (fx+ 2 iij)      (flvr lhs (fx+ 2 iij)))
          (flvs! lhsp (fx+ 3 iij) (fl+ (flvr lhs (fx+ 3 iij)) d+))
          (flvs! lhsp (fx+ 4 iij)      (flvr lhs (fx+ 4 iij)))

          (flvs! lhsm (fx+ 0 iij)      (flvr lhs (fx+ 0 iij)))
          (flvs! lhsm (fx+ 1 iij) (fl+ (flvr lhs (fx+ 1 iij)) d-))
          (flvs! lhsm (fx+ 2 iij)      (flvr lhs (fx+ 2 iij)))
          (flvs! lhsm (fx+ 3 iij) (fl- (flvr lhs (fx+ 3 iij)) d+))
          (flvs! lhsm (fx+ 4 iij)      (flvr lhs (fx+ 4 iij)))))

;;;//---------------------------------------------------------------------
;;;//                          FORWARD ELIMINATION  
;;;//---------------------------------------------------------------------
;;;
;;;//---------------------------------------------------------------------
;;;//      perform the Thomas algorithm; first, FORWARD ELIMINATION     
;;;//---------------------------------------------------------------------
      (for ([ii (in-range nii2)])
        (let* ([iij (fx* ii jsize4)]
               [ii1j (fx* (add1 ii) jsize4)]
               [ii2j (fx* (fx+ ii 2) jsize4)]
               [fac1 (fl/ 1.0 (flvr lhs (fx+ 2 iij)))]
               [idx (fx+ JKIDX (fx* ii iisize1))]
               [idx1 (fx+ idx iisize1)]
               [idx2 (fx+ idx1 iisize1)])
          (flvs!* lhs (fx+ 3 iij) fac1)
          (flvs!* lhs (fx+ 4 iij) fac1)
          (for ([m (in-range 3)])
            (let ([midx (fx+ m idx)])
              (flvs!* rhs midx fac1)))
          (flvs!- lhs (fx+ 2 ii1j) (fl* (flvr lhs (fx+ 1 ii1j)) (flvr lhs (fx+ 3 iij))))
          (flvs!- lhs (fx+ 3 ii1j) (fl* (flvr lhs (fx+ 1 ii1j)) (flvr lhs (fx+ 4 iij))))
          (for ([m (in-range 3)])
            (let ([midx (fx+ m idx)]
                  [midx1 (fx+ m idx1)])
              (flvs!- rhs midx1 (fl* (flvr lhs (fx+ 1 ii1j)) (flvr rhs midx)))))
          (flvs!- lhs (fx+ 1 ii2j) (fl* (flvr lhs ii2j) (flvr lhs (fx+ 3 iij))))
          (flvs!- lhs (fx+ 2 ii2j) (fl* (flvr lhs ii2j) (flvr lhs (fx+ 4 iij))))
          (for ([m (in-range 3)])
            (let ([midx (fx+ m idx)]
                  [midx2 (fx+ m idx2)])
              (flvs!- rhs midx2 (fl* (flvr lhs ii2j) (flvr rhs midx)))))))

;;;//---------------------------------------------------------------------
;;;//      The last two rows in this grid block are a bit different, 
;;;//      since they do not have two more rows available for the
;;;//      elimination of off-diagonal entries
;;;//---------------------------------------------------------------------
      (let* ([ii nii2]
             [iij (fx* ii jsize4)]
             [ii1j (fx* (add1 nii2) jsize4)]
             [fac1 (fl/ 1.0 (flvr lhs (fx+ 2 iij)))]
             [idx (fx+ JKIDX (fx* ii iisize1))]
             [idx1 (fx+ idx iisize1)])
        (flvs!* lhs (fx+ 3 iij) fac1)  
        (flvs!* lhs (fx+ 4 iij) fac1)  

        (for ([m (in-range 3)])
          (let ([midx (fx+ m idx)])
            (flvs!* rhs midx fac1)))

        (flvs!- lhs (fx+ 2 ii1j) (fl* (flvr lhs (fx+ 1 ii1j)) (flvr lhs (fx+ 3 iij))))
        (flvs!- lhs (fx+ 3 ii1j) (fl* (flvr lhs (fx+ 1 ii1j)) (flvr lhs (fx+ 4 iij))))

        (for ([m (in-range 3)])
          (let ([midx (fx+ m idx)]
                [midx1 (fx+ m idx1)])
            (flvs!- rhs midx1 (fl* (flvr lhs (fx+ 1 ii1j)) (flvr rhs midx)))))

;;;//---------------------------------------------------------------------
;;;//            scale the last row immediately 
;;;//---------------------------------------------------------------------
        (let ([fac2 (/ 1.0 (flvr lhs (fx+ 2 ii1j)))]
              [idx1 (fx+ idx iisize1)])
            (for ([m (in-range 3)])
            (let ([midx1 (fx+ m idx1)])
              (flvs!* rhs midx1 fac2)))))

;;;//---------------------------------------------------------------------
;;;//      do the u+c and the u-c factors                 
;;;//---------------------------------------------------------------------
      (for ([ii (in-range nii2)])
        (let* ([iij (fx* ii jsize4)]
               [ii1j (fx* (add1 ii) jsize4)]
               [ii2j (fx* (fx+ ii 2) jsize4)]
               [idx (fx+ JKIDX (fx* ii iisize1))]
               [idx1 (fx+ idx iisize1)]
               [idx2 (fx+ idx1 iisize1)])
          (let* ([fac1 (fl/ 1.0 (flvr lhsp (+ 2 iij)))]
                 [m 3]
                 [midx (+ m idx)]
                 [midx1 (+ m idx1)]
                 [midx2 (+ m idx2)])
            (flvs!* lhsp (+ 3 iij) fac1)
            (flvs!* lhsp (+ 4 iij) fac1)
            (flvs!* rhs midx fac1)
            (flvs!- lhsp (+ 2 ii1j) (fl* (flvr lhsp (+ 1 ii1j)) (flvr lhsp (+ 3 iij))))
            (flvs!- lhsp (+ 3 ii1j) (fl* (flvr lhsp (+ 1 ii1j)) (flvr lhsp (+ 4 iij))))
            (flvs!- rhs midx1       (fl* (flvr lhsp (+ 1 ii1j)) (flvr rhs midx)))
            (flvs!- lhsp (+ 1 ii2j) (fl* (flvr lhsp ii2j) (flvr lhsp (+ 3 iij))))
            (flvs!- lhsp (+ 2 ii2j) (fl* (flvr lhsp ii2j) (flvr lhsp (+ 4 iij))))
            (flvs!- rhs midx2       (fl* (flvr lhsp ii2j) (flvr rhs midx))))

          (let* ([fac1 (fl/ 1.0 (flvr lhsm (+ 2 iij)))]
                 [m 4]
                 [midx (+ m idx)]
                 [midx1 (+ m idx1)]
                 [midx2 (+ m idx2)])
            (flvs!* lhsm (+ 3 iij) fac1)
            (flvs!* lhsm (+ 4 iij) fac1)
            (flvs!* rhs midx fac1)
            (flvs!- lhsm (+ 2 ii1j) (fl* (flvr lhsm (+ 1 ii1j)) (flvr lhsm (+ 3 iij))))
            (flvs!- lhsm (+ 3 ii1j) (fl* (flvr lhsm (+ 1 ii1j)) (flvr lhsm (+ 4 iij))))
            (flvs!- rhs midx1       (fl* (flvr lhsm (+ 1 ii1j)) (flvr rhs midx)))
            (flvs!- lhsm (+ 1 ii2j) (fl* (flvr lhsm ii2j) (flvr lhsm (+ 3 iij))))
            (flvs!- lhsm (+ 2 ii2j) (fl* (flvr lhsm ii2j) (flvr lhsm (+ 4 iij))))
            (flvs!- rhs midx2       (fl* (flvr lhsm ii2j) (flvr rhs midx))))))

;;;//---------------------------------------------------------------------
;;;//         And again the last two rows separately
;;;//---------------------------------------------------------------------
      (let* ([ii nii2]
             [iij (fx* ii jsize4)]
             [ii1j (fx* (add1 ii) jsize4)]
             [idx (fx+ JKIDX (fx* ii iisize1))]
             [idx1 (fx+ idx iisize1)])
        (let* ([fac1 (/ 1.0 (flvr lhsp (+ 2 iij)))]
               [m 3]
               [midx (fx+ m idx)]
               [midx1 (fx+ m idx1)])
          (flvs!* lhsp (fx+ 3 iij) fac1)
          (flvs!* lhsp (fx+ 4 iij) fac1)
          (flvs!* rhs midx       fac1)
          (flvs!- lhsp (fx+ 2 ii1j) (fl* (flvr lhsp (fx+ 1 ii1j)) (flvr lhsp (fx+ 3 iij))))
          (flvs!- lhsp (fx+ 3 ii1j) (fl* (flvr lhsp (fx+ 1 ii1j)) (flvr lhsp (fx+ 4 iij))))
          (flvs!- rhs midx1       (fl* (flvr lhsp (fx+ 1 ii1j)) (flvr rhs midx))))
        (let* ([fac1 (fl/ 1.0 (flvr lhsm (fx+ 2 iij)))]
               [m 4]
               [midx (fx+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
               [midx1 (fx+ m idx1)])
          (flvs!* lhsp (fx+ 3 iij) fac1)
          (flvs!* lhsm (fx+ 3 iij) fac1)
          (flvs!* lhsm (fx+ 4 iij) fac1)
          (flvs!* rhs midx       fac1)
          (flvs!- lhsm (fx+ 2 ii1j) (fl* (flvr lhsm (fx+ 1 ii1j)) (flvr lhsm (fx+ 3 iij))))
          (flvs!- lhsm (fx+ 3 ii1j) (fl* (flvr lhsm (fx+ 1 ii1j)) (flvr lhsm (fx+ 4 iij))))
          (flvs!- rhs midx1       (fl* (flvr lhsm (fx+ 1 ii1j)) (flvr rhs midx))))

;;;//---------------------------------------------------------------------
;;;//               Scale the last row immediately
;;;//---------------------------------------------------------------------
        (flvs!/ rhs (fx+ 3 idx1) (flvr lhsp (fx+ 2 ii1j)))
        (flvs!/ rhs (fx+ 4 idx1)  (flvr lhsm (fx+ 2 ii1j))))


;;;//---------------------------------------------------------------------
;;;//                         BACKSUBSTITUTION 
;;;//---------------------------------------------------------------------
      (let* ([ii nii2]
             [iij (fx* ii jsize4)]
             [ii1j (fx* (add1 ii) jsize4)]
             [idx (fx+ JKIDX (fx* ii iisize1))]
             [idx1 (fx+ idx iisize1)])
        (for ([m (in-range 3)])
          (let ([midx (fx+ m idx)]
                [midx1 (fx+ m idx1)])
            (flvs!- rhs midx  (fl* (flvr lhs (fx+ 3 ii1j)) (flvr rhs midx1)))))

        (flvs!- rhs (fx+ 3 idx) (fl* (flvr lhsp (fx+ 3 iij)) (flvr rhs (fx+ 3 idx1))))
        (flvs!- rhs (fx+ 4 idx) (fl* (flvr lhsm (fx+ 3 iij)) (flvr rhs (fx+ 4 idx1)))))

;;;//---------------------------------------------------------------------
;;;//      The first three factors
;;;//---------------------------------------------------------------------
      (for ([ii (in-range (fx- nii2 1) -1 -1)])
        (let* ([iij (fx* ii jsize4)]
               [idx (fx+ JKIDX (fx* ii iisize1))]
               [idx1 (fx+ idx iisize1)]
               [idx2 (fx+ idx1 iisize1)])
          (for ([m (in-range 3)])
            (let ([midx (fx+ m idx)])
              (flvs!- rhs midx  (fl* (flvr lhs (fx+ 3 iij)) (flvr rhs (fx+ m idx1)))
                                (fl* (flvr lhs (fx+ 4 iij)) (flvr rhs (fx+ m idx2))))))

          (flvs!- rhs (fx+ 3 idx) (fl* (flvr lhsp (fx+ 3 iij)) (flvr rhs (fx+ 3 idx1)))
                                (fl* (flvr lhsp (fx+ 4 iij)) (flvr rhs (fx+ 3 idx2))))
          (flvs!- rhs (fx+ 4 idx) (fl* (flvr lhsm (fx+ 3 iij)) (flvr rhs (fx+ 4 idx1)))
                                (fl* (flvr lhsm (fx+ 4 iij)) (flvr rhs (fx+ 4 idx2)))))))))
      (CG-B cg)
      inverter)))

; (- nz2 2)  = grid_points[2]-4
; (- nz2 1)  = grid_points[2]-3
; nz2        = grid_points[2]-2
; (+ nz2 1)  = grid_points[2]-1

(define (x_solve cg nz2 ny2 nx2 us rhon 
dx2 dx5 dx1 dttx1 dttx2 c2dttx1
isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
rho_i cv con43 c1c5 c3c4 dxmax bt
rhs lhs speed lhsp lhsm
comz1 comz4 comz5 comz6
)
(define-values (i j k) (values 0 0 0))
  (__solve cg 1 k j i nz2 ny2 nx2 i j k 
    (ninvr cg nz2 ny2 nx2 isize1 jsize1 ksize1 rhs bt)
    us rhon rhs lhs speed lhsp lhsm
    dx2 dx5 dx1 dttx1 dttx2 c2dttx1
    isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
    rho_i cv con43 c1c5 c3c4 dxmax bt
    comz1 comz4 comz5 comz6
))
(define (y_solve cg nz2 ny2 nx2 vs rhoq
dy3 dy5 dy1 dtty1 dtty2 c2dtty1
isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
rho_i cv con43 c1c5 c3c4 dymax bt
rhs lhs speed lhsp lhsm
comz1 comz4 comz5 comz6
)
(define-values (i j k) (values 0 0 0))
  (__solve cg 2 k i j nz2 ny2 nx2 i j k 
    (pinvr cg nz2 ny2 nx2 isize1 jsize1 ksize1 rhs bt)
    vs rhoq rhs lhs speed lhsp lhsm
    dy3 dy5 dy1 dtty1 dtty2 c2dtty1
    isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
    rho_i cv con43 c1c5 c3c4 dymax bt
    comz1 comz4 comz5 comz6
))

(define (z_solve cg nz2 ny2 nx2 ws rhos
dz4 dz5 dz1 dttz1 dttz2 c2dttz1
isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
rho_i cv con43 c1c5 c3c4 dzmax bt 
rhs lhs speed lhsp lhsm
comz1 comz4 comz5 comz6
u us vs qs
)
(define c2iv 2.5)
(define-values (i j k) (values 0 0 0))
  (__solve cg 3 i j k ny2 nx2 nz2 i j k 
    (tzetar cg u nz2 ny2 nx2 isize1 jsize1 ksize1 jsize2 ksize2 us vs ws qs speed rhs bt c2iv)
    ws rhos rhs lhs speed lhsp lhsm
    dz4 dz5 dz1 dttz1 dttz2 c2dttz1
    isize1 jsize1 ksize1 jsize2 ksize2 jsize4 
    rho_i cv con43 c1c5 c3c4 dzmax bt
    comz1 comz4 comz5 comz6
))

(define (checkSum arr nz2 ny2 nx2 isize1 jsize1 ksize1)
  (for*/fold ([csum 0.0]) ([k (in-range (add1 nz2))]
         [j (in-range (add1 ny2))]
         [i (in-range (add1 nx2))]
         [m (in-range 5)])
    (let* ([offset (+ m (* i isize1) (* j jsize1) (* k ksize1))]
           [arro (flvr arr offset)])
      (+ csum (/ (sqr arro) (* (+ 2 nz2) (+ 2 ny2) (+ 2 nx2) 5))))))

(define (get-input-pars maxlevel)
  (define fn "sp.input")
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
