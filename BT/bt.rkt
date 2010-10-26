#lang racket/base

(provide main)
  
(require "../bm-args.rkt") 
(require "../bm-results.rkt") 
(require "../rand-generator.rkt")
(require "../timer.rkt")
(require "../parallel-utils.rkt")
(require "../debug.rkt")
(require racket/match)
(require racket/math)
(require (for-syntax scheme/base))

(require (only-in scheme/flonum make-flvector 
                                make-shared-flvector 
                                shared-flvector
                                flvector-length
                                flvector
                                flmax
                                flvector-set!
                                flvector-ref
                                fl+
                                fl-
                                fl*
                                fl/)
          (rename-in scheme/flonum [flvector-ref fr] 
                                   [flvector-set! f!]))
(require (only-in scheme/fixnum fx+ 
                                fx-
                                fx*))

#|
(require (only-in scheme/flonum make-flvector 
                                make-shared-flvector 
                                shared-flvector
                                flvector-length
                                flvector
                                flmax))
(require (rename-in scheme/unsafe/ops
                    [unsafe-vector-ref vr] 
                    [unsafe-vector-set! vs!]
                    [unsafe-flvector-ref fr] 
                    [unsafe-flvector-set! f!]
                    [unsafe-fl+ fl+]
                    [unsafe-fl- fl-]
                    [unsafe-fl* fl*]
                    [unsafe-fl/ fl/]
                    [unsafe-fx+ fx+]
                    [unsafe-fx- fx-]
                    [unsafe-fx* fx*]
))
|#
 
(define-syntax-rule (vidx3 i1 i2 i3 n1 n2) (fx+ i1 (fx* n1 (fx+ i2 (fx* n2 i3)))))
(define-syntax-rule (vr3 v i1 i2 i3 n1 n2) (vr v (vidx3 i1 i2 i3 n1 n2)))
(define-syntax-rule (vidx off i1 i2 i3 n1 n2) (fx+ off (vidx3 i1 i2 i3 n1 n2)))
(define-syntax-rule (vro3 v off i1 i2 i3 n1 n2) (vr v (fx+ off (vidx3 i1 i2 i3 n1 n2))))

(define-syntax-rule (f!+ v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (+ (fr v idx) val ...))))
(define-syntax-rule (f!- v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (- (fr v idx) val ...))))
(define-syntax-rule (f!* v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (fl* (fr v idx) val ...))))
(define-syntax-rule (f!/ v idx_ val ...)
  (let ([idx idx_])
    (f! v idx (fl/ (fr v idx) val ...))))
(define-syntax (flmax* stx)
  (syntax-case stx ()
    [(_ a b) #'(flmax a b)]
    [(_ a b ...) #'(flmax a (flmax* b ...))]))

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values  12  0.01    60)]
    [(#\W) (values  24  0.0008 200)]
    [(#\A) (values  64  0.0008 200)]
    [(#\B) (values  102 0.0003 200)]
    [(#\C) (values  162 0.0001 200)]
    [else (error "Unknown class")]))

(define (get-verify-values class)
  (case class
    [(#\S) (values
        (flvector
       1.7034283709541311E-01
       1.2975252070034097E-02
       3.2527926989486055E-02
       2.6436421275166801E-02
       1.9211784131744430E-01)
        (flvector
       4.9976913345811579E-04
       4.5195666782961927E-05
       7.3973765172921357E-05
       7.3821238632439731E-05
       8.9269630987491446E-04)
       0.01)]
    [(#\W) (values
        (flvector
       0.1125590409344E+03
       0.1180007595731E+02
       0.2710329767846E+02
       0.2469174937669E+02
       0.2638427874317E+03)

        (flvector
       0.4419655736008E+01
       0.4638531260002
       0.1011551749967E+01
       0.9235878729944
       0.1018045837718E+02)
       0.0008)]
    [(#\A) (values
        (flvector
       1.0806346714637264E+02
       1.1319730901220813E+01
       2.5974354511582465E+01
       2.3665622544678910E+01
       2.5278963211748344E+02)
        (flvector
       4.2348416040525025
       4.4390282496995698E-01
       9.6692480136345650E-01
       8.8302063039765474E-01)
       9.7379901770829278
       0.0008)]
    [(#\B) (values
        (flvector
       1.4233597229287254E+03
       9.9330522590150238E+01
       3.5646025644535285E+02
       3.2485447959084092E+02
       3.2707541254659363E+03)
        (flvector
       5.2969847140936856E+01
       4.4632896115670668
       1.3122573342210174E+01
       1.2006925323559144E+01
       1.2459576151035986E+02)
       0.0003)]
    [(#\C) (values
        (flvector
       0.62398116551764615E+04
       0.50793239190423964E+03
       0.15423530093013596E+04
       0.13302387929291190E+04
       0.11604087428436455E+05)
        (flvector
       0.16462008369091265E+03
       0.11497107903824313E+02
       0.41207446207461508E+02
       0.37087651059694167E+02
       0.36211053051841265E+03)
       0.0001)]
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

(define make-fxvector make-vector)

(define (run-benchmark args) 
  (let ([bmname "BT"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(problem_size dt_default niter_default) (get-class-size CLASS)])
    (let* (
          [niter niter_default]
          ;[niter 10]
          [dt dt_default]
          [IMAX problem_size]
          [JMAX problem_size]
          [KMAX problem_size]
          ;[grid_points (make-fxvector 3 problem_size)]
          [nx problem_size]
          [ny problem_size]
          [nz problem_size]
          [nx2 (- problem_size 2)]
          [ny2 (- problem_size 2)]
          [nz2 (- problem_size 2)]
          [dnxm1 (/ 1.0 (- problem_size 1))]
          [dnym1 (/ 1.0 (- problem_size 1))]
          [dnzm1 (/ 1.0 (- problem_size 1))]

          [isize (- problem_size 1)]
          [jsize (- problem_size 1)]
          [ksize (- problem_size 1)]

          [jsize1 (add1 IMAX)]
          [ksize1 (* (add1 IMAX) (add1 JMAX))]

          [isize2 5]
          [jsize2 (* 5 (add1 IMAX))]
          [ksize2 (* 5 (add1 IMAX) (add1 JMAX))]

          [jsize3 (+ problem_size 2)]
          [isize4 5]
          [jsize4 (* 5 5)]
          [ksize4 (* 5 5 3)]
          [s1 (* ksize1 KMAX)]
          [s2 (* ksize2 KMAX)]
          [s3 (* ksize4 (add1 problem_size))]
          [s4 (* jsize4 (add1 problem_size))]
          [us      (make-shared-flvector s1 0.0)]
          [vs      (make-shared-flvector s1 0.0)]
          [ws      (make-shared-flvector s1 0.0)]
          [qs      (make-shared-flvector s1 0.0)]
          [rho_i   (make-shared-flvector s1 0.0)]
          [square  (make-shared-flvector s1 0.0)]

          [u       (make-shared-flvector s2 0.0)]
          [rhs     (make-shared-flvector s2 0.0)]
          [forcing (make-shared-flvector s2 0.0)]

          [lhs     (make-shared-flvector s3 0.0)]
          [fjac    (make-shared-flvector s4 0.0)]
          [njac    (make-shared-flvector s4 0.0)]

          [cv      (make-shared-flvector (+ problem_size 2) 0.0)]
          [cuf     (make-shared-flvector (+ problem_size 2) 0.0)]
          [q       (make-shared-flvector (+ problem_size 2) 0.0)]

          [ue      (make-shared-flvector (* 5 jsize3) 0.0)]
          [buf     (make-shared-flvector (* 5 jsize3) 0.0)]

          [c1      1.4]
          [c2      0.4]
          [c3      0.1]
          [c4      1.0]
          [c5      1.4]
          [c1c2    (* 1.4 0.4)]
          [c1c5    (* 1.4 1.4)]
          [c3c4    (* 0.1 1.0)]
          [c1345   (* c1 c3 c4 c5)]
          [tx1     (/ 1.0 (sqr dnxm1))]
          [ty1     (/ 1.0 (sqr dnym1))]
          [tz1     (/ 1.0 (sqr dnzm1))]
          [tx2     (/ 1.0 (* 2.0 dnxm1))]
          [ty2     (/ 1.0 (* 2.0 dnym1))]
          [tz2     (/ 1.0 (* 2.0 dnzm1))]
          [tx3     (/ 1.0 dnxm1)]
          [ty3     (/ 1.0 dnym1)]
          [tz3     (/ 1.0 dnzm1)]
          [c3c4tx3 (* c3c4 tx3)]
          [c3c4ty3 (* c3c4 ty3)]
          [c3c4tz3 (* c3c4 tz3)]
          [con43   (/ 4.0 3.0)]
          [conz1   (- 1.0 c1c5)]
          [con16   (/ 1.0 6.0)]
          [xxcon1  (* c3c4tx3 con43 tx3)]
          [xxcon2  (* c3c4tx3       tx3)]
          [xxcon3  (* c3c4tx3 conz1 tx3)]
          [xxcon4  (* c3c4tx3 con16 tx3)]
          [xxcon5  (* c3c4tx3 c1c5  tx3)]
          [yycon1  (* c3c4ty3 con43 ty3)]
          [yycon2  (* c3c4ty3       ty3)]
          [yycon3  (* c3c4ty3 conz1 ty3)]
          [yycon4  (* c3c4ty3 con16 ty3)]
          [yycon5  (* c3c4ty3 c1c5  ty3)]
          [zzcon1  (* c3c4tz3 con43 tz3)]
          [zzcon2  (* c3c4tz3       tz3)]
          [zzcon3  (* c3c4tz3 conz1 tz3)]
          [zzcon4  (* c3c4tz3 con16 tz3)]
          [zzcon5  (* c3c4tz3 c1c5  tz3)]
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
          [dtdx1tx1  (* dt dx1tx1)]
          [dtdx2tx1  (* dt dx2tx1)]
          [dtdx3tx1  (* dt dx3tx1)]
          [dtdx4tx1  (* dt dx4tx1)]
          [dtdx5tx1  (* dt dx5tx1)]
          [dtdy1ty1  (* dt dy1ty1)]
          [dtdy2ty1  (* dt dy2ty1)]
          [dtdy3ty1  (* dt dy3ty1)]
          [dtdy4ty1  (* dt dy4ty1)]
          [dtdy5ty1  (* dt dy5ty1)]
          [dtdz1tz1  (* dt dz1tz1)]
          [dtdz2tz1  (* dt dz2tz1)]
          [dtdz3tz1  (* dt dz3tz1)]
          [dtdz4tz1  (* dt dz4tz1)]
          [dtdz5tz1  (* dt dz5tz1)]
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
          [dssp    (* 0.25 (flmax* dx1 dy1 dz1))]
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
          [civ 2.5]
)
(define (compute_rhs_thunk)
(compute_rhs (CGSingle) isize2 jsize2 ksize2 jsize1 ksize1 u us vs ws rho_i square qs 
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
      (printf "No input file inputbt.data, Using compiled defaults\n")
      (printf "Size: ~a X ~a X ~a\n" nx ny nz)
      (printf "Iterations: ~a dt: ~a\n" niter dt)
      (initialize u nx ny nz isize2 jsize2 ksize2 dnxm1 dnym1 dnzm1)
      (exact_rhs nx2 ny2 nz2 isize2 jsize2 ksize2 jsize3 forcing dnxm1 dnym1 dnzm1 ue buf cuf q
        rhs u c1 c2 0.25
        tx2 ty2 tz2
        xxcon1 xxcon2 xxcon3 xxcon4 xxcon5
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        yycon1 yycon2 yycon3 yycon4 yycon5
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        zzcon1 zzcon2 zzcon3 zzcon4 zzcon5
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

      (CGspawn (if serial 0 num-threads) bt-body
        u us vs ws rho_i square qs rhs forcing lhs fjac njac
        nx ny nz 
        nx2 ny2 nz2 
        jsize1 ksize1 
        isize2 jsize2 ksize2 
        isize4 jsize4 ksize4 
        isize jsize ksize
        dnxm1 dnym1 dnzm1
        c1 c2 c1c2 c3c4 c1345 con43 dt dssp niter
        dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
        xxcon2 xxcon3 xxcon4 xxcon5
        dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
        yycon2 yycon3 yycon4 yycon5
        dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
        zzcon2 zzcon3 zzcon4 zzcon5
        tx1 tx2 dtdx1tx1 dtdx2tx1 dtdx3tx1 dtdx4tx1 dtdx5tx1
        ty1 ty2 dtdy1ty1 dtdy2ty1 dtdy3ty1 dtdy4ty1 dtdy5ty1
        tz1 tz2 dtdz1tz1 dtdz2tz1 dtdz3tz1 dtdz4tz1 dtdz5tz1)

      (let* ([verified (verify CLASS niter dt compute_rhs_thunk
        nx2 ny2 nz2 isize2 jsize2 ksize2 u rhs dnzm1 dnym1 dnxm1)])
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

(define (bt-body cg  u us vs ws rho_i square qs rhs forcing lhs_ fjac_ njac_
      nx ny nz 
      nx2 ny2 nz2 
      jsize1 ksize1 
      isize2 jsize2 ksize2 
      isize4 jsize4 ksize4 
      isize jsize ksize
      dnxm1 dnym1 dnzm1
      c1 c2 c1c2 c3c4 c1345 con43 dt dssp niter
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
      xxcon2 xxcon3 xxcon4 xxcon5
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
      yycon2 yycon3 yycon4 yycon5
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
      zzcon2 zzcon3 zzcon4 zzcon5
      tx1 tx2 dtdx1tx1 dtdx2tx1 dtdx3tx1 dtdx4tx1 dtdx5tx1
      ty1 ty2 dtdy1ty1 dtdy2ty1 dtdy3ty1 dtdy4ty1 dtdy5ty1
      tz1 tz2 dtdz1tz1 dtdz2tz1 dtdz3tz1 dtdz4tz1 dtdz5tz1
)
   (define s3 (* ksize4 (add1 nx)))
   (define s4 (* jsize4 (add1 nx)))
   (define lhs     (make-flvector s3 0.0))
   (define fjac    (make-flvector s4 0.0))
   (define njac    (make-flvector s4 0.0))
   (define (adi)
    (compute_rhs cg isize2 jsize2 ksize2 jsize1 ksize1 u us vs ws rho_i square qs 
      c1c2 rhs forcing nx2 ny2 nz2 c1 c2 dssp
      tx2 ty2 tz2 con43 dt
      dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
      xxcon2 xxcon3 xxcon4 xxcon5
      dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
      yycon2 yycon3 yycon4 yycon5
      dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
      zzcon2 zzcon3 zzcon4 zzcon5
    )
    
    (CG-B cg)
  
    (x_solve cg nz2 ny2 nx2
      jsize1 ksize1 
      isize2 jsize2 ksize2 
      isize4 jsize4 ksize4 
      isize jsize ksize isize
      u square rhs lhs fjac njac rho_i qs 
      c1 c2 c3c4 c1345 con43 dt 
      tx1 tx2 dtdx1tx1 dtdx2tx1 dtdx3tx1 dtdx4tx1 dtdx5tx1)

    (CG-B cg)

    (y_solve cg nz2 nx2 ny2
      jsize1 ksize1 
      isize2 jsize2 ksize2 
      isize4 jsize4 ksize4 
      isize jsize ksize jsize
      u square rhs lhs fjac njac rho_i qs 
      c1 c2 c3c4 c1345 con43 dt 
      ty1 ty2 dtdy1ty1 dtdy2ty1 dtdy3ty1 dtdy4ty1 dtdy5ty1)

    (CG-B cg)

    (z_solve cg ny2 nx2 nz2
      jsize1 ksize1 
      isize2 jsize2 ksize2 
      isize4 jsize4 ksize4 
      isize jsize ksize ksize
      u square rhs lhs fjac njac rho_i qs 
      c1 c2 c3c4 c1345 con43 dt 
      tz1 tz2 dtdz1tz1 dtdz2tz1 dtdz3tz1 dtdz4tz1 dtdz5tz1)

    (CG-B cg)

    (add cg nz2 ny2 nx2 jsize1 ksize1 isize2 jsize2 ksize2 u rhs)
    )

;;;//---------------------------------------------------------------------
;;;//      do one time step to touch all code, and reinitialize
;;;//---------------------------------------------------------------------
      (adi)

      (CG-n0-only cg
        (initialize u nx ny nz isize2 jsize2 ksize2 dnxm1 dnym1 dnzm1)

        (timer-start 1))

      (for ([step (in-range 1 (add1 niter))])
        (CG-n0-only cg
          (when (or (zero? (modulo step 20)) (= step 1) (= step niter))
            (printf "Time step ~a\n" step)))
        (adi))

      (CG-n0-only cg
        (timer-stop 1))

)

(define (get-mflops total-time niter nx ny nz)
  (if (not (= total-time 0.0))
    (let* ([n3 (* nx ny nz)]
           [t  (/ (+ nx ny nz) 3.0)])
      (/ (* (+ (* 3478.8 n3)
               (* -17655.7(* t t))
               (* 28023.7 t))
            niter)
          (* total-time 1000000.0)))
      0.0))

(define (exact_solution xi eta zeta dtemp offset)
  (for ([m (in-range 5)])
    (f! dtemp (+ m offset) 
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
                                    (* zeta (fr ce (+ m (* 12 5))))))))))))))
  
(define (initialize u nx ny nz isize1 jsize1 ksize1 dnxm1 dnym1 dnzm1)
  (for* ([k (in-range nz)]
         [j (in-range ny)]
         [i (in-range nx)]
         [m (in-range 5)])
    (let ([midx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
      (f! u midx 1.0)))

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
                        [pxi   (+ (* xi           (fr Pface (+ m (* 0 5) (* 1 15))))
                                  (* (- 1.0 xi)   (fr Pface (+ m (* 0 5) (* 0 15)))))]
                        [peta  (+ (* eta          (fr Pface (+ m (* 1 5) (* 1 15))))
                                  (* (- 1.0 eta)  (fr Pface (+ m (* 1 5) (* 0 15)))))]
                        [pzeta (+ (* zeta         (fr Pface (+ m (* 2 5) (* 1 15))))
                                  (* (- 1.0 zeta) (fr Pface (+ m (* 2 5) (* 0 15)))))])
                    (f! u idx  (+ (- (+ pxi peta pzeta) 
                                        (* pxi peta) 
                                        (* pxi pzeta)
                                        (* peta pzeta))
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
                (f! u idx  (fr temp m))      ;west face
                (f! u idx2 (fr temp2 m)))))) ;east face
            
        (for ([i (in-range nx)])
          (let ([xi (* i dnxm1)])
            (exact_solution xi 0.0 zeta temp 0)
            (exact_solution xi 1.0 zeta temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* i isize1) (* k ksize1))]
                     [idx2 (+ (* j2 jsize1) idx)])
                (f! u idx  (fr temp m))        ;south face
                (f! u idx2 (fr temp2 m)))))))) ;north face
  
    (for ([j (in-range ny)])
      (let ([eta (* j dnym1)])
        (for ([i (in-range nx)])
          (let ([xi (* i dnxm1)])
            (exact_solution xi eta 0.0 temp 0)
            (exact_solution xi eta 1.0 temp2 0)
            (for ([m (in-range 5)])
              (let* ([idx (+ m (* i isize1) (* j jsize1))]
                     [idx2 (+ (* k2 ksize1) idx)])
                (f! u idx  (fr temp m))          ;bottom face
                (f! u idx2 (fr temp2 m)))))))))) ;top face

(define (lhsinit lhs size isize4 jsize4 ksize4)
  (for ([i (in-range 0 (add1 size) size)])
    (let ([ik (* i ksize4)]
          [js41 jsize4]
          [js42 (* 2 jsize4)])
      (for ([m (in-range 5)])
        (for ([n (in-range 5)])
          (let ([idx (+ m (* n isize4) ik)])
            (f! lhs (+ (* 0 jsize4) idx) 0.0)
            (f! lhs (+ (* 1 jsize4) idx) 0.0)
            (f! lhs (+ (* 2 jsize4) idx) 0.0)))
        (f! lhs (+ m (* m isize4) (* 1 jsize4) ik) 1.0)))))

(define (add cg nz2 ny2 nx2 jsize1 ksize1 isize2 jsize2 ksize2 u rhs)
  (CGfor cg ([k (in-range 1 (add1 nz2))])
    (for* ([j (in-range 1 (add1 ny2))]
           [i (in-range 1 (add1 nx2))]
           [m (in-range 5)])
    (let ([idx (+ m (fx* i isize2) (fx* j jsize2) (fx* k ksize2))])
      (f!+ u idx (fr rhs idx))))))

(define-syntax-rule (matvec_sub ablock blkoffst avect avcoffst bvect bvcoffst)
  (for ([i (in-range 5)])
    (f!+ bvect (fx+ i bvcoffst)
      (for/fold ([S 0.0]) ([n (in-range 0 21 5)]
                           [j (in-range 5)])
        (fl- S (fl* (fr ablock (+ i n blkoffst)) (fr avect (+ j avcoffst))))))))

(define-syntax-rule (matmul_sub ablock ablkoffst bblock bblkoffst cblock cblkoffst)
  (for ([nj (in-range 0 21 5)])
    (for ([i (in-range 5)])
      (f! cblock (+ i nj cblkoffst)
        (for/fold ([S (fr cblock (+ i nj cblkoffst))]) ([m (in-range 0 21 5)]
                             [k (in-range 5)])
        (fl- S (fl* (fr ablock (+ i m ablkoffst)) (fr bblock (+ k nj bblkoffst)))))))))
 
(define (binvcrhs lhss lhsoffst c coffst r roffst)
  (define-syntax-rule (V- V I N C J O) (f!- V (+ I N O) (fl* C (fr V (+ J N O)))))
  (define-syntax-rule (lhss- I N C J)  (V- lhss I N C J lhsoffst))
  (define-syntax-rule (c-   I N C J)   (V- c I N C J coffst))
  (define-syntax-rule (r-   I C J)     (f!- r (+ I roffst) (fl* C (fr r (+ J roffst)))))
  (define-syntax-rule (lhss* I N C)    (f!* lhss (+ I N lhsoffst) C))
  (define-syntax-rule (c*   I N C)     (f!* c (+ I N coffst) C))
  (define-syntax-rule (r*   I C)       (f!* r (+ I roffst) C))

  (for ([P (in-range 5)]
        [nP (in-range 0 21 5)]
        [nP1 (in-range 5 26 5)])
    (let ([pivot (/ 1.0 (fr lhss (+ P nP lhsoffst)))])
        (for ([n (in-range nP1 21 5)]) (lhss* P n pivot))
        (for ([n (in-range 0   21 5)]) (c*    P n pivot))
        (r* P pivot))
    (for ([i (in-range 0  5)] #:when (not (= i P)))
      (let ([coeff (fr lhss (+ i nP lhsoffst))])
        (for ([n (in-range nP1 21 5)]) (lhss- i n coeff P))
        (for ([n (in-range 0   21 5)]) (c-    i n coeff P))
        (r- i coeff P)))))

(define (binvrhs lhss lhsoffst r roffst)
  (define-syntax-rule (V- V I N C J O) (f!- V (+ I N O) (fl* C (fr V (+ J N O)))))
  (define-syntax-rule (lhss- I N C J) (V- lhss I N C J lhsoffst))
  (define-syntax-rule (r-   I C J) (f!- r (+ I roffst) (fl* C (fr r (+ J roffst)))))
  (define-syntax-rule (lhss* I N C) (f!* lhss (+ I N lhsoffst) C))
  (define-syntax-rule (r*   I C)   (f!* r (+ I roffst) C))

  (for ([P (in-range 5)]
        [nP (in-range 0 21 5)]
        [nP1 (in-range 5 26 5)])
    (let ([pivot (/ 1.0 (fr lhss (+ P nP lhsoffst)))])
        (for ([n (in-range nP1 21 5)]) (lhss* P n pivot))
        (r* P pivot))
    (for ([i (in-range 0  5)] #:when (not (= i P)))
      (let ([coeff (fr lhss (+ i nP lhsoffst))])
        (for ([n (in-range nP1 21 5)]) (lhss- i n coeff P))
        (r- i coeff P)))))




(define (error-norm rms nx2 ny2 nz2 isize1 jsize1 ksize1 u dnzm1 dnym1 dnxm1)
  (for ([m (in-range 5)]) (f! rms m 0.0))

  (let ([u-exact (make-flvector 5 0.0)])
    (for ([k (in-range (+ nz2 2))])
      (let ([zeta (* k dnzm1)])
        (for ([j (in-range (+ ny2 2))])
          (let ([eta (* j dnym1)])
            (for ([i (in-range (+ nx2 2))])
              (let ([xi (* i dnxm1)])
                (exact_solution xi eta zeta u-exact 0)
                (for ([m (in-range 5)])
                  (let* ([idx (+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
                         [add (fl- (fr u idx) (fr u-exact m))])
                    (f!+ rms m (sqr add)))))))))))

  
  (for ([m (in-range 5)])
    (f! rms m (sqrt (/ (fr rms m) nx2 ny2 nz2)))))

(define (rhs-norm rms nz2 ny2 nx2 isize1 jsize1 ksize1 rhs)
  (for ([m (in-range 5)]) (f! rms m 0.0))

  (for* ([k (in-range 1 (add1 nz2))]
         [j (in-range 1 (add1 ny2))]
         [i (in-range 1 (add1 nx2))]
         [m (in-range 5)])
    (let* ([idx (+ m (fx* i isize1) (fx* j jsize1) (fx* k ksize1))]
           [add (fr rhs idx)])
      (f!+ rms m (sqr add))))
  
  (for ([m (in-range 5)])
    (f! rms m (sqrt (/ (fr rms m) nx2 ny2 nz2)))))

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
        (f!+ V didx (- (fl* dssp (+ (fl*  5.0 (fr V2 midx))
                                     (fl* -4.0 (fr V2 midx+))
                                             (fr V2 midx+2))))))))
  (let* ([ii 2]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx+  MIDX+]
             [midx+2 MIDX+2]
             [midx-  MIDX-]
             [didx   DIDX])
        (f!+ V didx (- (fl* dssp (+ (fl* -4.0 (fr V2 midx-)) 
                                     (fl*  6.0 (fr V2 midx))
                                     (fl* -4.0 (fr V2 midx+))
                                             (fr V2 midx+2))))))))
  (for ([ii (in-range 3 (sub1 nII2))])
    (let ([idx IDX])
      (for ([m (in-range 5)])
        (let* ([midx   MIDX]
               [midx+  MIDX+]
               [midx+2 MIDX+2]
               [midx-  MIDX-]
               [midx-2 MIDX-2]
               [didx   DIDX])
          (f!+ V didx (- (fl* dssp (+         (fr V2 midx-2)
                                       (fl* -4.0 (fr V2 midx-))
                                       (fl*  6.0 (fr V2 midx))
                                       (fl* -4.0 (fr V2 midx+))
                                               (fr V2 midx+2)))))))))
  (let* ([ii (sub1 nII2)]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx+  MIDX+]
             [midx-  MIDX-]
             [midx-2 MIDX-2]
             [didx   DIDX])
        (f!+ V didx (- (* dssp (+         (fr V2 midx-2)
                                     (fl* -4.0 (fr V2 midx-))
                                     (fl*  6.0 (fr V2 midx))
                                     (fl* -4.0 (fr V2 midx+)))))))))
  (let* ([ii nII2]
         [idx IDX])
    (for ([m (in-range 5)])
      (let* ([midx   MIDX]
             [midx-  MIDX-]
             [midx-2 MIDX-2]
             [didx   DIDX])
        (f!+ V didx (- (* dssp (+         (fr V2 midx-2)
                                     (fl* -4.0 (fr V2 midx-))
                                     (fl*  5.0 (fr V2 midx)))))))))))

(define (exact_rhs nx2 ny2 nz2 isize1 jsize1 ksize1 jsize3 forcing dnxm1 dnym1 dnzm1 ue buf cuf q
  rhs u c1 c2 dssp
  tx2 ty2 tz2
  xxcon1 xxcon2 xxcon3 xxcon4 xxcon5
  dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
  yycon1 yycon2 yycon3 yycon4 yycon5
  dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
  zzcon1 zzcon2 zzcon3 zzcon4 zzcon5
  dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1)

  (for* ([k (in-range (+ 2 nz2))]
         [j (in-range (+ 2 ny2))]
         [i (in-range (+ 2 nx2))]
         [m (in-range 5)])
    (let ([idx (+ m (* i isize1) (* j jsize1) (* k ksize1))])
      (f!+ forcing idx 0.0)))

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
                (let ([dtpp (/ 1.0 (fr dtemp 0))])
                  (for ([m (in-range 5)])
                    (let ([idx (+ ii (* m jsize3))]
                          [dtv (fr dtemp m)])
                    (f! ue idx dtv)
                    (f! buf idx (* dtpp dtv)))))
                (let* ([i1j  (+ ii (* 1 jsize3))]
                       [i2j (+ ii (* 2 jsize3))]
                       [i3j (+ ii (* 3 jsize3))]
                       [bufij (fr buf i1j)]
                       [bufi2j (fr buf i2j)]
                       [bufi3j (fr buf i3j)]
                       [ueij (fr ue i1j)]
                       [uei2j (fr ue i2j)]
                       [uei3j (fr ue i3j)]
                       [bufij2 (sqr bufij)]
                       [bufi2j2 (sqr bufi2j)]
                       [bufi3j2 (sqr bufi3j)])
                (f! cuf ii (sqr (fr buf (+ ii (* hh jsize3)))))
                (f! buf ii (+ bufij2 bufi2j2 bufi3j2))
                (f! q ii (* 0.5 (+ (* bufij ueij) (* bufi2j uei2j) (* bufi3j uei3j)))))))

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

              (define-syntax-rule (citA C C1 C2 C3) (fl* C (+ (- C1 (* 2.0 C2)) C3)))
              (define-syntax-rule (citS3 C V I1 I2 I3) (citA C (fr V I1) (fr V I2) (fr V I3)))
              (define-syntax-rule (citS C V) (citA C (fr V ip1) (fr V ii) (fr V im1)))

              (define-syntax-rule (t_2it l r) (- (fl* t_2 (-  l r))))
              (define-syntax-rule (t_2it3 l r o) (- (fl* t_2 (+ (fl- l r) o))))
              (define-syntax-rule (t_2itlr UI ZSI) (* (fr ue UI) (fr buf ZSI)))
              (define-syntax-rule (t_2ito) (- (fl* c2 (fl- (fr ue A4J+) (fr q ip1)))
                                              (fl* c2 (fl- (fr ue A4J-) (fr q im1)))))

              (define-syntax-rule (mid d__t_1 __con2X A t_2itother)
                   (let* ([AJ (+ ii (* A jsize3))]
                          [AJ+ (+ AJ 1)]
                          [AJ- (- AJ 1)])
                (f!+ forcing (+ didx A)
                  (citS3 d__t_1 ue AJ+ AJ AJ-)
                  (citS3 __con2X buf AJ+ AJ AJ-)
                  (t_2it3 (t_2itlr AJ+ idx2+)
                         (t_2itlr AJ- idx2-)
                         (t_2itother (t_2ito))))))

              (f!+ forcing didx 
                   (t_2it (fr ue idx2+)
                          (fr ue idx2-))
                   (citS d_1t_1 ue))

              (mid d_2t_1 _con_0 1 t_2m1)
              (mid d_3t_1 _con_1 2 t_2m2)
              (mid d_4t_1 _con_2 3 t_2m3)

              (define-syntax-rule (t4clause I1 I2 I3)
                (* (fr buf I1) (- (* c1 (fr ue I2)) (* c2 (fr q I3)))))
              (f!+ forcing (+ didx 4)
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
      (f!* forcing (+ m (* i isize1) (* j jsize1) (* k ksize1)) -1.0))
)

(define (compute_rhs cg isize2 jsize2 ksize2 jsize1 ksize1 u us vs ws rho_i square qs 
c1c2 rhs forcing nx2 ny2 nz2 c1 c2 dssp
    tx2 ty2 tz2 con43 dt
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
    xxcon2 xxcon3 xxcon4 xxcon5
   dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
    yycon2 yycon3 yycon4 yycon5
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
    zzcon2 zzcon3 zzcon4 zzcon5
)
  (CGfor cg ([k (in-range 0 (+ nz2 2))])
    (for* ([j (in-range (+ ny2 2))]
           [i (in-range (+ nx2 2))])
    (let* ([idx (+ (* i isize2) (* j jsize2) (* k ksize2))]
           [idx2 (+ i (* j jsize1) (* k ksize1))]
           [rho_inv (/ 1.0 (fr u idx))]
           [u1 (fr u (+ idx 1))]
           [u2 (fr u (+ idx 2))]
           [u3 (fr u (+ idx 3))]
           [u4 (fr u (+ idx 4))]
           [sq (* 0.5 (+ (sqr u1) (sqr u2) (sqr u3)) rho_inv)])
      (f! rho_i idx2 rho_inv)
      (f! us idx2 (fl* rho_inv u1))
      (f! vs idx2 (fl* rho_inv u2))
      (f! ws idx2 (fl* rho_inv u3))
      (f! square idx2 sq)
      (f! qs idx2 (fl* rho_inv sq))

      (for* ([m (in-range 5)])
        (f! rhs (+ m idx) (fr forcing (+ m idx)))))))

  (CG-B cg)

  (define-syntax-rule (KZERO a ...) 0)
  (define-syntax-rule (KIDENT a ...) (begin a ...))

  (define-syntax-rule (DISSIP kk jj ii nkk2 njj2 nii2 u rhs
    i j k A zs
    idxz+ idxz- idxz+2 idxz-2 idx2z+ idx2z-
    d_1t_1 d_2t_1 d_3t_1 d_4t_1 d_5t_1
    _con_0 _con_1 _con_2 __con3 __con4 __con5
    t_2 t_2m1 t_2m2 t_2m3)

    (CGfor cg ([kk (in-range 1 (add1 nkk2))])
      (for ([jj (in-range 1 (add1 njj2))])
        (for ([ii (in-range 1 (add1 nii2))])
          (let* ([idx (+ (* i isize2) (* j jsize2) (* k ksize2))]
                 [idx2 (+ i (* j jsize1) (* k ksize1))]
                 [idx4 (+ idx 4)]
                 [idxz+4 (+ idxz+ 4)]
                 [idxz-4 (+ idxz- 4)])
            (define-syntax-rule (citA C C1 C2 C3) (fl* C (+ (- C1 (fl* 2.0 C2)) C3)))
            (define-syntax-rule (cit3S C V I1 I2 I3) (citA C (fr V I1) (fr V I2) (fr V I3)))
            (define-syntax-rule (cit2 C V)  (citA C (fr V idx2z+) (fr V idx2) (fr V idx2z-)))
            (define-syntax-rule (cit3 C F V)(citA C (F (fr V idx2z+)) (F (fr V idx2)) (F (fr V idx2z-))))
            (define-syntax-rule (t_2m)
              (* (+    (fr u idxz+4) (-  (fr square idx2z+))
                    (- (fr u idxz-4))    (fr square idx2z-)) c2))
            (define-syntax-rule (t_2it l r) (- (* t_2 (-  l r))))
            (define-syntax-rule (t_2it3 l r o) (- (* t_2 (+ (-  l r) o))))
            (define-syntax-rule (t_2itlr UI ZSI) (fl* (fr u UI) (fr zs ZSI)))

            (define-syntax-rule (mid d__t_1 __con2X ZS AA t_2m_)
              (let ([idxA (+ idx AA)]
                    [idxz+A (+ idxz+ AA)]
                    [idxz-A (+ idxz- AA)])
                (f!+ rhs (+ idx AA)
                  (cit3S d__t_1 u idxz+A idxA idxz-A)
                  (cit2 __con2X ZS)
                  (t_2it3 (t_2itlr idxz+A idx2z+)
                          (t_2itlr idxz-A idx2z-)
                           (t_2m_ (t_2m))))))
        
            (f!+ rhs idx
              (cit3S d_1t_1 u idxz+ idx idxz-)
              (t_2it (fr u (+ idxz+ A)) 
                     (fr u (+ idxz- A))))

            (mid d_2t_1 _con_0 us 1 t_2m1)
            (mid d_3t_1 _con_1 vs 2 t_2m2)
            (mid d_4t_1 _con_2 ws 3 t_2m3)

            (define-syntax-rule (CONS5W UI RI) (fl* (fr u UI) (fr rho_i RI)))
            (define-syntax-rule (T_25W UI RI) (fl* (- (fl* c1 (fr u UI))
                                                    (fl* c2 (fr square RI)))
                                                 (fr zs RI)))
            (f!+ rhs (+ idx 4)
              (cit3S d_5t_1 u idxz+4 idx4 idxz-4)
              (cit2 __con3 qs)
              (cit3 __con4 sqr zs)
              (citA __con5 (CONS5W idxz+4 idx2z+)
                           (CONS5W idx4 idx2)
                           (CONS5W idxz-4 idx2z-))
              (t_2it (T_25W idxz+4 idx2z+)
                     (T_25W idxz-4 idx2z-)))))

        (fourth-order-dissipation ii nii2 rhs u i j k m midx idx 
        (+ m idx) (+ m idxz+) (+ m idxz+2) (+ m idxz-) (+ m idxz-2)
        midx (+ (* i isize2) (* j jsize2) (* k ksize2)) dssp))))

  (DISSIP k j i nz2 ny2 nx2 u rhs
    i j k 1 us 
    (+ (* (+ i 1) isize2) (* j jsize2) (* k ksize2)) 
    (+ (* (- i 1) isize2) (* j jsize2) (* k ksize2))
    (+ (* (+ i 2) isize2) (* j jsize2) (* k ksize2)) 
    (+ (* (- i 2) isize2) (* j jsize2) (* k ksize2))
    (+ (+ i 1) (* j jsize1) (* k ksize1)) 
    (+ (- i 1) (* j jsize1) (* k ksize1))
    dx1tx1 dx2tx1 dx3tx1 dx4tx1 dx5tx1
    (* xxcon2 con43) xxcon2 xxcon2 xxcon3 xxcon4 xxcon5
    tx2 KIDENT KZERO KZERO)

  (CG-B cg)

  (DISSIP k i j nz2 nx2 ny2 u rhs
    i j k 2 vs 
    (+ (* i isize2) (* (+ j 1) jsize2) (* k ksize2)) 
    (+ (* i isize2) (* (- j 1) jsize2) (* k ksize2))
    (+ (* i isize2) (* (+ j 2) jsize2) (* k ksize2)) 
    (+ (* i isize2) (* (- j 2) jsize2) (* k ksize2))
    (+ i (* (+ j 1) jsize1) (* k ksize1)) 
    (+ i (* (- j 1) jsize1) (* k ksize1))
    dy1ty1 dy2ty1 dy3ty1 dy4ty1 dy5ty1
    yycon2 (* yycon2 con43) yycon2 yycon3 yycon4 yycon5
    ty2 KZERO KIDENT KZERO)

  (CG-B cg)

  (DISSIP i j k nx2 ny2 nz2 u rhs
    i j k 3 ws 
    (+ (* i isize2) (* j jsize2) (* (+ k 1) ksize2)) 
    (+ (* i isize2) (* j jsize2) (* (- k 1) ksize2))
    (+ (* i isize2) (* j jsize2) (* (+ k 2) ksize2)) 
    (+ (* i isize2) (* j jsize2) (* (- k 2) ksize2))
    (+ i (* j jsize1) (* (+ k 1) ksize1)) 
    (+ i (* j jsize1) (* (- k 1) ksize1))
    dz1tz1 dz2tz1 dz3tz1 dz4tz1 dz5tz1
    zzcon2 zzcon2 (* zzcon2 con43) zzcon3 zzcon4 zzcon5
    tz2 KZERO KZERO KIDENT)

  (CG-B cg)
  (CGfor cg ([k (in-range 1 (add1 nz2))])
    (for* ([j (in-range 1 (add1 ny2))]
           [i (in-range 1 (add1 nx2))]
           [m (in-range 5)])
    (let* ([idx (+ m (* i isize2) (* j jsize2) (* k ksize2))])
      (f!* rhs idx dt)))))

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
      (f!/ xcr m dt))

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
      (let ([xcrr (fr xcrref m)]
            [xcer (fr xceref m)])
        (f! xcrdif m (abs (/ (- (fr xcr m) xcrr) xcrr)))
        (f! xcedif m (abs (/ (- (fr xce m) xcer) xcer)))))

  (define  epsilon 1.0E-8)
  (begin0
    (if (not (equal? class #\U))
      (let ([verified ((abs (- dt dtref)) . <= . epsilon)])
        (printf "Verification being performed for class ~a\n" class)
        (printf "Accuracy setting for epsilon = ~a\n" epsilon)
        (unless verified (printf "DT does not match the reference value of ~a\n" dtref))
        verified)
      (begin
        (printf " Unknown CLASS")
        (printf " RMS-norms of residual")
        -1))
    (printf "Comparison of RMS-norms of residual\n")
    (for ([m (in-range (flvector-length xcr))])
      (printf "~a. ~a ~a ~a\n" m (fr xcr m) (fr xcrref m) (fr xcrdif m)))
    (printf "Comparison of RMS-norms of solution error\n")
    (for ([m (in-range (flvector-length xce))])
      (printf "~a. ~a ~a ~a\n" m (fr xce m) (fr xceref m) (fr xcedif m)))))

(define-syntax-rule (define-syntax-case (N a ...) b ...)
  (define-syntax (N stx)
    (syntax-case stx ()
      [(N a ...) b ...])))

(define-syntax-case (__solve NAME R i j k kk jj ii nkk2 njj2 nii2)

#'(define (NAME cg nkk2 njj2 nii2 jsize1 ksize1 isize2 jsize2 ksize2 isize4 jsize4 ksize4 isize jsize ksize IISIZE u square rhs lhs fjac njac rho_i qs c1 c2 c3c4 c1345 con43 dt t_1 t_2 dtd_1t_1 dtd_2t_1 dtd_3t_1 dtd_4t_1 dtd_5t_1)
  
  (CGfor cg ([kk (in-range 1 (add1 nkk2))])
    (for ([jj (in-range 1 (add1 njj2))])
      (define-syntax-case (MODDET R MACRO N1 N2 N3 V1 V2 V3 BODY (... ...))
        (with-syntax ([LETBINDINGS (case R
                                [(1) #'([N1 (MACRO V1)] [N2 V2] [N3 V3])]
                                [(2) #'([N1 V1] [N2 (MACRO V2)] [N3 V3])]
                                [(3) #'([N1 V1] [N2 V2] [N3 (MACRO V3)])]
                                [else (raise "invalid MODDET value")])])
        #'(let LETBINDINGS
          BODY (... ...)
        )))

      (define-syntax-case (REPDET R N1 N2 N3 V1 V2 V3 VV1 VV2 VV3 BODY (... ...))
        (with-syntax ([LETBINDINGS (case R
                                [(1) #'([N1 VV1] [N2 V2] [N3 V3])]
                                [(2) #'([N1 V1] [N2 VV2] [N3 V3])]
                                [(3) #'([N1 V1] [N2 V2] [N3 VV3])]
                                [else (raise "invalid REPDET value")])])
        #'(let LETBINDINGS
          BODY (... ...)
        )))

      (for ([ii (in-range (+ nii2 2))])
        (let* ([i0jk4 (* ii jsize4)]
               [i1jk4  (+ (* 1 isize4) (* ii jsize4))]
               [i2jk4  (+ (* 2 isize4) (* ii jsize4))]
               [i3jk4  (+ (* 3 isize4) (* ii jsize4))]
               [i4jk4  (+ (* 4 isize4) (* ii jsize4))]
               [ijk1   (+ i (* j jsize1) (* k ksize1))]
               [ijk2   (+ (* i isize2) (* j jsize2) (* k ksize2))]
               [ijk21  (+ 1 ijk2)]
               [ijk22  (+ 2 ijk2)]
               [ijk23  (+ 3 ijk2)]
               [ijk24  (+ 4 ijk2)]
               [tmp1 (fr rho_i ijk1)]
               [tmp2 (sqr tmp1)]
               [tmp3 (* tmp1 tmp2)])

          (define-syntax-case (DIAG0 V A V1 V2 V3 V4 V5)
            (if (= (syntax->datum #'R) 4)
              #'(begin
                (printf "~a ~a ~a\n" A (+ A i0jk4) V1)
                (printf "~a ~a ~a\n" A (+ A i1jk4) V2)
                (printf "~a ~a ~a\n" A (+ A i2jk4) V3)
                (printf "~a ~a ~a\n" A (+ A i3jk4) V4)
                (printf "~a ~a ~a\n" A (+ A i4jk4) V5)
                (f! V (+ A i0jk4) V1)
                (f! V (+ A i1jk4) V2)
                (f! V (+ A i2jk4) V3)
                (f! V (+ A i3jk4) V4)
                (f! V (+ A i4jk4) V5))
              #'(begin
              (f! V (+ A i0jk4) V1)
              (f! V (+ A i1jk4) V2)
              (f! V (+ A i2jk4) V3)
              (f! V (+ A i3jk4) V4)
              (f! V (+ A i4jk4) V5))))

          (define-syntax-case (RETDET R2 N1 N2 N3 V1 V2 V3 BODY (... ...))
            (with-syntax ([LETBINDINGS (case (syntax->datum #'R2)
                                    [(1) #'([N1 V1] [N2 V2] [N3 V3])]
                                    [(2) #'([N1 V2] [N2 V1] [N3 V3])]
                                    [(3) #'([N1 V3] [N2 V2] [N3 V1])]
                                    [else (raise (format "invalid RETDET value ~a" #'R1))])])
            #'(let LETBINDINGS
              BODY (... ...)
            )))

          (define-syntax-case (ROTASN V A R1 V1 V2 V3 V4 V5)
            #'(RETDET R1 N2 N3 N4 V2 V3 V4 (DIAG0 V A V1 N2 N3 N4 V5)))

          (ROTASN fjac 0 R 0.0 1.0 0.0 0.0 0.0)
          (RETDET R A1 A2 A3 1 2 3
              (ROTASN fjac A1 R
              (- (fl* c2 (fr qs ijk1)) (fl* tmp2 (sqr (fr u (+ A1 ijk2)))))
              (/ (* (- 2.0 c2) (fr u (+ A1 ijk2))) (fr u ijk2))
              (- (* c2 tmp1 (fr u (+ A2 ijk2))))
              (- (* c2 tmp1 (fr u (+ A3 ijk2))))
              c2)
            (ROTASN fjac A2 R
              (- (* tmp2 (fr u (+ A2 ijk2)) (fr u (+ R ijk2))))
              (* tmp1 (fr u (+ A2 ijk2)))
              (* tmp1 (fr u (+ R ijk2)))
              0.0
              0.0)
            (ROTASN fjac A3 R
              (- (* tmp2 (fr u (+ A3 ijk2)) (fr u (+ R ijk2))))
              (* tmp1 (fr u (+ A3 ijk2)))
              0.0
              (* tmp1 (fr u (+ R ijk2)))
              0.0))
          (RETDET R A1 A2 A3 1 2 3 
            (ROTASN fjac 4 R
              (* tmp2 (fr u (+ R ijk2)) (- (* c2 2.0 (fr square ijk1))
                                           (* c1 (fr u (+ 4 ijk2)))))
              (- (* c1 tmp1 (fr u (+ 4 ijk2)))
                 (* c2 tmp2 (sqr (fr u (+ A1 ijk2))))
                 (fl* c2 (fr qs ijk1)))
              (- (* c2 tmp2 (fr u (+ A2 ijk2)) (fr u (+ R ijk2))))
              (- (* c2 tmp2 (fr u (+ A3 ijk2)) (fr u (+ R ijk2))))
              (* c1 tmp1 (fr u (+ R ijk2))))
              )

          (define-syntax-case (R43IT A a (... ...))
            (let ([aa (syntax->datum #'A)] 
                  [rr (syntax->datum #'R)])
              (if (= aa rr) 
                (let ([z #'(* con43 a (... ...))])
                  ;(printf "* ~a ~a ~a\n" aa rr z)
                  z)
                (let ([z #'(* a (... ...))])
                  ;(printf "  ~a ~a ~a\n" aa rr z)
                  z))))
                
           
          (define-syntax-case (NJACM A)
            #'
            (ROTASN njac A A
              (- (R43IT A c3c4 tmp2 (fr u (+ A ijk2))))
              (R43IT A c3c4 tmp1)
              0.0
              0.0
              0.0))


          (DIAG0 njac 0 0.0 0.0 0.0 0.0 0.0)
          (NJACM 1)
          (NJACM 2)
          (NJACM 3)
          (DIAG0 njac 4
            (+ (- (* (- (R43IT 1 c3c4) c1345) tmp3 (sqr (fr u ijk21))))
               (- (* (- (R43IT 2 c3c4) c1345) tmp3 (sqr (fr u ijk22))))
               (- (* (- (R43IT 3 c3c4) c1345) tmp3 (sqr (fr u ijk23))))
               (- (* c1345 tmp2 (fr u ijk24))))
            (* (- (R43IT 1 c3c4) c1345) tmp2 (fr u ijk21))
            (* (- (R43IT 2 c3c4) c1345) tmp2 (fr u ijk22))
            (* (- (R43IT 3 c3c4) c1345) tmp2 (fr u ijk23))
            (* c1345 tmp1))))
      
      (lhsinit lhs IISIZE isize4 jsize4 ksize4)
      
      (let ([tmp1 (fl* dt t_1)]
            [tmp2 (fl* dt t_2)]
            [dtd_1t_1*2+1 (+ 1.0 (fl* dtd_1t_1 2.0))]
            [dtd_2t_1*2+1 (+ 1.0 (fl* dtd_2t_1 2.0))]
            [dtd_3t_1*2+1 (+ 1.0 (fl* dtd_3t_1 2.0))]
            [dtd_4t_1*2+1 (+ 1.0 (fl* dtd_4t_1 2.0))]
            [dtd_5t_1*2+1 (+ 1.0 (fl* dtd_5t_1 2.0))])
      (for ([ii (in-range 1 IISIZE)])
        (let ([di (+ (* 0 jsize4) (* ii ksize4))]
              [si- (* (sub1 ii) jsize4)])
          (for ([m (in-range 5)])
            (let* ([mi4 (* m isize4)]
                   [dmi4 (+ di mi4)]
                   [smi4 (+ si- mi4)])
            (for ([n (in-range 5)])
              (f! lhs (+ n dmi4) (- (+ (fl* tmp2 (fr fjac (+ n smi4)))
                                       (fl* tmp1 (fr njac (+ n smi4))))))
            )))
          (f!- lhs (+ 0 (* 0 isize4) di) dtd_1t_1)
          (f!- lhs (+ 1 (* 1 isize4) di) dtd_2t_1)
          (f!- lhs (+ 2 (* 2 isize4) di) dtd_3t_1)
          (f!- lhs (+ 3 (* 3 isize4) di) dtd_4t_1)
          (f!- lhs (+ 4 (* 4 isize4) di) dtd_5t_1))

        (let ([di (+ (* 1 jsize4) (* ii ksize4))]
              [si= (* ii jsize4)])
          (for ([m (in-range 5)])
            (let* ([mi4 (* m isize4)]
                   [dmi4 (+ mi4 di)]
                   [smi4 (+ mi4 si=)])
            (for ([n (in-range 5)])
              (f! lhs (+ n dmi4) (* 2.0 tmp1 (fr njac (+ n smi4))))

            )))
          (f!+ lhs (+ 0 (* 0 isize4) di) dtd_1t_1*2+1)
          (f!+ lhs (+ 1 (* 1 isize4) di) dtd_2t_1*2+1)
          (f!+ lhs (+ 2 (* 2 isize4) di) dtd_3t_1*2+1)
          (f!+ lhs (+ 3 (* 3 isize4) di) dtd_4t_1*2+1)
          (f!+ lhs (+ 4 (* 4 isize4) di) dtd_5t_1*2+1))

        (let ([di (+ (* 2 jsize4) (* ii ksize4))]
              [si+ (* (add1 ii) jsize4)])
          (for ([m (in-range 5)])
            (let* ([mi4 (* m isize4)]
                   [dmi4 (+ mi4 di)]
                   [smi4 (+ mi4 si+)])
            (for ([n (in-range 5)])
              (f! lhs (+ n dmi4) (- (* tmp2 (fr fjac (+ n smi4)))
                                    (* tmp1 (fr njac (+ n smi4))))))))
          (f!- lhs (+ 0 (* 0 isize4) di) dtd_1t_1)
          (f!- lhs (+ 1 (* 1 isize4) di) dtd_2t_1)
          (f!- lhs (+ 2 (* 2 isize4) di) dtd_3t_1)
          (f!- lhs (+ 3 (* 3 isize4) di) dtd_4t_1)
          (f!- lhs (+ 4 (* 4 isize4) di) dtd_5t_1))))


      ; gaussian elimination 
      (define-syntax-rule (IDX2 I J K) (+ (* I isize2) (* J jsize2) (* K ksize2)))
      (define-syntax-rule (MINUSIT e) (- e 1))
      (define-syntax-rule (ZERO2 R) (REPDET R N1 N2 N3 i j k 0 0 0 (IDX2 N1 N2 N3)))
      (define-syntax-rule (REPIDX2 R V1 V2 V3 VV1 VV2 VV3) (REPDET R N1 N2 N3 V1 V2 V3 VV1 VV2 VV3 (IDX2 N1 N2 N3)))
    
      (define aa 0) 
      (define bb 1) 
      (define cc 2) 
      (binvcrhs lhs (* bb jsize4)
                lhs (* cc jsize4)
                rhs (ZERO2 R))

      (for ([ii (in-range 1 IISIZE)])
        (define-syntax-rule (MINUS2 R) (MODDET R MINUSIT N1 N2 N3 i j k (IDX2 N1 N2 N3)))
        (matvec_sub lhs (+ (* aa jsize4) (* ii ksize4))
                    rhs (MINUS2 R)
                    rhs (+ (* i isize2) (* j jsize2) (* k ksize2)))

        (matmul_sub lhs (+ (* aa jsize4) (* ii ksize4))
                    lhs (+ (* cc jsize4) (* (- ii 1) ksize4))
                    lhs (+ (* bb jsize4) (* ii ksize4)))

        (binvcrhs lhs (+ (* bb jsize4) (* ii ksize4))
                  lhs (+ (* cc jsize4) (* ii ksize4))
                  rhs (+ (* i isize2) (* j jsize2) (* k ksize2))))


      (matvec_sub lhs (+ (* aa jsize4) (* IISIZE ksize4))
                  rhs (REPIDX2 R i j k (- isize 1) (- jsize 1) (- ksize 1))
                  rhs (REPIDX2 R i j k isize jsize ksize))

      (matmul_sub lhs (+ (* aa jsize4) (* IISIZE ksize4))
                  lhs (+ (* cc jsize4) (* (- IISIZE 1) ksize4))
                  lhs (+ (* bb jsize4) (* IISIZE ksize4)))

      (binvrhs lhs (+ (* bb jsize4) (* IISIZE ksize4))
                rhs (REPIDX2 R i j k isize jsize ksize))

      (for ([ii (in-range (sub1 IISIZE) -1 -1)])
        (let ([BLOCK_SIZE 5])
        (for ([m (in-range BLOCK_SIZE)])
          (let ([midx (+ m (* i isize2) (* j jsize2) (* k ksize2))])
            (define-syntax-rule (PLUSIT e) (+ e 1))
            (define-syntax-rule (PLUS2 R) (MODDET R PLUSIT N1 N2 N3 i j k (IDX2 N1 N2 N3)))
            (f! rhs midx
              (for/fold ([x (fr rhs midx)]) ([n (in-range BLOCK_SIZE)])
                (let ([NIidx (+ m (* n isize4) (* cc jsize4) (* ii ksize4))])
                  (- x (* (fr lhs NIidx) (fr rhs (+ n (PLUS2 R))))))))))))))))
(define-values (i j k) (values 0 0 0))
(__solve x_solve 1 i j k k j i nz ny nx)
(__solve y_solve 2 i j k k i j nz nx ny)
(__solve z_solve 3 i j k j i k ny nx nz)

(define (checkSum arr nz2 ny2 nx2 isize1 jsize1 ksize1)
  (for*/fold ([csum 0.0]) ([k (in-range (add1 nz2))]
         [j (in-range (add1 ny2))]
         [i (in-range (add1 nx2))]
         [m (in-range 5)])
    (let* ([offset (+ m (* i isize1) (* j jsize1) (* k ksize1))]
           [arro (fr arr offset)])
      (+ csum (/ (sqr arro) (* (+ 2 nz2) (+ 2 ny2) (+ 2 nx2) 5))))))

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
