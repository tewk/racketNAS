#lang racket/base

(provide main)
  
(require "../bm-args.rkt") 
(require "../bm-results.rkt") 
(require "../rand-generator.rkt")
(require "../timer.rkt")
(require "../parallel-utils.rkt")
(require racket/future)
(require racket/list)
(require racket/match)
(require racket/math)
(require (for-syntax racket/base))

#;(require scheme/fixnum scheme/flonum)

(require (only-in scheme/flonum make-flvector make-shared-flvector)
         racket/require
         (filtered-in (lambda (name) (regexp-replace #rx"unsafe-" name "")) scheme/unsafe/ops))

(define DOPLACES #t)
(define mkflvec (if DOPLACES make-shared-flvector make-flvector))

(define-syntax (pfor stx)
  (syntax-case stx (in-range)
    [(_ cg serial nproc ([loopvar (in-range imax)]) ([letname letval] ... ) body ...)
       (begin
         (unless (identifier? #'lv) (raise-syntax-error 'parray "expected an identifier" stx #'lv))
         #'(if (not serial)
           (if (not DOPLACES)
             (let ([per-future (ceiling (/ imax nproc))])
               (for ([f (in-list
                         (for/list ([i (in-range 0 nproc)])
                           (let ([letname (vector-ref letval i)]...)
                             (future
                              (λ ()
                                  (for ([loopvar (in-range (fx* i per-future) (min imax (fx* (fx+ i 1) per-future)))])
                                    body ...))))))])
                 (touch f)))
             (begin
               (CG-B cg)
               (for ([loopvar (p-range cg (in-range imax))])
                 body ...)
               (CG-B cg)))
           (for ([loopvar (in-range imax)]) body ...)))]))

(define-syntax (timer stx)
  (syntax-case stx ()
    [(_ id body ...)
     #'(begin body ...)]))

(define (->inexact x)
  (if (exact? x)
    (exact->inexact x)
    x))
;Constants 
(define REAL 0)
(define IMAG 1)  

(define (get-class-size CLASS)
  (case CLASS 
    [(#\S) (values 64 64 64 6)]
    [(#\W) (values 128 128 32 6)] 
    [(#\A) (values 256 256 128 6)]
    [(#\B) (values 512 256 256 20)]
    [(#\C) (values 512 512 512 20)]
    [else (error "Unknown class")]))

(define (main . argv) 
  (let ([args (parse-cmd-line-args argv "Fourier Transform")]) 
    (run-benchmark args)))

(define (run-benchmark args) 
  (let ([bmname "FT"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])

  (let-values ([(nx ny nz niter-default) (get-class-size CLASS)])
    (let ([maxdim (max nx (max ny nz))]
          [checksum (make-vector (* 2 niter-default) 0.0)]) ;isize2 = 2 

  (print-banner "Fourier Transform" args) 
  (printf "Size = ~a X ~a X ~a niter = ~a~n" nx ny nz niter-default) 

  (let ([exp1 (mkflvec (* 2 nx) 0.0)] 
        [exp2 (mkflvec (* 2 ny) 0.0)] 
        [exp3 (mkflvec (* 2 nz) 0.0)]
        [xtr  (mkflvec (* (* (* 2 (+ ny 1)) nx) nz) 0.0)]
        [xnt  (mkflvec (* (* (* 2 (+ ny 1)) nz) nx) 0.0)])

    (CGspawn (if (or serial (not DOPLACES)) 0 num-threads)
             fft-body nx ny nz exp1 exp2 exp3 xnt xtr maxdim niter-default checksum serial num-threads))

  (timer 14
    (let ([verified (verify bmname CLASS niter-default checksum)])
      (if verified 
          (printf "Verification succeeded~n") 
          (printf "Verification failed~n"))
      (timer-stop 1) 
      (let* ([time (/ (read-timer 1) 1000)]
             [results (new-BMResults bmname CLASS nx ny nz niter-default time 
                                     (get-mflops time nx ny nz) 
                                     "floating point" 
                                     verified 
                                     serial 
                                     num-threads 
                                     0)]) 
        (print-results results))))))))

;;ilog2 : int -> int
(define (ilog2 n) 
  (if (= n 1) 
      0
      (let loop ([lg 1]
                 [nn 2])
        (if (< nn n)
          (loop (+ lg 1) (* 2 nn))
          lg))))

;;comp-exp : int vector-of-double -> void
(define (comp-exp n exponent) 
  (let ([nu n]
        [m (ilog2 n)] 
        [eps 1.0E-16] 
        [ku 1] 
        [ln 1]) 
    (flvector-set! exponent 0 (->inexact m))
    (for/fold ([ku 1] [ln 1]) ([j (in-range 1 (+ m 1))]) 
      (let ([t (/ pi ln)]) 
        (for ([i (in-range ln)])
          (let ([ti (* i t)] 
                [idx (* (+ i ku) 2)]) 
            (flvector-set! exponent (+ REAL idx) (->inexact (cos ti)))
            (flvector-set! exponent (+ IMAG idx) (->inexact (sin ti)))
            (when (< (abs (flvector-ref exponent (+ REAL idx))) eps) (flvector-set! exponent (+ REAL idx) 0.0))
            (when (< (abs (flvector-ref exponent (+ IMAG idx))) eps) (flvector-set! exponent (+ IMAG idx) 0.0))
)))
      (values (+ ku ln) (* 2 ln)))))

(define-syntax-rule (helper1 i it ap n12 n22 n32 nx ny nz isize4 jsize4 ksize4 isize3 jsize3 ksize3 xnt xtr)
          (let* ([ii (fx- i (fx* (fxquotient i n12) nx))] 
                 [ii2 (fx* ii ii)]) 
            (for ([k (in-range nz)]) 
              (let* ([kk (fx- k (fx* (fxquotient k n32) nz))] 
                     [ik2 (fx+ ii2 (fx* kk kk))]) 
                (for ([j (in-range ny)]) 
                  (let* ([jj (fx- j (fx* (fxquotient j n22) ny))]
                        [idx1 (fx+ (fx+ (fx* j isize4) (fx* k jsize4)) (fx* i ksize4))]
                        [idx2 (fx+ (fx+ (fx* j isize3) (fx* i jsize3)) (fx* k ksize3))]
                        [flval (flexp (fl* ap (exact->inexact (fx* (fx+ (fx* jj jj) ik2) (fx+ it 1)))))])
                    (flvector-set! xnt (fx+ REAL idx1) (fl* (flvector-ref xtr (fx+ REAL idx2)) flval)) 
                    (flvector-set! xnt (fx+ IMAG idx1) (fl* (flvector-ref xtr (fx+ IMAG idx2)) flval))))))))

(define-syntax-rule (helper2 k sign log ic jc isize3 jsize3 ksize3 isize1 jsize1 x plane exp scr)
  (begin
    (define-syntax-rule (x->plane x x-idx plane plane-idx RI)
      (flvector-set! plane (fx+ RI plane-idx) (flvector-ref x (fx+ RI x-idx))))
    (define-syntax-rule (plane->x x x-idx plane plane-idx RI)
      (flvector-set! x (fx+ RI x-idx) (flvector-ref plane (fx+ RI plane-idx))))
    (define-syntax-rule (helper3 k ic jc isize3 jsize3 ksize3 isize1 jsize1 x plane copy)
      (for ([i (in-range ic)]) 
        (for ([j (in-range jc)]) 
           (let ([plane-idx (fx+ (fx* j isize1) (fx* i jsize1))]
                 [x-idx     (fx+ (fx+ (fx* i isize3) (fx* j jsize3)) (fx* k ksize3))])
             (copy x x-idx plane plane-idx REAL)
             (copy x x-idx plane plane-idx IMAG)))))


    (helper3 k ic jc isize3 jsize3 ksize3 isize1 jsize1 x plane x->plane)
    (swarztrauber sign log jc ic plane 0 jc exp scr) 
    (helper3 k ic jc isize3 jsize3 ksize3 isize1 jsize1 x plane plane->x)))

(define-syntax-rule (fft-xyz cg sign x scr plane exp1 exp2 exp3 n1 n2 n3 log1 log2 log3 isize3 jsize3 ksize3 isize1 jsize11 jsize12 serial np) 
  (begin
    (pfor cg serial np ([k (in-range n3)]) 
      ([scr scr])
      (swarztrauber sign log2 n1 n2 x (fx* k ksize3) n1 exp2 scr))

    (pfor cg serial np ([k (in-range n3)])
      ([scr scr]
       [plane plane])
      (helper2 k sign log1 n1 n2 isize3 jsize3 ksize3 isize1 jsize11 x plane exp1 scr))

    (pfor cg serial np ([k (in-range n2)])
      ([scr scr]
       [plane plane])
      (helper2 k sign log3 n3 n1 ksize3 isize3 jsize3 isize1 jsize12 x plane exp3 scr))))

(define (fft-body cg nx ny nz exp1 exp2 exp3 xnt xtr maxdim niter-default checksum serial np)
  (define alpha .000001)
  (define-syntax-rule (palloc serial nproc body ...)
    (if (and (not serial) (not DOPLACES))
      (apply vector (for/list ([i (in-range nproc)]) body ...))
      (begin body ...)))

  (let ([isize1 2]
        [isize2 2]
        [isize3 2]
        [isize4 2]
        [jsize11 (fx* 2 (fx+ nx 1))]
        [jsize12 (fx* 2 (fx+ ny 1))]
        [jsize13 (fx* 2 (fx+ nz 1))]
        [jsize3 (* 2 (+ ny 1))]
        [jsize4 (* 2 (+ ny 1))]
        [ksize3 (* (* 2 (+ ny 1)) nx)]
        [ksize4 (* (* 2 (+ ny 1)) nz)]
        [ap (fl* (fl* -4.0 alpha) (expt pi 2))] 
        [n12 (quotient nx 2)] 
        [n22 (quotient ny 2)] 
        [n32 (quotient nz 2)]
        [logx (ilog2 nx)]
        [logy (ilog2 ny)]
        [logz (ilog2 nz)]
        [scr (palloc serial np (mkflvec (* (* 2 (+ maxdim 1)) maxdim) 0.0))]
        [plane (palloc serial np (mkflvec (* (* 2 (+ maxdim 1)) maxdim) 0.0))])

  (timer 2
    (CG-n0-only cg
      (initial-conditions xtr ny nx nz maxdim) 
      (comp-exp nx exp1) 
      (comp-exp ny exp2) 
      (comp-exp nz exp3))
    (fft-xyz cg 1.0 xtr scr plane exp2 exp1 exp3 ny nx nz logy logx logz isize3 jsize3 ksize3 isize1 jsize11 jsize12 serial np))

    (timer-start 1) 
    (timer 12
      (CG-n0-only cg
        (initial-conditions xtr ny nx nz maxdim))
    (fft-xyz cg 1.0 xtr scr plane exp2 exp1 exp3 ny nx nz logy logx logz isize3 jsize3 ksize3 isize1 jsize11 jsize12 serial np))
  (timer 15
    (for ([it (in-range niter-default)]) 
      (timer 11
      (pfor cg serial np ([i (in-range nx)]) ()
        (helper1 i it ap n12 n22 n32 nx ny nz isize4 jsize4 ksize4 isize3 jsize3 ksize3 xnt xtr)))
      (timer 15
      (fft-xyz cg -1.0 xnt scr plane exp2 exp3 exp1 ny nz nx logy logz logx isize4 jsize4 ksize4 isize1 jsize13 jsize12 serial np))
      (timer 10
        (CG-n0-only cg
            (calculate-checksum checksum (fx+ REAL (fx* it isize2)) it xnt ny nz nx)))))))

(define-syntax-rule (swarztrauber-worker is n exponent block-start block-end l j x jsize1 scr lk li kbt-offset1 kbt-offset1o kbt-offset2)
  (let* ([n1 (fxquotient n 2)]
;         [lk lk]
;         [li li]
         [lj (fx* 2 lk)]
         [ku li]) 
    (for ([i (in-range li)]) 
      (let* ([i11 (fx* i lk)]
             [i12 (fx+ i11 n1)]
             [i21 (fx* i lj)]
             [i22 (fx+ i21 lk)]
             [kbt-idx (fx* (fx+ ku i) 2)]
             [u1-REAL (flvector-ref exponent (fx+ REAL kbt-idx))]
             [u1-IMAG (fl* is (flvector-ref exponent (fx+ IMAG kbt-idx)))])
        (for ([k (in-range lk)]) 
          (for ([j (in-range block-start (fx+ block-end 1))]) 
            (let* (;[kbt-offset1 kbt-offset1]
;                   [kbt-offset1o kbt-offset1o]
;                   [kbt-offset2 kbt-offset2]
                   [kbt-REALoff (fx+ REAL kbt-offset1)]
                   [kbt-IMAGoff (fx+ IMAG kbt-offset1)]
                   [i11-k-jsize1 (fx* (fx+ i11 k) jsize1)]
                   [i12-k-jsize1 (fx* (fx+ i12 k) jsize1)])
              (let ([x11-REAL (flvector-ref x (fx+ i11-k-jsize1 kbt-REALoff))]
                    [x11-IMAG (flvector-ref x (fx+ i11-k-jsize1 kbt-IMAGoff))]
                    [x21-REAL (flvector-ref x (fx+ i12-k-jsize1 kbt-REALoff))]
                    [x21-IMAG (flvector-ref x (fx+ i12-k-jsize1 kbt-IMAGoff))])
                (let ([j-isize-i21-k-jsize (fx+ kbt-offset1o (fx* (fx+ i21 k) jsize1))])
                   (flvector-set! scr (fx+ REAL j-isize-i21-k-jsize) (fl+ x11-REAL x21-REAL))
                   (flvector-set! scr (fx+ IMAG j-isize-i21-k-jsize) (fl+ x11-IMAG x21-IMAG)))
                (let ([j-2-i22-k-jsize1 (fx+ kbt-offset2 (fx* (fx+ i22 k) jsize1))]
                      [x11-x21-REAL (fl- x11-REAL x21-REAL)]
                      [x11-x21-IMAG (fl- x11-IMAG x21-IMAG)])
                   (flvector-set! scr (fx+ REAL j-2-i22-k-jsize1) (fl- (fl* u1-REAL x11-x21-REAL) (fl* u1-IMAG x11-x21-IMAG)))
                   (flvector-set! scr (fx+ IMAG j-2-i22-k-jsize1) (fl+ (fl* u1-IMAG x11-x21-REAL) (fl* u1-REAL x11-x21-IMAG))))))))))))

;;swarztrauber : int int int int vector-of-double int int vector-of-double vector-of-double -> void
(define (swarztrauber is m len n x xoffset xd1 exponent scr)
  (define fftblock-default (fx* 4 4096)) ;Size of L1 cache on SGI 02K
  (let ([isize1 2] 
        [jsize1 (fx* 2 (fx+ xd1 1))])
    (timer 4
      ;Perform one variant of the Stockham FFT 
      (let ([fftblock (fxmax 8 (fxquotient fftblock-default n))])
      (for ([block-start (in-range 0 len fftblock)])
        (let ([block-end (fxmin (fx- (fx+ block-start fftblock) 1) (fx- len 1))])
          (for ([l (in-range 1 (fx+ m 1) 2)]) 

            (swarztrauber-worker is n exponent block-start block-end l j x jsize1 scr
              (fxlshift 1 (fx- l 1))
              (fxlshift 1 (fx- m l)) 
              (fx+ (fx* j 2) xoffset) 
              (fx* j isize1)
              (fx* j 2))

            (if (fx= l m) 
              (for ([k (in-range n)]) 
                (for ([j (in-range block-start (fx+ block-end 1))]) 
                  (let ([dest-off (fx+ (fx+ (fx* j 2) (fx* k jsize1)) xoffset)]
                        [src-off  (fx+ (fx* j isize1) (fx* k jsize1))])
                    (flvector-set! x (fx+ REAL dest-off) (flvector-ref scr (fx+ REAL src-off)))
                    (flvector-set! x (fx+ IMAG dest-off) (flvector-ref scr (fx+ IMAG src-off))))))

              (swarztrauber-worker is n exponent block-start block-end l j scr jsize1 x
                (fxlshift 1 l) 
                (fxlshift 1 (fx- (fx- m l) 1)) 
                (fx* j isize1) 
                (fx+ (fx* j 2) xoffset) 
                (fx+ (fx* j 2) xoffset)))
  )))))))

(define (get-mflops total-time nx ny nz) 0)

(define (initial-conditions u0 d1 d2 d3 maxdim) 
  (let* ([tmp (make-flvector (* 2 maxdim) 0.0)] 
         [ran-starts (make-vector maxdim 0.0)] 
         [seed 314159265.0] 
         [a (expt 5.0 13)] 
         [start seed]
         [isize3 2] 
         [jsize3 (* isize3 (+ d1 1))] 
         [ksize3 (* jsize3 d2)])
    ;Jump to the starting element for our first plane 
    (random-init seed) 
    (let ([an (ipow46 a 0)]) 
      (randlc/a seed an) 
      (set! an (ipow46 a (* (* 2 d1) d2))) 
      ;Go through by z planes filling in one square at a time 
      (vector-set! ran-starts 0 start) 
      (for/fold ([start start]) ([k (in-range 1 d3)])
        (let ([newseed (randlc/a start an)])
          (vector-set! ran-starts k newseed)
          newseed))
 
      (for ([k (in-range d3)]) 
        (let ([x0 (vector-ref ran-starts k)]) 
          (for ([j (in-range d1)]) 
            
            (set! x0 (vranlc (* 2 d2) x0 a tmp 0)) 
            (for ([i (in-range d2)])
              (let ([dest-idx (+ (+ (* j isize3) (* i jsize3)) (* k ksize3))]
                    [iidx (* i 2)])
                (flvector-set! u0 (+ REAL dest-idx) (flvector-ref tmp (+ REAL iidx)))
                (flvector-set! u0 (+ IMAG dest-idx) (flvector-ref tmp (+ IMAG iidx)))))))))))

(define (calculate-checksum csum csmffst iterN u d1 d2 d3) 
  (let* (
         [isize3 2] 
         [jsize3 (* isize3 (+ d1 1))] 
         [ksize3 (* jsize3 d2)]
         [d123 (* d1 (* d2 d3))])
    (let-values ([(csumr csumi)
      (for/fold ([csumr 0.0]
                 [csumi 0.0])
                ([i (in-range 1 1025)]) 
        (let* ([ii (fxmodulo (fx* 1 i) d3)]
               [ji (fxmodulo (fx* 3 i) d1)] 
               [ki (fxmodulo (fx* 5 i) d2)]
               [idx (fx+ (fx+ (fx* ji isize3) (fx* ki jsize3)) (fx* ii ksize3))])
          (values (fl+ csumr (flvector-ref u (fx+ REAL idx)))
                  (fl+ csumi (flvector-ref u (fx+ IMAG idx))))))])
    (vector-set! csum (+ REAL csmffst) (/ csumr d123)) 
    (vector-set! csum (+ IMAG csmffst) (/ csumi d123)))))

;;verify : int int int int int vector-of-double -> boolean
(define (verify bmname CLASS niter-default cksum) 
  (let ([cexpd
    (case CLASS
       ;Class S reference values
    [(#\S) 
         #(
         554.6087004964 484.5363331978;IMAG 0
         554.6385409189 486.5304269511;IMAG 1  
         554.6148406171 488.3910722336;IMAG 2
         554.5423607415 490.1273169046;IMAG 3
         554.4255039624 491.7475857993;IMAG 4
         554.2683411902 493.2597244941;IMAG 5
        )]

    [(#\W)
          #( 
            567.3612178944 529.3246849175;IMAG 0
            563.1436885271 528.2149986629;IMAG 1
            559.4024089970 527.0996558037;IMAG 2
            556.0698047020 526.0027904925;IMAG 3
            553.0898991250 524.9400845633;IMAG 4
            550.4159734538 523.9212247086;IMAG 5
          )]

    [(#\A) 
          #( 
            504.6735008193 511.4047905510;IMAG 0
            505.9412319734 509.8809666433;IMAG 1
            506.9376896287 509.8144042213;IMAG 2
            507.7892868474 510.1336130759;IMAG 3
            508.5233095391 510.4914655194;IMAG 4
            509.1487099959 510.7917842803;IMAG 5
          )]
    [(#\B) 
     ;Class B reference values 
          #(
            517.7643571579 507.7803458597;IMAG 0
            515.4521291263 508.8249431599;IMAG 1
            514.6409228649 509.6208912659;IMAG 2
            514.2378756213 510.1023387619;IMAG 3
            513.9626667737 510.3976610617;IMAG 4
            513.7423460082 510.5948019802;IMAG 5
            513.5547056878 510.7404165783;IMAG 6
            513.3910925466 510.8576573661;IMAG 7
            513.2470705390 510.9577278523;IMAG 8
            513.1197729984 511.0460304483;IMAG 9
            513.0070319283 511.1252433800;IMAG 10
            512.9070537032 511.1968077718;IMAG 11
            512.8182883502 511.2616233064;IMAG 12
            512.7393733383 511.3203605551;IMAG 13
            512.6691062020 511.3735928093;IMAG 14
            512.6064276004 511.4218460548;IMAG 15
            512.5504076570 511.4656139760;IMAG 16
            512.5002331720 511.5053595966;IMAG 17
            512.4551951846 511.5415130407;IMAG 18
            512.4146770029 511.5744692211;IMAG 19
          )]
    [(#\C) 
     ;Class C reference values 
          #(
            519.5078707457  514.9019699238;IMAG 0   
            515.5422171134  512.7578201997;IMAG 2
            514.4678022222  512.2251847514;IMAG 3
            514.0150594328  512.1090289018;IMAG 4
            513.7550426810  512.1143685824;IMAG 5
            513.5811056728  512.1496764568;IMAG 6
            513.4569343165  512.1870921893;IMAG 7
            513.3651975661  512.2193250322;IMAG 8
            513.2955192805  512.2454735794;IMAG 9
            513.2410471738  512.2663649603;IMAG 10
            513.1971141679  512.2830879827;IMAG 11
            513.1605205716  512.2965869718;IMAG 12
            513.1290734194  512.3075927445;IMAG 13
            513.1012720314  512.3166486553;IMAG 14
            513.0760908195  512.3241541685;IMAG 15
            513.0528295923  512.3304037599;IMAG 16
            513.0310107773  512.3356167976;IMAG 17
            513.0103090133  512.3399592211;IMAG 18
            512.9905029333  512.3435588985;IMAG 19
            512.9714421109  512.3465164008;IMAG 20
          )])])
  (let* ([epsilon 1.0E-12]
        [verified 
    (if (<= niter-default 0) 
        -1
        (for/fold ([verified #t]) ([it (in-range niter-default)]) 
          (let* ([rit2 (+ REAL (* it 2))]
                 [iit2 (+ IMAG (* it 2))]
                 [csumr (/ (- (vector-ref cksum rit2) (vector-ref cexpd rit2)) (vector-ref cexpd rit2))] 
                 [csumi (/ (- (vector-ref cksum iit2) (vector-ref cexpd iit2)) (vector-ref cexpd iit2))])
            (if (or (<= (abs csumr) epsilon)
                    (<= (abs csumi) epsilon)) 
                (and verified #t)
                (begin 
                  (printf "Verification failure: ~n") 
                  (printf "epsilon: ~a~n" epsilon)
                  (printf "csumr: ~a~n" csumr) 
                  (printf "csumi: ~a~n" csumi)
                  (and verified #f))))))])
    (print-verification-status CLASS 
      (case verified
        [(-1) -1]
        [(#t)  1]
        [(#f)  0])
       bmname)
    verified)))
