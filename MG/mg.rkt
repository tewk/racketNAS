#lang racket/base

(provide main)
  
(require "../bm-args.rkt") 
(require "../bm-results.rkt") 
(require "../rand-generator.rkt")
(require "../timer.rkt")
(require "../parallel-utils.rkt")
(require "../macros.rkt")
(require "../debug.rkt")
(require racket/match)
(require racket/math)
(require (for-syntax racket/base
                     racket/list))

(require (only-in scheme/flonum make-flvector make-shared-flvector shared-flvector flvector-length))
;(require (only-in scheme/flonum flvector-set! flvector-ref))
;(define vr vector-ref)
;(define vs! vector-set!)
;(define flvs! flvector-set!)
;(define flvr flvector-ref)
;(require (only-in scheme/fixnum fx+ fx- fx*))

#|
(require scheme/fixnum scheme/flonum)

(require (only-in scheme/flonum make-flvector make-shared-flvector)
         scheme/require (for-syntax scheme/base)
   (filtered-in
    (lambda (name) (regexp-replace #rx"unsafe-" name ""))
    scheme/unsafe/ops))
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
                    [unsafe-fl+ fl+op]
                    [unsafe-fl- fl-op]
                    [unsafe-fl* fl*op]
                    [unsafe-fl/ fl/]
                    [unsafe-fx+ fx+op]
                    [unsafe-fx- fx-]
                    [unsafe-fx* fx*op]
                    [unsafe-fx= fx=]
))

(define-syntax (*it stx)
  (syntax-case stx ()
    [(_ id)
    (with-syntax ([NAME (datum->syntax #'id (string->symbol (string-append (symbol->string (syntax->datum #'id)) "*")))])
      #'(define-syntax (NAME stx)
        (syntax-case stx ()
          [(_ a b) #'(id a b)]
          [(_ a b c) #'(id a (id b c))]
          [(_ a (... ...)) (let-values ([(b c) (split-at #'(a (... ...)) (quotient (length #'(a (... ...)))))])
            (with-syntax ([(d (... ...)) b]
                          [(e (... ...)) c])
            #'(id (NAME d (... ...)) (NAME e (... ...)))))])))]))
(*it fl+)
(*it fl-)

(define-syntax (fx** stx)
  (define (halves lst)
    (values (split-at lst (quotient (length lst) 2))))
  (syntax-case stx ()
    [(_ a b) #'(fx* a b)]
    [(_ a b c) #'(fx* a (fx* b c))]
    [(_ a ...) (let-values ([(b c) (halves #'(a ...))])
      (with-syntax ([(d ...) b]
                    [(e ...) c])
      #'(fx* (fx** d ...) (fx** e ...))))]))
    
 
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

(define make-fxvector make-vector)

(define (run-benchmark args) 
  (define maxlevel 11)
  (let ([bmname "MG"]
        [CLASS (BMArgs-class args)]
        [num-threads (BMArgs-num-threads args)]
        [serial (BMArgs-serial args)])
  
  (let-values ([(nx_default ny_default nz_default nit_default lm lt_default ndim1 ndim2 ndim3) (get-class-size CLASS)])
    (let* (
          [nx (make-vector maxlevel 0)]
          [ny (make-vector maxlevel 0)]
          [nz (make-vector maxlevel 0)]
          [ir (make-vector maxlevel 0)]
          [m1 (make-vector maxlevel 0)]
          [m2 (make-vector maxlevel 0)]
          [m3 (make-vector maxlevel 0)]
          [lt lt_default]
          [lt1 (fx-- lt)] 
          [nit nit_default]
          [nm (+ (arithmetic-shift 1 lm) 2)]
          [nv (* (+ 2 (arithmetic-shift 1 ndim1)) 
                 (+ 2 (arithmetic-shift 1 ndim2)) 
                 (+ 2 (arithmetic-shift 1 ndim3)))]
          [nr (floor (/ (* 8 (+ nv (* nm nm) (* 5 nm) (* 7 lm))) 7))]
          [nm2 (* 2 nm nm)]
          [r (make-shared-flvector nr 0.0)]
          [v (make-shared-flvector nv 0.0)]
          [u (make-shared-flvector nr 0.0)]
          [a (shared-flvector (/ -8.0 3.0) 0.0 (/ 1.0 6.0) (/ 1.0 12.0))]
          [c (case CLASS 
               [(#\A #\S #\W) 
                  (shared-flvector (/ -3.0 8.0) (/ 1.0 32.0) (/ -1.0 64.0) 0.0)]
               [else
                  (shared-flvector (/ -3.0 17.0) (/ 1.0 33.0) (/ -1.0 61.0) 0.0)])])
        (vs! nx lt1 nx_default)
        (vs! ny lt1 ny_default)
        (vs! nz lt1 nz_default)
        (let-values ([(is1 is2 is3 ie1 ie2 ie3 n1 n2 n3) (setup lt ir nx ny nz m1 m2 m3 )])
      
      (print-banner "MG" args) 
      (get-input-pars maxlevel)
      (printf " Size:  ~ax~ax~a Iterations:   ~a~n" (vr nx lt1) (vr ny lt1) (vr nz lt1)  nit) 

      
      (CGspawn (if serial 0 num-threads) mg-body c a u v r is1 is2 is3 ie1 ie2 ie3 n1 n2 n3 ir m1 m2 m3 nm lt nit nx ny)

      (printf " Initialization time: ~a seconds\n" 0)
      (let* ([verified (verify CLASS 
              (norm2u3 r n1 n2 n3 void (vr nx lt1) (vr ny lt1) (vr nz lt1)))])
        (print-verification-status CLASS verified bmname)
        (let* ([time (/ (read-timer 1) 1000)]
               [results (new-BMResults bmname CLASS (vr nx lt1) (vr ny lt1) (vr nz lt1)  nit time 
                                       (get-mflops time nit n1 n2 n3)
                                       "floating point" 
                                       (if verified 1 0)
                                       serial 
                                       num-threads 
                                       -1)]) 
            (print-results results))))))))

(define (mg-body cg c a u v r is1 is2 is3 ie1 ie2 ie3 n1 n2 n3 ir m1 m2 m3 nm lt nit nx ny)
;;;//--------------------------------------------------------------------
;;;//    One iteration for startup
;;;//--------------------------------------------------------------------
  (define u1 (make-flvector (fx++ nm) 0.0))
  (define u2 (make-flvector (fx++ nm) 0.0))
  (define m 535)
  (define z1 (make-flvector m 0.0))
  (define z2 (make-flvector m 0.0))
  (define z3 (make-flvector m 0.0))


  (CG-n0-only cg
    (zero3 u 0 n1 n2 n3)
    (zran3 cg v n1 n2 n3 (vr nx (fx-- lt)) (vr ny (fx-- lt)) is1 is2 is3 ie1 ie2 ie3))

  (resid cg a u v r 0 n1 n2 n3 nm u1 u2)

  (mg3P cg c a u v r n1 n2 n3 ir m1 m2 m3 nm lt u1 u2 z1 z2 z3)
  (resid cg a u v r 0 n1 n2 n3 nm u1 u2)

;;;//--------------------------------------------------------------------
;;;//    Main Loop
;;;//--------------------------------------------------------------------

  (CG-n0-only cg
    (zero3 u 0 n1 n2 n3)
    (zran3 cg v n1 n2 n3 (vr nx (fx-- lt)) (vr ny (fx-- lt)) is1 is2 is3 ie1 ie2 ie3)
    (timer-start 1))

  (resid cg a u v r 0 n1 n2 n3 nm u1 u2)

  (for ([it (in-range 1 (fx++ nit))])
    (mg3P cg c a u v r n1 n2 n3 ir m1 m2 m3 nm lt u1 u2 z1 z2 z3)
    (resid cg a u v r 0 n1 n2 n3 nm u1 u2))

  (CG-n0-only cg
    (timer-stop 1)))

  
  

(define-syntax-rule (vidx3 i1 i2 i3 n1 n2) (fx+ i1 (fx* n1 (fx+ i2 (fx* n2 i3)))))
(define-syntax-rule (vr3 v i1 i2 i3 n1 n2) (flvr v (vidx3 i1 i2 i3 n1 n2)))
(define-syntax-rule (vidx off i1 i2 i3 n1 n2) (fx+ off (vidx3 i1 i2 i3 n1 n2)))
(define-syntax-rule (vro3 v off i1 i2 i3 n1 n2) (flvr v (fx+ off (vidx3 i1 i2 i3 n1 n2))))

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
        (printf " The correct L2 Norm is ~a\n" verify-value)
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
          (printf " lt=~a Maximum allowable=~a\n" lt maxlevel)
          (exit 0))
        (values nit lt lnx lnz)]
      [else 
        (printf " Error reading from file mg.input\n")
        (exit 0)])
    (printf " No input file mg.input, Using compiled defaults\n")))

(define (setup lt ir nx ny nz m1 m2 m3 )
  (define lt1 (fx-- lt))

  (for ([k (in-range (- lt 2) -1 -1)])
    (let ([k1 (fx++ k)])
        (vs! nx k (/ (vr nx k1) 2))
        (vs! ny k (/ (vr ny k1) 2))
        (vs! nz k (/ (vr nz k1) 2))))
    
  (for ([k (in-range (- lt 1) -1 -1)])
    (let ([k1 (fx++ k)])
        (vs! m1 k (+ 2 (vr nx k)))
        (vs! m2 k (+ 2 (vr ny k))) 
        (vs! m3 k (+ 2 (vr nz k)))))

  (for ([j (in-range (- lt 2) -1 -1)])
    (let ([j1 (fx++ j)])
    (vs! ir j (+ (vr ir j1) (* (vr m1 j1) (vr m2 j1) (vr m3 j1))))))

  (values 2 2 2 (fx++ (vr nx lt1)) (fx++ (vr ny lt1)) (fx++ (vr nz lt1))
                (+ 2 (vr nx lt1)) (+ 2 (vr ny lt1)) (+ 2 (vr nz lt1))))

(define (zero3 z off n1 n2 n3)
  (for* ([i3 (in-range n3)]
         [i2 (in-range n2)]
         [i1 (in-range n1)])
    (flvs! z (fx+ off i1 (fx* n1 (fx+ i2 (fx* n2 i3)))) 0.0)))

(define (zran3 cg z n1 n2 n3 nx ny is1 is2 is3 ie1 ie2 ie3)
  (define mm 10)
  (define j1 (make-fxvector (* 2 mm) 0))
  (define j2 (make-fxvector (* 2 mm) 0))
  (define j3 (make-fxvector (* 2 mm) 0))
  (define jg (make-fxvector (* 4 2 mm) 0))
  (define ten (make-flvector (* 2 mm) 0.0))
  (zero3 z 0 n1 n2 n3)

  (let* ([i (vidx3 (- is1 2) (- is2 2) (- is3 2) nx ny)]
        [d1 (fx++ (- ie1 is1))]
        [e1 (+ 2 (- ie1 is1))]
        [e2 (+ 2 (- ie2 is2))]
        [e3 (+ 2 (- ie3 is3))]
        [a (expt 5.0 13)]
        [rng (new-rand)]
        [a1 (power/r rng a nx)]
        [a2 (power/r rng a (* nx ny))]
        [ai (power/r rng a i)])
    (for/fold ([x0 (randlc/a 314159265.0 ai)])  ([i3 (in-range 2 (fx++ e3))])
      (for/fold ([x1 x0]) ([i2 (in-range 2 (fx++ e2))])
        (vranlc d1 x1 a z (vidx3 1 (fx-- i2) (fx-- i3) n1 n2))
        (randlc/a x1 a1))
      (randlc/a x0 a2)))

  (for ([i (in-range mm)])
    (let ([imm (+ i mm)])
      (flvs! ten imm 0.0)
      (vs! j1 imm 0)
      (vs! j2 imm 0)
      (vs! j3 imm 0)
      (flvs! ten i 1.0)
      (vs! j1 i 0)
      (vs! j2 i 0)
      (vs! j3 i 0)))
  
  (for* ([i3 (in-range 1 (fx-- n3))]
         [i2 (in-range 1 (fx-- n2))]
         [i1 (in-range 1 (fx-- n1))])
    (let ([zv (vr3 z i1 i2 i3 n1 n2)])
      (when (zv . > . (flvr ten mm))
        (flvs! ten mm zv)
        (vs! j1 mm i1)
        (vs! j2 mm i2)
        (vs! j3 mm i3)
        (bubble ten j1 j2 j3 mm 1))
      (when (zv . < . (flvr ten 0))
        (flvs! ten 0 zv)
        (vs! j1 0 i1)
        (vs! j2 0 i2)
        (vs! j3 0 i3)
        (bubble ten j1 j2 j3 mm 0))))

  
  (define-syntax-rule (zidx i1 i2 i3) (fx+ i1 (fx* n1 (fx+ i2 (fx* n2 i3)))))
  (define-syntax-rule (vrzj ii) (flvr z (zidx (vr j1 ii) (vr j2 ii) (vr j3 ii))))
  (let-values ([(i0 i1)
    (for/fold ([i0 mm]
               [i1 mm]) ([i (in-range (fx-- mm) -1 -1 )])
      (let* ([ii (+ (fx-- i1) mm)]
             [im (+ i mm)]
             [im4 (* 4 im)]
             [best (vrzj ii)])
        (vs! jg im4 0)
        (vs! jg (+ im4 1) (+ (- is1 2) (vr j1 ii)))
        (vs! jg (+ im4 2) (+ (- is2 2) (vr j2 ii)))
        (vs! jg (+ im4 3) (+ (- is3 2) (vr j3 ii)))
        (flvs! ten im best))

      (let* ([ii (fx-- i0)]
             [i4 (* 4 i)]
             [best (vrzj ii)])
        (vs! jg i4 0)
        (vs! jg (+ i4 1) (+ (- is1 2) (vr j1 ii)))
        (vs! jg (+ i4 2) (+ (- is2 2) (vr j2 ii)))
        (vs! jg (+ i4 3) (+ (- is3 2) (vr j3 ii)))
        (flvs! ten i best))
      (values (fx-- i0) (fx-- i1)))])
  
    (for* ([i3 (in-range 1 n3)]
           [i2 (in-range 1 n2)]
           [i1 (in-range 1 n1)])
      (flvs! z (zidx i1 i2 i3) 0.0))

    (for ([i (in-range mm i0 -1)])
      (let ([i1 (fx-- i)])
        (flvs! z (zidx (vr j1 i1) (vr j2 i1) (vr j3 i1)) -1.0)))
    (for ([i (in-range mm i1 -1)])
      (let ([i1 (fx+ (fx-- i) mm)])
        (flvs! z (zidx (vr j1 i1) (vr j2 i1) (vr j3 i1)) 1.0))))

  (comm3 (CGSingle) z 0 n1 n2 n3))


(define (norm2u3 r n1 n2 n3 rnmu nx ny nz)
  (let-values ([(rnmu rnm2)
    (for*/fold ([rnmu 0.0]
                [rnm2 0.0]) 
          ([i3 (in-range 1 (fx-- n3))]
           [i2 (in-range 1 (fx-- n2))]
           [i1 (in-range 1 (fx-- n1))])
      (let ([rv (vr3 r i1 i2 i3 n1 n2)])
        (values (max rnmu (abs rv)) (+ rnm2 (* rv rv)))))])
    (sqrt (/ rnm2 (fx** nx ny nz)))))

(define (bubble ten j1 j2 j3 m ind)
  (define-syntax-rule (swapv v i0 i1 v0 v1)
    (begin
      (flvs! v i0 v1)
      (flvs! v i1 v0)))
  (define-syntax-rule (swap v i0 i1)
    (let ([t (vr v i0)])
      (vs! v i0 (vr v i1))
      (vs! v i1 t)))
  (define-syntax-rule (swapem cmp)
    (let ([l (fx-- m)])
      (let loop ([i 0])
        (when (i . < . l)
          (let* ([i0 (fx+ i (fx* m ind))]
                 [i1 (fx++ i0)]
                 [v0 (flvr ten i0)]
                 [v1 (flvr ten i1)])
            (when (v0 . cmp . v1)
              (swapv ten i0 i1 v0 v1)
              (swap j1 i0 i1)
              (swap j2 i0 i1)
              (swap j3 i0 i1)
              (loop (fx++ i))))))))

  (if (= ind 1)
    (swapem >)
    (swapem <)))

(define (stencil1 u ii1 i2 i3 n1 n3)
  (fl+ (fl+ (vr3 u ii1 (fx- i2 1) i3         n1 n3)
            (vr3 u ii1 (fx+ i2 1) i3         n1 n3))
       (fl+ (vr3 u ii1 i2         (fx- i3 1) n1 n3)
            (vr3 u ii1 i2         (fx+ i3 1) n1 n3))))

(define (stencil2 u ii1 i2 i3 n1 n3)
  (fl+ (fl+ (vr3 u ii1 (fx- i2 1) (fx- i3 1) n1 n3)
            (vr3 u ii1 (fx+ i2 1) (fx- i3 1) n1 n3))
       (fl+ (vr3 u ii1 (fx- i2 1) (fx+ i3 1) n1 n3)
            (vr3 u ii1 (fx+ i2 1) (fx+ i3 1) n1 n3))))

(define (resid cg a u v r off n1 n2 n3 nm u1 u2)
  (CGfor cg ([i3 (in-range 1 (fx-- n3))])
    (for ([i2 (in-range 1 (fx-- n2))])
      (for ([i1 (in-range n1)])
        (let ([ii1 (fx+ off i1)])
          (flvs! u1 i1 (stencil1 u ii1 i2 i3 n1 n3))
          (flvs! u2 i1 (stencil2 u ii1 i2 i3 n1 n3))))

    (for ([i1 (in-range 1 (fx-- n1))])
      (flvs! r (vidx off i1 i2 i3 n1 n3) 
             (fl- (vro3 v off i1 i2 i3 n1 n3)
                (fl+ (fl* (flvr a 0) (vro3 u off i1 i2 i3 n1 n3))
                (fl+ (fl* (flvr a 2) (fl+ (flvr u2 i1) (fl+ (flvr u1 (fx-- i1)) (flvr u1 (fx++ i1)))))
                     (fl* (flvr a 3) (fl+ (flvr u2 (fx-- i1)) (flvr u2 (fx++ i1)))))))))))
  (comm3 cg r off n1 n2 n3))

(define (psinv cg c r u off n1 n2 n3 nm r1 r2)
  (CGfor cg ([i3 (in-range 1 (fx-- n3))])
    (for ([i2 (in-range 1 (fx-- n2))])
      (for ([i1 (in-range n1)])
        (let ([ii1 (fx+ off i1)])
          (flvs! r1 i1 (stencil1 r ii1 i2 i3 n1 n2))
          (flvs! r2 i1 (stencil2 r ii1 i2 i3 n1 n2))))
      (for ([i1 (in-range 1 (fx-- n1))])
        (flvs! u (vidx off i1 i2 i3 n1 n2) 
          (fl+ 
             (fl+ (flvr u (vidx off i1 i2 i3 n1 n2))
                  (fl* (flvr c 0) (vro3 r off i1 i2 i3 n1 n2)))
             (fl+ (fl* (flvr c 1) (fl+ (vro3 r off (fx-- i1) i2 i3 n1 n2)
                                  (fl+ (vro3 r off (fx++ i1) i2 i3 n1 n2)
                                       (flvr r1 i1))))
                  (fl* (flvr c 2) (fl+ (flvr r2 i1) (fl+ (flvr r1 (fx-- i1)) (flvr r1 (fx++ i1)))))))))))

  (comm3 cg u off n1 n2 n3))

(define (mg3P cg c a u v r n1 n2 n3 ir m1 m2 m3 nm lt u1 u2 z1 z2 z3)
  (define lb 1)
  (for ([k (in-range (fx-- lt) (fx-- lb) -1)])
    (let ([j (fx- k 1 )])
      (rprj3 cg r (vr ir k) (vr m1 k) (vr m2 k) (vr m3 k) (vr ir j) (vr m1 j) (vr m2 j) (vr m3 j) nm u1 u2)))
  (let* ([k (fx- lb 1 )]
         [irk (vr ir k)]
         [m1k (vr m1 k)]
         [m2k (vr m2 k)]
         [m3k (vr m3 k)])
    (CG-n0-only cg
      (zero3 u irk m1k m2k m3k))
    (psinv cg c r u irk m1k m2k m3k nm u1 u2))
  (for ([k (in-range lb (fx-- lt))])
    (let* ([j (fx- k 1 )]
           [irk (vr ir k)]
           [m1k (vr m1 k)]
           [m2k (vr m2 k)]
           [m3k (vr m3 k)]
           [irj (vr ir j)]
           [m1j (vr m1 j)]
           [m2j (vr m2 j)]
           [m3j (vr m3 j)])
      (CG-n0-only cg
        (zero3 u irk m1k m2k m3k))
      (interp cg u irj m1j m2j m3j irk m1k m2k m3k z1 z2 z3)
      (CG-B cg)
      (resid  cg a u r r irk m1k m2k m3k nm u1 u2)
      (psinv  cg c r u irk m1k m2k m3k nm u1 u2)))
  (let ([j (fx- lt 2 )])
    (interp cg u (vr ir j) (vr m1 j) (vr m2 j) (vr m3 j) 0 n1 n2 n3 z1 z2 z3)
      (CG-B cg)
    (resid cg a u v r 0 n1 n2 n3 nm u1 u2)
    (psinv cg c r u 0 n1 n2 n3 nm u1 u2)))

(define (rprj3 cg r roff m1k m2k m3k soff m1j m2j m3j nm x1 y1)
  (let ([d1 (if (fx= m1k 3) 2 1)]
        [d2 (if (fx= m2k 3) 2 1)]
        [d3 (if (fx= m3k 3) 2 1)])
    (CGfor cg ([j3 (in-range 2 m3j)])
      (let ([i3 (fx- (fx- (fx* 2 j3) d3) 1)])
        (for ([j2 (in-range 2 m2j)])
          (let ([i2 (fx- (fx- (fx* 2 j2) d2) 1)])
            (for ([j1 (in-range 2 (fx++ m1j))])
              (let* ([i1 (fx- (fx- (fx* 2 j1) d1) 1)]
                     [ii1 (fx- (fx+ roff i1) 1)])
                (flvs! x1 (fx-- i1) (stencil1 r ii1 i2 i3 m1k m2k))
                (flvs! y1 (fx-- i1) (stencil2 r ii1 i2 i3 m1k m2k))))

            (for ([j1 (in-range 2 m1j)])
              (let* ([i1 (fx- (fx- (fx* 2 j1) d1) 1)]
                     [ii1 (fx+ roff i1)]
                     [x2 (stencil1 r ii1 i2 i3 m1k m2k)]
                     [y2 (stencil2 r ii1 i2 i3 m1k m2k)])

                (let ([a1 (fl* 0.5       (vro3 r roff i1        i2 i3 m1k m2k))]
                      [a2 (fl* 0.25   (fl+ (vro3 r roff (fx-- i1) i2 i3 m1k m2k)
                                       (fl+ (vro3 r roff (fx++ i1) i2 i3 m1k m2k)
                                       x2)))]
                      [a3 (fl* 0.125  (fl+ (flvr x1 (fx-- i1)) 
                                       (fl+ (flvr x1 (fx++ i1)) 
                                       y2)))]
                      [a4 (fl* 0.0625 (fl+ (flvr y1 (fx-- i1)) 
                                       (flvr y1 (fx++ i1))))])

                (flvs! r (vidx soff (fx-- j1) (fx-- j2) (fx-- j3) m1j m2j)
                         (fl+ (fl+ a1 a2) (fl+ a3 a4)))))))))))

  (comm3 cg r soff m1j m2j m3j))

(define (interp cg u zoff mm1 mm2 mm3 uoff n1 n2 n3 z1 z2 z3)
  (define-syntax-rule (vs!+ u idx v) (flvs! u idx (fl+ (flvr u idx) v)))
  (if (and (not (= n1 3))
           (not (= n2 3))
           (not (= n3 3)))
    (CGfor cg ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 (fx++ mm1))])
          (let* ([si1 (fx-- i1)]
                 [si2 (fx-- i2)]
                 [si3 (fx-- i3)]
                 [ii (fx+ zoff si1)]
                 [u_s (vr3 u ii  i2 si3 mm1 mm2)]
                 [uss (vr3 u ii si2 si3 mm1 mm2)]
                 [us_ (vr3 u ii si2  i3 mm1 mm2)]
                 [u__ (vr3 u ii  i2  i3 mm1 mm2)])
            (flvs! z1 si1 (fl+ u_s uss))
            (flvs! z2 si1 (fl+ us_ uss))
            (flvs! z3 si1 (fl+ u__ (fl+ us_ (flvr z1 (fx-- i1)))))))
        (for ([i1 (in-range 1 mm1)])
          (let* ([si1 (fx-- i1)]
                 [si2 (fx-- i2)]
                 [si3 (fx-- i3)]
                 [xi1 (fx- (fx* 2 i1) 2)]
                 [xi2 (fx- (fx* 2 i2) 2)]
                 [xi3 (fx- (fx* 2 i3) 2)]
                 [ii (fx+ uoff xi1)] 
                 [jj (fx+ zoff si1)] 
                 [uz1 (vr3 u jj        si2 si3 mm1 mm2)]
                 [uz2 (vr3 u (fx++ jj) si2 si3 mm1 mm2)])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) uz1)
            (vs!+ u (vidx3 (fx++ ii) xi2 xi3 n1 n2) (fl* 0.5 (+ uz2 uz1)))))
        (for ([i1 (in-range 1 mm1)])
          (let* ([xi1 (fx- (fx* 2 i1) 2)]
                 [xi2 (fx- (fx* 2 i2) 1)]
                 [xi3 (fx- (fx* 2 i3) 2)]
                 [ii (fx+ uoff xi1)] 
                 [zvs (flvr z1 (fx-- i1))])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) (fl* 0.5 zvs))
            (vs!+ u (vidx3 (fx++ ii) xi2 xi3 n1 n2) (fl* 0.25 (fl+ zvs (flvr z1 i1))))))
        (for ([i1 (in-range 1 mm1)])
          (let* ([xi1 (fx- (fx* 2 i1) 2)]
                 [xi2 (fx- (fx* 2 i2) 2)]
                 [xi3 (fx- (fx* 2 i3) 1)]
                 [ii (fx+ uoff xi1)] 
                 [zvs (flvr z2 (fx-- i1))])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) (fl* 0.5 zvs))
            (vs!+ u (vidx3 (fx++ ii) xi2 xi3 n1 n2) (fl* 0.25 (fl+ zvs (flvr z2 i1))))))
        (for ([i1 (in-range 1 mm1)])
          (let* ([xi1 (fx- (fx* 2 i1) 2)]
                 [xi2 (fx- (fx* 2 i2) 1)]
                 [xi3 (fx- (fx* 2 i3) 1)]
                 [ii (fx+ uoff xi1)] 
                 [zvs (flvr z3 (fx-- i1))])
            (vs!+ u (vidx3 ii        xi2 xi3 n1 n2) (fl* 0.25 zvs))
            (vs!+ u (vidx3 (fx++ ii) xi2 xi3 n1 n2) (fl* 0.125 (fl+ zvs (flvr z3 i1))))))))
    (let ([d1 (if (= n1 3) 2 1)]
          [t1 (if (= n1 3) 1 0)]
          [d2 (if (= n2 3) 2 1)]
          [t2 (if (= n2 3) 1 0)]
          [d3 (if (= n3 3) 2 1)]
          [t3 (if (= n3 3) 1 0)])
    (define-syntax-rule (u!idx uoff i1 i2 i3 d1 d2 d3 n1 n2)
      (vidx3 (fx+ uoff (fx- (fx* 2 i1) 1 d1)) (fx- (fx* 2 i2) 1 d2 ) (fx- (fx* 2 i3) 1 d3) n1 n2))
    (CGfor cg ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 d2 d3 n1 n2)
                  (vro3 u zoff (fx-- i1) (fx-- i2) (fx-- i3) mm1 mm2)))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 d2 d3 n1 n2)
                  (fl* 0.5  (fl+ (vro3 u zoff i1  si2 si3 mm1 mm2)
                             (vro3 u zoff si1 si2 si3 mm1 mm2)))))))
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 t2 d3 n1 n2)
                  (fl* 0.5  (fl+ (vro3 u zoff si1 i2  si3 mm1 mm2)
                             (vro3 u zoff si1 si2 si3 mm1 mm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 t2 d3 n1 n2)
                  (fl* 0.25 (fl+ (vro3 u zoff i1  i2  si3 mm1 mm2)
                             (vro3 u zoff i1  si2 si3 mm1 mm2)
                             (vro3 u zoff si1 i2  si3 mm1 mm2)
                             (vro3 u zoff si1 si2 si3 mm1 mm2))))))))

    (CGfor cg ([i3 (in-range 1 mm3)])
      (for ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 d2 t3 n1 n2)
                  (fl* 0.5   (fl+ (vro3 u zoff si1 si2 i3   mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 d2 t3 n1 n2)
                  (fl* 0.25  (fl+ (vro3 u zoff i1  si2 i3  mm1 mm2)
                              (vro3 u zoff si1 si2 i3  mm1 mm2)
                              (vro3 u zoff i1  si2 si3 mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2)))))))
      (CGfor cg ([i2 (in-range 1 mm2)])
        (for ([i1 (in-range 1 mm1)])
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 d1 t2 t3 n1 n2)
                  (fl* 0.25  (fl+ (vro3 u zoff si1 i2  i3  mm1 mm2)
                              (vro3 u zoff si1 si2 i3  mm1 mm2)
                              (vro3 u zoff si1 i2  si3 mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2))))))
        (for ([i1 (in-range 1 mm1)])
          (let ([si1 (fx-- i1)] [si2 (fx-- i2)] [si3 (fx-- i3)])
          (vs!+ u (u!idx uoff i1 i2 i3 t1 t2 t3 n1 n2)
                  (fl* 0.125 (fl+ (vro3 u zoff i1  i2  i3  mm1 mm2)
                              (vro3 u zoff i1  si2 i3  mm1 mm2)
                              (vro3 u zoff si1 i2  i3  mm1 mm2)
                              (vro3 u zoff si1 si2 i3  mm1 mm2)
                              (vro3 u zoff i1  i2  si3 mm1 mm2)
                              (vro3 u zoff i1  si2 si3 mm1 mm2)
                              (vro3 u zoff si1 i2  si3 mm1 mm2)
                              (vro3 u zoff si1 si2 si3 mm1 mm2))))))))))))

(define (comm3 cg u off n1 n2 n3)
  (CGfor cg ([i3 (in-range 1 (fx-- n3))])
    (for ([i2 (in-range 1 (fx-- n2))])
      (flvs! u (vidx off 0         i2 i3 n1 n2) (vro3 u off (fx- n1 2) i2 i3 n1 n2))
      (flvs! u (vidx off (fx-- n1) i2 i3 n1 n2) (vro3 u off 1          i2 i3 n1 n2))))

  (CGfor cg ([i3 (in-range 1 (fx-- n3))])
    (for ([i1 (in-range n1)])
      (flvs! u (vidx off i1 0         i3 n1 n2) (vro3 u off i1 (fx- n2 2) i3 n1 n2))
      (flvs! u (vidx off i1 (fx-- n2) i3 n1 n2) (vro3 u off i1 1          i3 n1 n2))))

  (CG-n0-only cg
    (for* ([i2 (in-range n2)]
           [i1 (in-range n1)])
      (flvs! u (vidx off i1 i2 0         n1 n2) (vro3 u off i1 i2 (fx- n3 2) n1 n2))
      (flvs! u (vidx off i1 i2 (fx-- n3) n1 n2) (vro3 u off i1 i2 1          n1 n2)))))
