#lang racket/base 
(require racket/flonum)
(provide random-init 
         randlc 
         ;randlc/a 
         ;randlc/2
         vranlc 
         ipow46
         print-seed
         krandlc/a
         krandlc/2
         (rename-out [krandlc/a randlc/a]
                     [krandlc/2 randlc/2])

         new-rand
         power/r)

(define-struct rand (seed tran) #:mutable)
(define (new-rand) (make-rand 0.0 314159265.0))

(define seed 0.0) 

;Default seed
(define tran 314159265.0)  ;First 9 digits of PI  

;Constants 
(define d2m46 (expt 0.5 46)) 
(define i246m1 (- (floor (expt 2 46)) 1))

(define (random-init sd) (set! seed sd))
(define (print-seed) (printf "~a~n" seed))

(define r23 (expt 0.5 23))
(define t23 (expt 2.0 23))
(define r46 (expt r23 2))
(define t46 (expt t23 2))

;;randlc/a : double double -> double
(define (randlc/a x a) 
  (let* ( ;Break A into two parts such that A = 2^23 * A1 + A2 
         [a1 (floor (* r23 a))]
         [a2 (- a (* t23 a1))]
         ;Break X into two parts such that X = 2^23 * X1 + X2, compute 
         ;Z = A1 * X2 + A2 * X1 (mod 2^23), and then 
         ;X = 2^23 * Z + A2 * X2 (mod 2^46)
         [x1 (floor (* r23 x))]
         [x2 (- x (* t23 x1))]
         [t1 (+ (* a1 x2) (* a2 x1))]
         [t2 (floor (* r23 t1))]
         [z (- t1 (* t23 t2))]
         [t3 (+ (* t23 z) (* a2 x2))]
         [t4 (floor (* r46 t3))]
         [r (- t3 (* t46 t4))])
    r))

(define (krandlc/a x a)
  (exact->inexact (bitwise-and (* (inexact->exact (floor x)) (inexact->exact (floor a))) #x3FFFFFFFFFFF)))


(define (krandlc/2 x a)
  (define nx (krandlc/a x a))
  (values nx (* r46 nx)))

;(define (gen-randlc ))
(define (randlc/2 x a)
  (define nx (randlc/a x a))
  (values nx (* r46 nx)))

;;randlc : double -> double
(define (randlc a) 
  (let ([r (randlc/a tran a)])
    (set! tran r)
    (* r46 r)))

;;vranlc : double double double vector-of-double int -> double
(define (vranlc n x a y offset) 
  (let ([Lx (inexact->exact (floor x))] 
        [La (inexact->exact (floor a))]) 
    (for ([i (in-range 0 n)]) 
      (set! Lx (bitwise-and (inexact->exact (* Lx La)) i246m1)) 
      (flvector-set! y 
                   (+ offset i) 
                   (* d2m46 Lx))) 
    (exact->inexact Lx)))

;;ipow46 : double int -> double
(define (ipow46 a exponent) 
  ;Use: 
  ; a^n = a^(n / 2) * a^(n / 2) if n is even, else 
  ; a^n = a * a^(n - 1)         if n is odd
  (if (= 0 exponent) 
      seed
      (let loop ([n exponent] 
                 [q a] 
                 [r 1])
        (if (<= n 1)
            (begin 
              (set! seed (randlc/a r q))
              seed)
            (if (even? n)
                (loop (quotient n 2) (randlc/a q q) r)
                (loop (- n 1) q (randlc/a r q)))))))

(define (power/r rng a n)
  (let loop ([pow 1.0]
             [seed (rand-seed rng)]
             [nj n]
             [aj a])
    (if (not (zero? nj))
      (let* ([njmod2 (= (modulo nj 2) 1)]
             [seed (if njmod2 (randlc/a pow aj) seed)]
             [pow  (if njmod2 seed pow)])
          (let ([seed (randlc/a aj aj)])
            (loop pow seed (quotient nj 2) seed)))
      (begin
        (set-rand-seed! rng seed)
        pow))))

(module+ test
  (time
    (for ([i 1000000])
      (krandlc/a 314159265.0 1220703125.0)))
  (time
    (for ([i 1000000])
      (randlc/a 314159265.0 1220703125.0)))
  (time
    (for ([i 1000000])
      (krandlc/2 314159265.0 1220703125.0)))
  (time
    (for ([i 1000000])
      (randlc/2 314159265.0 1220703125.0))))
