(module rand-generator scheme
  (provide random-init 
           randlc 
           randlc/a 
           vranlc 
           ipow46
           print-seed)
  
  (define seed 0.0) 
  
  ;Default seed
  (define tran 314159265.0)  ;First 9 digits of PI  
  
  ;Constants 
  (define d2m46 (expt 0.5 46)) 
  (define i246m1 (- (floor (expt 2 46)) 1))
  
  (define (random-init sd) (set! seed sd))
  (define (print-seed) (printf "~a~n" seed))
  
  ;;randlc/a : double double -> double
  (define (randlc/a x a) 
    (let* ([r23 (expt 0.5 23)]
           [r46 (expt r23 2)]
           [t23 (expt 2.0 23)]
           [t46 (expt t23 2)]
           ;Break A into two parts such that A = 2^23 * A1 + A2 
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
  
  ;;randlc : double -> double
  (define (randlc a) 
    (let ([r23 0.0] 
          [r46 0.0] 
          [t23 0.0] 
          [t46 0.0] 
          [t1 0.0]
          [t2 0.0] 
          [t3 0.0] 
          [t4 0.0] 
          [a1 0.0] 
          [a2 0.0] 
          [x1 0.0] 
          [x2 0.0] 
          [z 0.0]) 
      (set! r23 (expt 0.5 23)) 
      (set! r46 (expt r23 2)) 
      (set! t23 (expt 2.0 23))
      (set! t46 (expt t23 2)) 
      ;Break A into two parts such that A = 2^23 * A1 + A2 
      (set! t1 (* r23 a)) 
      (set! a1 (floor t1))
      (set! a2 (- a (* t23 a1))) 
      ;Break X into two parts such that X = 2^23 * X1 + X2, compute 
      ;Z = A1 * X2 + A2 * X1 (mod 2^23), and then 
      ;X = 2^23 * Z + A2 * X2 (mod 2^46) 
      (set! t1 (* r23 tran)) 
      (set! x1 (floor t1))
      (set! x2 (- tran (* t23 x1))) 
      (set! t1 (+ (* a1 x2) (* a2 x1))) 
      (set! t2 (floor (* r23 t1))) 
      (set! z (- t1 (* t23 t2))) 
      (set! t3 (+ (* t23 z) (* a2 x2))) 
      (set! t4 (floor (* r46 t3))) 
      (set! tran (- t3 (* t46 t4)))
      (* r46 tran)))
  
  ;;vranlc : double double double vector-of-double int -> double
  (define (vranlc n x a y offset) 
    (let ([Lx (inexact->exact (floor x))] 
          [La (inexact->exact (floor a))]) 
      (for ([i (in-range 0 n)]) 
        (set! Lx (bitwise-and (inexact->exact (* Lx La)) i246m1)) 
        (vector-set! y 
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
  
)
