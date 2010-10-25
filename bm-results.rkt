#lang racket
(provide (struct-out BMResults) 
         new-BMResults 
         print-verification-status 
         print-results) 

(define-struct BMResults 
  (name 
   machine-name 
   pr-lang 
   clss 
   n1 
   n2 
   n3 
   niter 
   time 
   acctime 
   wctime 
   mops 
   tm-sent 
   tm-received 
   rec-arr-size 
   optype 
   numthreads 
   serial 
   pid 
   verified) 
  #:mutable)

;Convenience method that matches the 
;BMResults ctor definition in the Java version.
(define (new-BMResults bmname 
                       clss 
                       bn1 
                       bn2 
                       bn3 
                       bniter 
                       btime 
                       bmops 
                       boptype 
                       passed-verification 
                       bserial 
                       num-threads 
                       bid) 
  (make-BMResults bmname 
                  "Machine Name?" 
                  "PLT Racket" 
                  clss 
                  bn1 
                  bn2 
                  bn3 
                  bniter 
                  btime 
                  0 
                  0 
                  bmops 
                  0 
                  0 
                  0 
                  boptype 
                  num-threads 
                  bserial 
                  bid 
                  passed-verification))
                       


;;print-verification-status : char number string -> void
(define (print-verification-status cls verified bmname) 
  (cond 
    [(or (char=? cls #\U) (= verified -1)) 
     (begin 
       (printf "Problem size unknown~n") 
       (printf "~a.~a: Verification Not Performed~n" bmname cls))] 
    [(= verified 1) 
     (printf "~a.~a: Verification Successful~n" bmname cls)] 
    [else 
     (printf "~a.~a: Verification Failed~n" bmname cls)]))

(define (widtho . r)
  (let* ([s (apply format r)]
         [sl (string-length s)])
    (printf "~a\n" (string-append s (make-string (- 63 sl) #\ ) "* "))))

;;print-results : BMResults -> void
(define (print-results results) 
  (save-results results)
  (printf "***** NAS Parallel Benchmarks PLT Racket version 0.1 ~a*********~n" 
          (BMResults-name results)) 
  (widtho "* Class             = ~a" (BMResults-clss results)) 
  (if (and (= (BMResults-n2 results) 0) (= (BMResults-n3 results) 0)) 
      (widtho "* Size              = ~a" (BMResults-n1 results)) 
      (widtho "* Size              = ~a X ~a X ~a" 
              (BMResults-n1 results) 
              (BMResults-n2 results) 
              (BMResults-n3 results))) 
  (widtho "* Iterations        = ~a" (BMResults-niter results)) 
  (widtho "* Time in seconds   = ~a" (BMResults-time results)) 
  (widtho "* ACCTime           = ~a" (BMResults-acctime results)) 
  (widtho "* Mops total        = ~a" (BMResults-mops results)) 
  (widtho "* Operation type    = ~a" (BMResults-optype results)) 
  (case (BMResults-verified results)
    [(#t) (widtho "* Verification      = Successful")] 
    [(#f) (widtho "* Verification      = Failed")] 
    [(1)  (widtho "* Verification      = Successful")] 
    [(0)  (widtho "* Verification      = Failed")] 
    [else (widtho "* Verification      = Not Performed")])
  (unless (BMResults-serial results) 
    (widtho "* Threads requested = ~a" (BMResults-numthreads results)))
  (widtho "*")
  (widtho "* Please send all errors/feedbacks to:")
  (widtho "* NPB Racket Working Team")
  (widtho "* tewk@cs.utah.edu")
  (printf "****************************************************************\n"))

(define (save-results results) 
  (call-with-output-file "racketNAS_Results.rktd" #:exists 'append (lambda (o)
    (write (list
      (BMResults-name results)
      (BMResults-clss results)
      (if (and (= (BMResults-n2 results) 0) (= (BMResults-n3 results) 0)) 
        (list (BMResults-n1 results))
        (list (BMResults-n1 results) 
              (BMResults-n2 results) 
              (BMResults-n3 results))) 
      (BMResults-niter results)
      (BMResults-time results)
      (BMResults-serial results)
      (BMResults-numthreads results)
      (BMResults-verified results)
      (BMResults-mops results)) o)
      (newline o))))
