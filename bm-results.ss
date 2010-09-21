(module bm-results scheme 
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
  
  ;;print-results : BMResults -> void
  (define (print-results results) 
    (printf "***** NAS Parallel Benchmarks PLT Racket version 0.1 ~a*****~n" 
            (BMResults-name results)) 
    (printf "* Class = ~a~n" (BMResults-clss results)) 
    (if (and (= (BMResults-n2 results) 0) (= (BMResults-n3 results) 0)) 
        (printf "* Size = ~a~n" (BMResults-n1 results)) 
        (printf "* Size = ~a X ~a X ~a~n" 
                (BMResults-n1 results) 
                (BMResults-n2 results) 
                (BMResults-n3 results))) 
    (printf "* Iterations = ~a~n" (BMResults-niter results)) 
    (printf "* Time in seconds = ~a~n" (BMResults-time results)) 
    (printf "* ACCTime = ~a~n" (BMResults-acctime results)) 
    (printf "* Mops total = ~a~n" (BMResults-mops results)) 
    (printf "* Operation type = ~a~n" (BMResults-optype results)) 
    (case (BMResults-verified results)
      [(#t) (printf "* Verification = Successful~n")] 
      [(#f) (printf "* Verification = Failed~n")] 
      [(1) (printf "* Verification = Successful~n")] 
      [(0) (printf "* Verification = Failed~n")] 
      [else (printf "* Verification = Not Performed~n")])
    (if (not (BMResults-serial results)) 
        (printf "* Threads requested = ~a~n" (BMResults-numthreads results)) 
        void))
  
  
  
  
  
  )
