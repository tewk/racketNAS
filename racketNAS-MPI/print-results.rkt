#lang racket/base

(provide print-results-fortran)

(define (print-results-fortran name class n1 n2 n3 niter nprocs-compiled nprocs-total
         t mops optype verified npbversion)

  (printf " ~a Benchmark Completed\n" name)
  (printf " Class           = ~a\n" class)
  (printf " Iterations      = ~a\n" niter)
  (printf " Time in seconds = ~a\n" t)
  (printf " Total Processes = ~a\n" nprocs-total)
  (printf " Compiled procs  = ~a\n" nprocs-compiled)
  (printf " Mop/s total     = ~a\n" mops)
  (printf " Mop/s/process   = ~a\n" (/ mops nprocs-total))
  (printf " Operation type  = ~a\n" optype)
  (printf " Verification    = ~a\n" (if verified "SUCCESSFFUL" "UNSUCCESSFUL"))
  (printf " Version         = ~a\n" npbversion)
  (printf " Please send feedbacks to\n")
  (printf " tewk@cs.utah.edu\n"))
