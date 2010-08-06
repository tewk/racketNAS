(module bm-args scheme  
  (provide (struct-out BMArgs) 
           parse-cmd-line-args 
           print-banner
           string-contains)
  
  (define-struct BMArgs 
    (class num-threads serial) 
    #:mutable)
  
  ;Helpers 
  (define (string-contains str terms) 
    (cond 
      [(empty? terms) #f] 
      [else 
       (if (regexp-match (first terms) str) 
           #t 
           (string-contains str (rest terms)))]))
  
  (define (print-command-line-error bmname) 
    (printf "Summarize command options here~n"))
  
  (define (print-args bmargs) 
    (printf "Class: ~a~n" (BMArgs-class bmargs)) 
    (printf "Num. Threads: ~a~n" (BMArgs-num-threads bmargs)) 
    (printf "Serial: ~a~n" (BMArgs-serial bmargs)))
  
  ;;print-banner : string BMArgs -> void
  (define (print-banner bmname args) 
    (printf "NAS Parallel Benchmarks PLT Scheme Version 0.1~n") 
    (if (BMArgs-serial args) 
        (printf "Serial Version ~a.~a~n" bmname (BMArgs-class args)) 
        (printf "Multithreaded Version ~a.~a np=~a~n" 
                bmname 
                (BMArgs-class args) 
                (BMArgs-num-threads args))))


  ;;parse-cmd-line-args : list-of-string string -> BMArgs
  (define (parse-cmd-line-args argv bmname) 
    (let ([ret (make-BMArgs #\U 1 #t)])
      (for ([arg argv]) 
        (cond 
          [(string-contains arg '("SERIAL" "serial" "-serial" "-SERIAL")) (set-BMArgs-serial! ret #t)]
          [(string-contains arg '("class=" "CLASS=" "-class" "-CLASS"))
           (begin
             (if (> (string-length arg) 6) 
                 (set-BMArgs-class! ret (char-upcase (string-ref arg 6))) 
                 void)
             (case (BMArgs-class ret) 
               [(#\A #\B #\C #\S #\W) void]
               [else 
                (begin 
                  (printf "Classes allowed are A, B, C, W, and S.~n") 
                  (print-command-line-error bmname))]))]
          [(string-contains arg '("np=" "NP=" "-NP" "-np")) 
           (if (> (string-length arg) 3) 
               (begin 
                 (set-BMArgs-num-threads! ret (string->number (substring arg 3))) 
                 (set-BMArgs-serial! ret #f)) 
               void)]
          [else 
           (printf "Invalid command line argument: ~a~n" arg)]))
      ret))
           

  )