;; EarSpring.lisp -- Solution of the EarSpring and Conductor equations
;; for restoration of sensioneural hearing impairment
;;
;; DM/RAL 01/07

;; -----------------------------------------------------------------------------
(defparameter *top-function* nil)

(defun chktype (val type-symbol)
  (unless (eq type-symbol (type-of val))
    (error "~A expected in ~A (got ~A)" type-symbol *top-function* val))
  val)

(defmacro chkarg ((arg type) &body body)
  `(let ((,arg (chktype ,arg ',type)))
     (declare (type ,type ,arg))
     ,@body))

(defun dfloat (x)
  (if (realp x)
      (the double-float (coerce x 'double-float))
    (error "Real number expected in ~A (got ~A)" *top-function* x)))

;; ----------------------------------------------------------
;; Bracket a root given two starting values
;;
(defun zbrac (func x1 x2 &key (factor 1.6d0) (ntry 50))
  (let ((*top-function* 'zbrac))
    (declare (special *top-function*))
    
    (chkarg
        (ntry fixnum)
      
      (multiple-value-bind (x1 x2)
          (values (dfloat (min x1 x2))
                  (dfloat (max x1 x2)))
        (declare (type double-float x1 x2))
        
        (if (>= x1 x2)
            (error " initial range in zbrac must be non-empty. ")) 
        
        (let ((factor (dfloat factor))
              (succes nil)
              (f1 (dfloat (funcall func x1)))
              (f2 (dfloat (funcall func x2))))
          (declare (type symbol succes))
          (declare (type double-float factor f1 f2))
          
          (do ((j ntry (1- j)))
              ((or (setf succes (minusp (* f1 f2)))
                   (not (plusp j))))
            (declare (type fixnum j))
            
            (cond
             ((< (abs f1) (abs f2))
              (setf f2 f1
                    x1 (+ x1 (* factor (- x1 (shiftf x2 x1))))
                    f1 (dfloat (funcall func x1))) )
             
             (t
              (setf f1 f2
                    x2 (+ x2 (* factor (- x2 (shiftf x1 x2))))
                    f2 (dfloat (funcall func x2))) )
             
             ))
          
          (values x1 x2 succes)
          )))))

;; -------------------------------------------------------------------------
;; Root by Bisection -- root needs to be bracketed
;;
(defun rtbis (func x1 x2 &key (xacc 0) (jmax 60))
  ;; double precision floats have 56 bit mantissas
  (let ((*top-function* 'rtbis))
    (declare (special *top-function*))
    
    (chkarg
        (jmax fixnum)
      
      (multiple-value-bind (x1 x2 xacc)
          (values (dfloat (min x1 x2))
                  (dfloat (max x1 x2))
                  (dfloat xacc))
        (declare (type double-float x1 x2 xacc))
        
        (let ((fmid (dfloat (funcall func x2))))
          (declare (type double-float fmid))
          
          (if (zerop fmid)
              (the double-float x2)
            
            (let ((f  (dfloat (funcall func x1))))
              (declare (type double-float f))
              
              (if (zerop f)
                  (the double-float x1)

                (progn
                  (unless (minusp (* f fmid))
                    (error " root must be bracketed for bisection in rtbis "))
                  
                  (multiple-value-bind (rtbis dx)
                      (if (minusp f)
                          (values x1 (- x2 x1))
                        (values x2 (- x1 x2)))
                    (declare (type double-float rtbis dx))
                    
                    (do ((j jmax (1- j)))
                        ((or (= rtbis (+ rtbis dx))
                             (< (abs dx) xacc)
                             (zerop fmid)
                             (and (not (plusp j))
                                  (error "Too many iterations in rtbis"))))
                      (declare (type fixnum j))
                      
                      (setf dx (* dx 0.5d0))
                      (let ((xmid (+ rtbis dx)))
                        (declare (type double-float xmid))
                        
                        (setf fmid (dfloat (funcall func xmid)))
                        
                        (unless (plusp fmid)
                          (setf rtbis xmid))
                        ))
                    
                    (the double-float rtbis)
                    
                    )))
              ))
          ))
      )))
  
;-------------------------------------------------------------------------
;; Root by Secant Method -- root does not need to be bracketed, but it helps...
;;
(defun rtsec (func x1 x2 &key (xacc 0) (maxit 30))
  (let* ((*top-function* 'rtsec))
    (declare (special *top-function*))

    (chkarg
        (maxit fixnum)
      
      (multiple-value-bind (x1 x2 xacc)
          (values (dfloat x1)
                  (dfloat x2)
                  (dfloat xacc))
        (declare (type double-float x1 x2 xacc))

        (let* ((fl    (dfloat (funcall func x1)))
               (f     (dfloat (funcall func x2)))
               (rtsec 0d0)
               (dx    (- x1 x2))
               (xl    0d0))
          (declare (type double-float fl f rtsec dx xl))
          
          (if (< (abs fl) (abs f))
              (let ((swap fl))
                (declare (type double-float swap))
                (setf rtsec x1
                      xl    x2
                      fl    f
                      f     swap))
            (setf xl x1
                  rtsec x2))
          
          (do ((j 1 (1+ j)))
              ((or (and (> j maxit)
                        (error "rtsec exceeded maximum iterations"))
                   (< (abs dx) xacc)
                   (zerop f)
                   (= rtsec (+ rtsec dx))))
            (declare (type fixnum j))

            (setf dx  (/ (* (- xl rtsec) f) (- f fl))
                  xl  rtsec
                  fl  f
                  rtsec (+ rtsec dx)
                  f  (dfloat (funcall func rtsec))
                  ))
          
          (the double-float rtsec)

          )))))

;-------------------------------------------------------------------------
;; Root by Brent's Method -- root must be bracketed
;;
(defun zbrent (func x1 x2 &key (tol 0) (itmax 100) (eps 1.0d-16))
  (let* ((*top-function* 'zbrent))
    (declare (special *top-function*))

    (chkarg
        (itmax fixnum)

      (multiple-value-bind (a b tol eps)
          (values (dfloat x1)
                  (dfloat x2)
                  (dfloat tol)
                  (dfloat eps))
        (declare (type double-float a b tol eps))

        (let* ((fa (dfloat (funcall func a)))
               (fb (dfloat (funcall func b))))
          (declare (type double-float fa fb))

          (when (plusp (* fa fb))
            (error "Root must be bracketed for zbrent"))

          (let* ((c    b)
                 (fc   fb)
                 (d    0d0)
                 (e    0d0)
                 (s    0d0)
                 (p    0d0)
                 (q    0d0)
                 (r    0d0)
                 (xm   (- a b))
                 (tol1 tol))
            (declare (type double-float c fc d e s xm tol1 p q r))
            
            (do ((iter 1 (1+ iter)))
                ((or (and (> iter itmax)
                          (error "Zbrent exceeded maximum iterations"))
                     (zerop fb)
                     (<= (abs xm) tol1)))
              (declare (type fixnum iter))
              
              (when (plusp (* fb fc))
                (setf c a
                      fc fa
                      d (- b a)
                      e d))
              
              (when (< (abs fc) (abs fb))
                (setf a b
                      b c
                      c a
                      fa fb
                      fb fc
                      fc fa))
              
              (setf tol1 (+ (* 2d0 eps (abs b)) (* 0.5d0 tol))
                    xm   (* 0.5d0 (- c b)))

              (cond ((and (>= (abs e) tol1)
                          (zerop fb))
                     (setf s (/ fb fa))
                     
                     (if (= a c)
                         (setf p (* 2d0 xm s)
                               q (- 1d0 s))
                       (setf q (/ fa fc)
                             r (/ fb fc)
                             p (* s (+ (* 2d0 xm q (- q r))
                                       (* (- a b) (- r 1d0))))
                             q  (* (- q 1d0) (- r 1d0) (- s 1d0))
                             ))
                     
                     (if (plusp p)
                         (setf q (- q)
                               p (abs p)))
                     
                     (if (< (* 2d0 p)
                            (min (+ (* 3d0 xm q)
                                    (- (abs (* tol1 q)) (abs (* e q)))
                                    )))
                         (setf e d
                               d (/ p q))
                       (setf d xm
                             e d)))

                    
                    (t (setf d xm
                             e d))
                    )
              
               (setf a b
                     fa fb)
               
               (if (> (abs d) tol1)
                   (setf b (+ b d))
                 (setf b (if (minusp xm)
                             (- b tol1)
                           (+ b tol1))
                       ))

               (setf fb (dfloat (funcall func b)))
               )
            (the double-float b)
            )))
      )))

;-------------------------------------------------------------------------
#|
;; not cleaned up yet...
(defun rtsafe (func funcd x1 x2 xacc &key (maxit 100))
  (declare (type fixnum maxit))
  
  (prog ((xd1 (dfloat x1))
         (xd2 (dfloat x2))
         (xdacc (dfloat xacc))
         (fl 0d0) (fh 0d0) (rtsafe 0d0) (dx 0d0) (dxold 0d0) 
         (xh 0d0) (xl 0d0) (f 0d0) (df 0d0) (temp 0d0))
    (declare (type double-float xd1 xd2 xdacc fl fh rtsafe dx dxold xh f xl df temp))
    
    (setq fl (dfloat (funcall func xd1)))
    (if (zerop fl)
        (return (the double-float xd1)))
    
    (setq fh (dfloat (funcall func xd2)))
    (if (zerop fh)
        (return (the double-float xd2)))
    
    (if (> (* fl fh) 0d0) (error " root must be bracketed for rtsafe ")) 

    (cond 
     ((< fl 0d0)
      (setf xl xd1) 
      (setf xh xd2))
     (t 
      (setf xh xd1) 
      (setf xl xd2))) 
    
    (setf rtsafe (* 0.5d0 (+ xd1 xd2))) 
    (setf dxold (abs (- xd2 xd1))) 
    (setf dx dxold) 
    (setq f  (dfloat (funcall func  rtsafe))
          df (dfloat (funcall funcd rtsafe)))
    (do ((j 1 (+ j 1)))
        ((> j maxit) t)
      (cond 
       ((or (>= (* (+ (* (+ rtsafe (- xh)) df) (- f))
                   (+ (* (+ rtsafe (- xl)) df) (- f)))
                0d0)
            (> (abs (* 2d0 f)) (abs (* dxold df))))
        (setf dxold dx)
        (setf dx (* 0.5d0 (+ xh (- xl)))) 
        (setf rtsafe (+ xl dx)) 
        (if (= xl rtsafe) (go end)))
       (t 
        (setf dxold dx)
        (setf dx (/ f df)) 
        (setf temp rtsafe) 
        (setf rtsafe (+ rtsafe (- dx)))
        (if (= temp rtsafe) (go end))))
      
      (if (< (abs dx) xdacc) (go end))
      (setq f  (dfloat (funcall func rtsafe))
            df (dfloat (funcall funcd rtsafe)))
      (if 
          (< f 0d0) 
          (setf xl rtsafe)
        (setf xh rtsafe))) 
    
    (error " rtsafe exceeding maximum iterations ")
end
    (return (the double-float rtsafe))))
|#

;-------------------------------------------------------------------------
(defun db20 (x)
  (if (plusp x)
      (* 20.0d0 (log (abs x) 10.0d0))
    -140.0))

(defun db10 (x)
  (if (plusp x)
      (* 10.0d0 (log (abs x) 10.0d0))
    -140.0))

(defun ampl20 (x)
  (expt 10.0d0 (/ x 20.0d0)))

(defun ampl10 (x)
  (expt 10.0d0 (/ x 10.0d0)))

(defun sqr (x)
  (* x x))

(defun abserror (a b)
  (abs (- a b)))

(defun relerror (a b)
  (if (zerop b)
      (abs a)
    (abserror (/ a b) 1.0d0)))

(defun iterate (p f x)
  (um:with-tail-pure-code
    (if (funcall p x)
        x
      (iterate p f (funcall f x))
      )))

(defun deriv (fn dx)
  (let ((half-dx (* 0.5d0 dx)))
    (lambda (x)
      (if (zerop x)
          (/ (- (funcall fn half-dx)
                (funcall fn (- half-dx)))
             dx)
        (/ (- (funcall fn (* x (+ 1.0d0 half-dx)))
              (funcall fn (* x (- 1.0d0 half-dx))))
           (* x dx))
        ))))

(defun newton (fn start dx eps)
  (labels ((is-ok (x)
             (< (abs (funcall fn x)) eps))
           (do-better (x)
             (let* ((fder (deriv fn dx)))
               (- x (/ (funcall fn x) (funcall fder x)))
               )))
    (iterate #'is-ok #'do-better start)))

(defun find-root (fn start &optional start2)
  #| |#
  (multiple-value-bind (start1 start2)
      (if start2
          (values start start2)
        (zbrac fn start (* start 1.1)))
    (zbrent fn start1 start2))
  #| |#
  ;; (newton fn start 0.01 1.0d-8)
  )

;;
;; Iterative solution for damping beta and BigGamma = 1/2*gamma*(abs a40)^2
;; Based on 75 cent detuning of tone pitch from 40 dBSPL to 90 dBSPL
;; and on the threshold Sones level of (1/22)^2
;;
(destructuring-bind (beta-sol bigGamma-sol s90-sol &rest rest)
    (let* ((pitchdev 75.0d0) ;; cents
           (f90sq    (sqr (expt 2.0d0 (/ pitchdev 1200.0d0))))
           (s0       (/ 22.0d0)))
      (labels ((is-ok (args)
                 (destructuring-bind (beta bigGamma s90
                                           prev-beta prev-bigGamma prev-s90)
                     args
                   (and (< (relerror beta prev-beta) 1.0d-8)
                        (< (relerror bigGamma prev-bigGamma) 1.0d-8)
                        (< (relerror s90 prev-s90) 1.0d-8))))
               
               (do-better (args)
                 (destructuring-bind (beta bigGamma s90 &rest others) args
                   (declare (ignore others))
                   (let* ((a  (/ (- f90sq 1.0d0) (- s90 f90sq))))
                     (labels ((b (beta)
                                (* a (- 1.0d0 (sqr beta))))
                              (betafn (beta)
                                (let ((b-of-beta (b beta)))
                                  (- (* 1.0d-4
                                        (/ (+ (* 4.0d0 (sqr beta))
                                              (sqr b-of-beta))
                                           (+ (* 4.0d0 (sqr beta))
                                              (sqr (* b-of-beta (sqr s0))))
                                           ))
                                     (sqr s0)))))
                       #|
                       (let ((win (plt:wset 'betafn :clear t)))
                         (plt:fplot win '(-0.0001 0.001) #'betafn :yrange '(-0.01 0.01))
                         (plt:plot  win (list beta beta) '(-100 100) :color :orange))
                       ;;(break)
                       |#
                       (let* ((beta-prime (find-root #'betafn beta)))
                         (labels ((s90fn (s90)
                                    (- s90
                                       (/ (* (ampl10 (- 90 40))
                                             (+ (* 4.0d0 (sqr beta-prime))
                                                (sqr (b beta-prime))))
                                          (+ (* 4.0d0 (sqr beta-prime))
                                             (sqr (* (b beta-prime) s90)))
                                          ))))
                           #|
                           (let ((win (plt:wset 's90fn :clear t)))
                             (plt:fplot win '(20 70) #'s90fn)
                             (plt:plot win (list s90 s90) '(-1000 1000) :color :orange))
                           ;; (break)
                           |#
                           (let* ((s90-prime (find-root #'s90fn s90)))
                             (list beta-prime (b beta-prime) s90-prime
                                   beta bigGamma s90)
                             )))
                       )))))
        (iterate #'is-ok #'do-better
                 (list 0.00022d0 0.00196d0 47.0d0 100.0d0 100.0d0 0.0d0))
        ))
  (declare (ignore rest))
  (defparameter beta     beta-sol)
  (defparameter bigGamma bigGamma-sol)
  (defparameter s90      s90-sol))

(defun cube (x)
  (* x x x))

(defun cubert (x)
  (expt x (/ 3.0d0)))

(let* ((a     (* 4.0d0 (sqr beta)))
       (b     (sqr bigGamma))
       (b27   (* 27.0d0 b b (+ a b)))
       (a108  (* 4.0d0 (cube (* 3.0d0 a b))))
       (den   (* 3.0d0 b (cubert 2.0d0)))
       (num   (* a (cubert 2.0d0))))
  
  (defun sones (dbphons)
    (let* ((p     (ampl10 (- dbphons 40.0d0)))
           (b27p  (* b27 p))
           (rad   (sqrt (+ a108 (sqr b27p))))
           (v     (cubert (+ b27p rad))))
      (- (/ v den) (/ num v)))))

  
(defun sdb (dbphons)
  (db10 (sones dbphons)))

#|
(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(0 115) #'sdb
             :title "EarSpring dBSones vs Phon"
             :xtitle "Input [Phon]"
             :ytitle "Loudness [dbSone]")
  (plt:fplot win '(0 115) (lambda (dbphons)
                            (db20 (max 0.01d0
                                       (- (ampl20 (sdb dbphons)) 1.5d0))))
             :color :red
             :watermarkfn nil))

(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(0 115) #'sdb
             :title "EarSpring dBSones vs Phon"
             :xtitle "Input [Phon]"
             :ytitle "Loudness [dbSone]"
             :watermarkfn nil
             :fullgrid nil)
  (plt:draw-text win "Linear Region"
                 '(:data 10 -20))
  (plt:draw-text win "Cube-Root Compression"
                 '(:data 60 5))
  )
|#

;; approximate average exponent of nonlinearity over 41 to 100 dBSPL
;; (based on LMS regression)
(defparameter alpha
  (let* ((dbphonss    (coerce
                   (loop for ix from 41 to 100 collect (float ix 1.0d0))
                   'vector))
         (dbsones (map 'vector #'sdb dbphonss))
         (num     (loop for y across dbphonss
                        and z across dbsones sum
                        (* (- y 40.0d0) z)))
         (den     (loop for y across dbphonss sum
                        (sqr (- y 40.0d0)))))
    (/ (* 2.0d0 num) den)))

(defun asones (dbphons &key (alpha alpha))
  ;; sones derived from assumed constant exponent of nonlinearity
  (expt (ampl20 (- dbphons 40.0d0)) alpha))

(defun asdb (dbphons &key (alpha alpha))
  (db10 (asones dbphons :alpha alpha)))

(defun frac-live-hc (dbelev)
  (/ (+ 1.0d0 (/ (sones dbelev) (sones 120.0d0)))))

(defun dbcorr (dbelev dbphons)
  ;; correction gain needed to hear dbphons as dbphons
  ;; with threshold elevation dbelev
  (if (zerop dbelev)
      0.0
    (let* ((frac     (frac-live-hc dbelev))
           (sonoff   (- (* frac (sones dbelev)) (sones 0.0d0)))
           (sdbphons (sones dbphons)))
      (labels ((sonfn (x)
                 (- (* frac (sones x)) sonoff sdbphons)))
        (- (find-root #'sonfn (max dbphons dbelev)) dbphons)
        ))))

(defun dbuncorr (dbelev dbphons)
  ;; apparent dbphons for actual input dbphons
  ;; given threshold elevation dbelev
  (cond ((= dbphons dbelev) 0.0d0)
        ((< dbphons dbelev) -100.0d0)
        (t  ;; dbphons > dbelev
            (let* ((frac    (frac-live-hc dbelev))
                   (sonoff  (- (* frac (sones dbelev)) (sones 0.0d0)))
                   (sonapp  (- (* frac (sones dbphons)) sonoff)))
              (labels ((sonfn (x)
                         (- sonapp (sones x))))
                (find-root #'sonfn dbphons)
                )))
        ))

(defun adbuncorr (dbelev dbphons &key (alpha alpha))
  ;; apparent dbphons for actual input dbphons
  ;; given threshold elevation dbelev
  (cond ((= dbphons dbelev) 0.0d0)
        ((< dbphons dbelev) -100.0d0)
        (t  ;; dbphons > dbelev
            (let* ((frac    (frac-live-hc dbelev))
                   (sonoff  (- (* frac (asones dbelev :alpha alpha))
                               (asones 0.0d0 :alpha alpha)))
                   (sonapp  (- (* frac (asones dbphons :alpha alpha)) sonoff)))
              (labels ((sonfn (x)
                         (- sonapp (asones x :alpha alpha))))
                (find-root #'sonfn dbphons)
                )))
        ))

(defun a-gain (pdb dbelev &key (alpha alpha))
  (* (/ 20 alpha)
     (log (+ 1 (expt 10 (* (/ alpha 20) (- dbelev pdb))))
          10)))

#|
(defun solve-incorrect (pdb dbelev)
  (let ((gn (a-gain pdb dbelev :alpha alpha)))
    (labels ((sonfn (x)
               (- (a-gain pdb x :alpha 0.4) gn)))
      (find-root #'sonfn dbelev))))


;; show excess gain from incorrect psychoacoustic model
;; with alpha = 0.2 at any given presentation level
;; as a function of true threshold elevation.
(let* ((dbelev  60)
       (pdbs    '(20 30 40 50 60 70 80 90))
       (dels    (loop for pdb in pdbs collect
                      (- (a-gain pdb dbelev :alpha 0.4)
                         (a-gain pdb dbelev)))))
  (plt:plot 'xxx pdbs dels :clear t))

;; show effective threshold elevation needed for
;; incorrect psychoacoustic model with alpha = 0.2
;; to be correct in gain at presentation level 60 dB
;; as a function of true threshold elevation
;;
;; looks like we need about:
;;
;; Eff-dBelev = 1.5 * (True-dBelev - 30)
;;
(let* ((dbpres  60)
       (pdbels  '(20 30 40 50 60 70 80 90))
       (dels    (loop for pdbelev in pdbels collect
                      (solve-incorrect dbpres pdbelev))))
  (plt:plot 'xxx pdbels dels :clear t))

;; adjusting the threshold elevation of a hearing aid
;; with too-low an exponent (0.2) to be correct at some
;; input presentation level.
(let ((dbelev 50)
      (dbpres 50))
  (plt:fplot 'plt '(20 100) #'identity
             :clear t
             :color :darkgreen
             :thick 2
             :title "Effect of Incorrect Exponent"
             :xtitle "Presentation Level [dBHL]"
             :ytitle "Perceived Level [dBHL]"
             :legend "Normal Hearing"
             )
  (plt:paramplot 'plt '(20 100)
                 (lambda (pdb)
                   (+ pdb (a-gain pdb dbelev)))
                 #'identity
                 :color :blue
                 :thick 2
                 :legend (format nil "~A dB Elev, correct alpha" dbelev))
  (plt:paramplot 'plt '(20 100)
                 (let ((dbelev (solve-incorrect dbpres dbelev)))
                   (lambda (pdb)
                     (+ pdb (a-gain pdb dbelev :alpha 0.4))))
                 #'identity
                 :color :red
                 :thick 2
                 :legend (format nil "~A dB Elev, alpha = 0.2" dbelev))
  )

;; showing the effect of too-low an exponent (0.2)
(let ((dbelev 50))
  (plt:fplot 'plt '(20 100) #'identity
             :clear t
             :color :darkgreen
             :thick 2
             :title "Effect of Incorrect Exponent"
             :xtitle "Presentation Level [dBHL]"
             :ytitle "Perceived Level [dBHL]"
             :legend "Normal Hearing"
             )
  (plt:paramplot 'plt '(20 100)
                 (lambda (pdb)
                   (+ pdb (a-gain pdb dbelev)))
                 #'identity
                 :color :blue
                 :thick 2
                 :legend "50 dB Elev, correct alpha")
  (plt:paramplot 'plt '(20 100)
                 (lambda (pdb)
                   (+ pdb (a-gain pdb dbelev :alpha 0.4)))
                 #'identity
                 :color :red
                 :thick 2
                 :legend "50 dB Elev, alpha = 0.2")
  (plt:plot 'plt '(40 60.3) '(40 40)
            :color :blue
            :alpha 0.5
            :thick 2
            :legend "Perceptual Excess = 13 dB")
  (plt:plot 'plt '(60.3 60.3) '(40 52.7)
            :color :blue
            :alpha 0.5
            :thick 2)
  (plt:plot 'plt '(61 68) '(45 45)
            :color :black)
  (plt:draw-text 'plt "Perceptual Excess" '(69 45))
  (plt:plot 'plt '(57 68) '(39 34)
             :color :black)
  (plt:draw-text 'plt "Overcorrection" '(69 33))
  )

;; Show true perception versus Phons presentation
(let ((dbelev 30))
  (plt:fplot 'plt '(0 100) #'sdb
             :thick 2
             :xtitle "Presented Sound [Phon]"
             :ytitle "Perceived Loudness [dB Sone]"
             :title  "Normal and Impaired Hearing"
             :clear t
             :watermarkfn nil)
  (plt:fplot 'plt '(0 100) #'asdb
             :alpha 0.5)
  (labels ((plot (dbelev)
             (plt:fplot 'plt '(0 100) (lambda (p)
                                        (sdb (dbuncorr dbelev p)))
                        :color :red)))
    (plot 30)
    (plot 50)
    (plot 70)
    (plot 90)))

(let ((dbelev 70))
  (labels ((plot (dbelev &rest args)
             (apply #'plt:fplot 'plt2 '(0.01 100) (lambda (p)
                                                    (dbcorr dbelev p))
                    args)))
    (plot 90
          :ytitle "dB correction"
          :xtitle "Phon"
          :thick 2
          :clear t)
    (plot 70)
    (plot 50)
    (plot 30)))
|#
#|
|#

;; -----------------------------------------------------------------
;; operator algebra -- functions all take dbphons as input
;;

(defun gain (lvldb)
  (um:rcurry #'+ lvldb))

(defun atten (lvldb)
  (um:rcurry #'- lvldb))

(defun recruitment (dbelev)
  ;; the apparent dbphons corresponding to input dbphons
  ;; with threshold elevation dbelev
  (um:curry #'dbuncorr dbelev))

(defun correction-gain (dbelev)
  ;; gain needed to correct hearing with
  ;; threshold elevation dbelev
  (um:curry #'dbcorr dbelev))

(defun corrected-input (dbelev)
  ;; input level needed to hear dbphons as dbphons
  ;; given threshold elevation dbelev
  (um:combine #'+ #'identity (correction-gain dbelev)))

#|
(defun arecruitment (dbelev &key (alpha alpha))
  ;; the apparent dbphons corresponding to input dbphons
  ;; with threshold elevation dbelev
  (lambda (pdb)
    (adbuncorr dbelev pdb :alpha alpha)))

(let ((win (plt:wset 'xx :clear t)))
  (plt:fplot win '(30 100) #'identity :color :darkgreen :thick 2) ;; normal hearing
  (plt:fplot win '(30 100) (recruitment 50) :color :red :thick 2)
  (plt:fplot win '(30 100) (arecruitment 50 :alpha 0.4) :color :magenta :thick 2)  
  )


;; demonstrate recruitment hearing
(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(20 100) #'identity :color :darkgreen
             :fullgrid nil
             :watermarkfn nil
             :title "Sensioneural Hearing Loss"
             :xtitle "Input Presentation Level [dB]"
             :ytitle "Perceived Level [dB]"
             :legend "Normal Hearing"
             )
  (plt:fplot win '(20 100) (recruitment 60) :color :red
             :legend "Recruitment Hearing"))

(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(20 100) #'identity :color :darkgreen)
  (plt:fplot win '(20 100) (recruitment 50) :color :red)
  (plt:plot win '(20 100) '(30 110) :color :orange) ;; EQ
  (plt:plot win '(20 100) '(60 100) :color :magenta)) ;; Linear Compression

;; NL Compression vs EQ
(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(20 100) #'identity :color :darkgreen
             :axis-values nil
             :thick 2
             :watermarkfn nil
             :title "Comparison of Nonlinear Compression with EQ"
             :xtitle "Input Level [dB]"
             :ytitle "Output Level [dB]"
             :legend "Normal Hearing")
  ;; (plt:fplot win '(20 100) (recruitment 70) :color :blue :thick 2)
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (dbcorr 70 p)))
             :color :blue :thick 2
             :legend "Anti-Recruitment Compression")
  (plt:plot win '(20 100) '(30 110) :color :red :thick 2
            :legend "Static EQ")) ;; EQ


;; NL Compression vs Linear Compression
(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(20 100) #'identity :color :darkgreen
             :axis-values nil
             :thick 2
             :watermarkfn nil
             :title "Comparison of Linear & Nonlinear Compression"
             :xtitle "Input Level [dB]"
             :ytitle "Output Level [dB]"
             :legend "Normal Hearing")
  ;; (plt:fplot win '(20 100) (recruitment 70) :color :blue :thick 2)
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (dbcorr 70 p)))
             :color :blue :thick 2
             :legend "Anti-Recruitment Compresssion")
  (plt:plot win '(20 100) '(73 100) :color :red :thick 2
            :legend "3:1 Linear Compression")) ;; Linear Compression

;; -----------------------------------------------------------
;; NL Compression vs NY Compression

(let* ((win (plt:wset 1 :clear t))
      (dbel    50)
      (thresh  (- dbel 20))
      (gmakeup 21.4)
      (ratio   4.6))
  (plt:fplot win '(20 100) #'identity :color :darkgreen
             :axis-values nil
             :thick 2
             :title "Comparison of Linear & NYC Compression"
             :xtitle "Input Level [dB]"
             :ytitle "Output Level [dB]")
  ;; (plt:fplot win '(20 100) (recruitment 70) :color :blue :thick 2)
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (dbcorr dbel p)))
             :color :blue :thick 2)
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (if (< p thresh)
                                      gmakeup
                                    (+ gmakeup (* (- p thresh) (- (/ ratio) 1))))))
             :color :orange
             :thick 2)
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (db20 (+ 1 (ampl20 (if (< p thresh)
                                                         gmakeup
                                                       (+ gmakeup (* (- p thresh) (- (/ ratio) 1)))))))))
             :color :red
             :thick 2))

;; show what happens when compression ratio = 3 and we add
;; the linear compression to the dry signal
(let* ((win (plt:wset 1 :clear t))
      (dbel    70)
      (thresh  40)
      (ratio   3)
      (gmakeup 30))
  (plt:fplot win '(20 100) #'identity :color :darkgreen
             :axis-values nil
             :thick 2
             :title "Comparison of NYC & Nonlinear Compression"
             :xtitle "Input Level [dB]"
             :ytitle "Output Level [dB]")
  ;; (plt:fplot win '(20 100) (recruitment 70) :color :blue :thick 2)
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (dbcorr dbel p)))
             :color :blue :thick 2)
  #|
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (if (< p thresh)
                                      gmakeup
                                    (+ gmakeup (* (- p thresh) (- (/ ratio) 1))))))
             :color :orange
             :thick 2)
  |#
  (plt:fplot win '(20 100) (lambda (p)
                             (+ p (db20 (+ 1 (ampl20 (if (< p thresh)
                                                         gmakeup
                                                       (+ gmakeup (* (- p thresh) (- (/ ratio) 1)))))))))
             :color :red
             :thick 2))

(defun gparcmpr (ratio gmakeup dpthr)
  ;; dB gain for parallel compression at dpthr above threshold
  (db20 (+ 0 (ampl20 (+ gmakeup (* dpthr (- (/ ratio) 1)))))))


(defun pargain-at-db (dbel pdbfs)
  ;; compute parallel compressor gain at pdb for given elevation threshold dbel.
  ;; (assumes 100 dBSPL = 0 dBFS)
  ;; e.g., dbel = 60, pdbfs = -50 => gpar = 13.57 dB
  ;; i.e., pdbfs at -50 is 50 dBSPL which is below threshold elevation of 60 dB.
  ;; So, to sound like 50 dBSPL we need parallel gain of 13.57 dB
  (let* ((ratio   4.6)
         (gmakeup 21.4)
         (thresh  (- dbel 20 100)))
    (+ gmakeup (* (- pdbfs thresh) (- (/ ratio) 1)))))
;; -----------------------------------------------------------
 
;; compute makeup gain, thrsholds for parallel linear compression
(let* ((aud '((1 28.6) ;; fkhz, db thrshold elevation DBM
              (2 44.7)
              (4 60.6)
              (8 73.3)))
       (db0  (mapcar (um:compose #'ath #'first) aud))
       (crest 6)
       (bwfact 6)
       (g0    21.4)
       (r     4.6)
       (gdb   (+ g0 (* crest (- (/ r) 1))))
       (spl0  77)
       (dbfs0 -20)
       (thrs (mapcar (lambda (pair db0)
                       (let ((dbthr (second pair)))
                         (+ (+ db0 dbthr -20 crest bwfact) (- dbfs0 spl0))))
                     aud db0)))
  (list g0 thrs))
       
(defun dbhc-vt (intcpt fkhz)
  (let ((fbk (cbr fkhz)))
    (* 3.575 (- fbk intcpt))))

(defun integrated-power (rolloff f1khz f2khz)
  (labels ((fn (f)
             (/ (* f (expt 10 (* 0.1 rolloff (log f 2))))
                (+ 1 (* rolloff 0.1 (log 10 2))))))
    (db10 (- (fn f2khz) (fn f1khz)))))

(mapcar (lambda (grp)
          (destructuring-bind (fckhz fkhz< fkhz>) grp
            (let* ((fbk  (cbr fckhz))
                   (pwr  (integrated-power -6 fkhz< fkhz>))
                   (bpwr (integrated-power -6
                                           (inverse-cbr (- fbk 0.5))
                                           (inverse-cbr (+ fbk 0.5)))))
              (list fckhz (- pwr bpwr)))))  ;; -> fctr khz  bwfactor
        '((1 0.75 1.5) ;; fctr khz  f< khz  f> khz
          (2 1.5  3)
          (4 3    6)
          (8 6    20)))

(defun bwfact (rolloff f0)
  (let* ((f2 (* f0 (sqrt 2))) ;; (* f0 1.5))
         (f1 (/ f0 (sqrt 2))) ;; (* f0 0.75))
         (fb  (cbr f0))
         (f1b (inverse-cbr (- fb 0.5)))
         (f2b (inverse-cbr (+ fb 0.5))))
    (- (integrated-power rolloff f1 f2)
       (integrated-power rolloff f1b f2b))))

(defun bwfact8 (rolloff)
  (let* ((f2 20)
         (f1  (* 4 (sqrt 2))) ;; 6)
         (fb  (cbr 8))
         (f1b (inverse-cbr (- fb 0.5)))
         (f2b (inverse-cbr (+ fb 0.5))))
    (- (integrated-power rolloff f1 f2)
       (integrated-power rolloff f1b f2b))))

(let ()
  (plt:fplot 'plt '(0 -24) (um:rcurry #'bwfact 1)
             :clear t
             :thick 2
             :yrange '(0 10)
             :title "BW Factor"
             :xtitle "Rolloff [dB/oct]"
             :ytitle "Factor [dB]"
             :legend "1 kHz")
  (plt:fplot 'plt '(0 -24) (um:rcurry #'bwfact 2) :color :red
             :thick 2
             :legend "2 kHz")
  (plt:fplot 'plt '(0 -24) (um:rcurry #'bwfact 4) :color :blue
             :thick 2
             :legend "4 kHz")
  (plt:fplot 'plt '(0 -24) #'bwfact8 :color :orange
             :thick 2
             :legend "8 kHz")
  (plt:draw-ellipse 'plt -6 6 2 1
                    :filled nil
                    :border-color :gray50
                    :border-thick 0.1))

(defun cthrs (vt &key (dbfs -20) (dbspl 77))
  ;; 0 <= vt <= 10
  (let* ((intcpt (- (* 3 vt) 10))
         (slope  3.575)) ;; dB/Bark
    (mapcar (lambda (f)
              (list f
                    (+ (* slope (- (cbr f) intcpt))
                       (ath f)
                       (- dbfs dbspl)
                       6  ;; bwfact
                       6  ;; crest fact
                       -20 ;; nom thrsh below elev thrsh
                       ) ))
            '(1 2 4 8))))
    

|#

;; ------------------------------------------------------
(defun hyper-recruitment (dbelev dbh)
  ;; apparent dbphons for threshold elevation dbelev
  ;; and hyper-recruitment dbh
  (um:compose
   (recruitment (+ dbelev dbh))
   (gain dbh)))

#|
;; demonstrate hyper-recruitment hearing
(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(20 100) #'identity :color :black)
  (plt:fplot win '(20 100) (hyper-recruitment 50 10)))
|#

(defun hyper-correction-gain (dbelev dbh)
  ;; gain needed to correct hyper-recruitment hearing
  ;; with threshold elevation dbelev and hyper-recruitment dbh
  (um:compose
   (atten dbh)
   (correction-gain (+ dbelev dbh))))

(defun hyper-corrected-input (dbelev dbh)
  ;; input level required to hear dbphons as dbphons
  ;; given threshold elevation dbelev and hyper-recruitment dbh
  (um:combine #'+ #'identity (hyper-correction-gain dbelev dbh)))

;; ----------------------------------------------
;;
(defun decruitment (dbelev dbd)
  ;; apparent dbphons for threshold elevation dbelev
  ;; and decruitment of dbd
  (um:compose
   (recruitment (- dbelev dbd))
   (atten dbd)))

#|
;; demonstrate decruitment hearing
(let ((win (plt:wset 1 :clear t)))
  (plt:fplot win '(20 100) #'identity :color :black)
  (plt:fplot win '(20 100) (decruitment 50 10)))
|#

(defun decruitment-correction-gain (dbelev dbd)
  ;; gain needed to correct decruitment hearing
  ;; with threshold elevation dbelev and decruitment dbd
  (um:compose
   (gain dbd)
   (correction-gain (- dbelev dbd))))

(defun decruitment-corrected-input (dbelev dbd)
  ;; input level required to hear dbphons as dbphons
  ;; given threshold elevation dbelev and decruitment dbd
  (um:combine #'+ #'identity (decruitment-correction-gain dbelev dbd)))

#|
;; demonstrate recruitment, decruitment, and hyper-recruitment hearing
(let ((win (plt:wset 1 :clear t)))
  
  (plt:fplot win '(0 100) #'identity :color :black :thick 1 ;; normal hearing
                     :fullgrid nil
                     :watermarkfn nil
                     :title "Variations on Recruitment"
                     :xtitle "Presentation Level [dB]"
                     :ytitle "Perceived Level [dB]")
             
  (plt:fplot win '(0 100) (recruitment 50) :color :darkgreen :thick 2)
  (plt:fplot win '(0 100) (decruitment 50 10)       :color :blue  :thick 2)
  (plt:fplot win '(0 100) (hyper-recruitment 50 10) :color :magenta :thick 2)
  (plt:draw-text win "Normal Recruitment" '(:data 88 90) :anchor :e :color :darkgreen)
  (plt:draw-text win "HyperRecruitment"   '(:data 70 80) :anchor :e :color :magenta)
  (plt:draw-text win "Decruitment"        '(:data 80 62) :anchor :w :color :blue))

;; demonstrate inadequacy of normal recruitment correction applied to decruitment
(let ((win (plt:wset 2 :clear t)))
  (plt:fplot win '(0 100) #'identity :color :black :thick 1) ;; normal hearing
  (plt:fplot win '(0 100) (decruitment 50 10) :color :magenta)
  (plt:fplot win '(0 100) (um:compose (decruitment 50 10)
                                      (corrected-input 50))
             :color :blue  :thick 2)
  (plt:fplot win '(0 100) (um:compose (hyper-recruitment 50 10)
                                      (corrected-input 50))
             :color :red  :thick 2)
  (plt:fplot win '(0 100) (um:compose (hyper-recruitment 50 10)
                                      (decruitment-corrected-input 50 10))
             :color :magenta  :thick 2)
  (plt:fplot win '(0 100) (um:compose (decruitment 50 10)
                                      (hyper-corrected-input 50 10))
             :color :magenta  :thick 2)
  (plt:fplot win '(0 100) (um:compose (recruitment 50)
                                      (corrected-input 55))
             :color :orange :thick 2)
  )

;; -------------------------------------------------------------------------------
;; Decruitment and HyperRecruitment can be corrected to a very close approximation
;; by applying pre-EQ to the signal before correction and presentation to the ear.
;;
;; When threshold audiology is tested, the results will appear to indicate a recruitment
;; for the measured threshold elevation. But proper correction for some degree of departure
;; from normal recruitment actually needs to use an adjusted recruitment threshold.
;;
;; Not knowing this adjusted threshold, we can approximate by pre-gain or pre-attenuation
;; since Corr(dBthr,dBsig+dBgain) = Corr(dBthr-dBgain,dBsig) to a very good approximation
;; over louder levels. The approximation breaks down at the lower near-threshold levels.
;;

;; show departure of correction approximation
;; In order to approximate a higher threshold we pre-attenuate the signal
;; In order to approximate a lower threshold we pre-boost the signal
(let ((win (plt:wset 2 :clear t)))
  (plt:fplot win '(0 100) (um:combine '-
                                      (um:compose (correction-gain 60)
                                                  (gain 10))
                                      (correction-gain 50))
             :thick 2
             :color :blue
             :legend "PreGain for Lower Thresh"
             :title  "PreEQ vs Proper Tresh Correction Gains"
             :xtitle "Input Level [dB]"
             :ytitle "PreEQ Error [dB]"
             :yrange '(-1 0.5))
  (plt:fplot win '(0 100) (um:combine '-
                                      (um:compose (correction-gain 40)
                                                  (atten 10))
                                      (correction-gain 50))
             :thick 2
             :color :red
             :legend "PreAtten for Higher Thresh")
  )

;; show how we can correct decruitment by using normal recruitment correction
;; after a boost gain applied to signal
(let ((win (plt:wset 2 :clear t)))
  (plt:fplot win '(0 100) #'identity
             :color :black
             :thick 1 ;; normal hearing
             :title "Decruitment Correction"
             :xtitle "Input Level [dBHL]"
             :ytitle "Perceived Level [dBHL]")
  ;; show impaired loudness response for decruitment
  (plt:fplot win '(0 100) (decruitment 50 10)
             :color :magenta
             :thick 2)
  (plt:fplot win '(0 100) (recruitment 50)
             :color :magenta
             :thick 2
             :alpha 0.5)
  ;; show inadequacy of normal recruitment correction
  ;; not loud enough almost everywhere
  (plt:fplot win '(0 100) (um:compose (decruitment 50 10)
                                      (corrected-input 50))
             :color :blue
             :thick 2
             :alpha 0.5)
  #|
  ;; show "proper" decruitment correction
  (plt:fplot win '(0 100) (um:compose (decruitment 50 10)
                                      (decruitment-corrected-input 50 10))
             :color :red  :thick 2)
  |#
  ;; show our pre-gain approximation for dercruitment correction
  ;; starts to break down at very low levels
  (plt:fplot win '(0 100) (um:compose (decruitment 50 10)
                                      (um:compose
                                       (corrected-input 50)
                                       (gain 10)))
             :color :darkgreen  :thick 2)
  #|
  ;; what happens if we just correct as though there were no decruitment,
  ;; then boosted the corrected signal before presentation to the ear?
  (plt:fplot win '(0 100) (um:compose (decruitment 50 10)
                                      (um:compose
                                       (gain 10)
                                       (corrected-input 50)
                                       ))
             :color :orange  :thick 2)
  |#
  )

;; show excess in decruitment correction approximation
(let ((win (plt:wset 2 :clear t)))
  (plt:fplot win '(0 100) (um:combine '-
                                      (um:compose
                                       (corrected-input 50)
                                       (gain 10))
                                      (decruitment-corrected-input 50 10))
             :thick 2
             :title "Decruitment Approx Excess"
             :xtitle "Input Level [dB]"
             :ytitle "Excess [dB]"
             ))

;; show how we can correct hyper-ruitment by using normal recruitment correction
;; after an attenuation applied to signal
(let ((win (plt:wset 2 :clear t)))
  (plt:fplot win '(0 100) #'identity
             :color :black
             :thick 1 ;; normal hearing
             :title "HyperRecruitment Correction"
             :xtitle "Input Level [dBHL]"
             :ytitle "Perceived Level [dBHL]")
  ;; show impaired loudness response for hyper-recruitment
  (plt:fplot win '(0 100) (hyper-recruitment 50 10)
             :color :magenta
             :thick 2)
  (plt:fplot win '(0 100) (recruitment 50)
             :color :magenta
             :thick 2
             :alpha 0.5)
  ;; show inadequacy of normal recruitment correction
  ;; too loud almost everywhere
  (plt:fplot win '(0 100) (um:compose (hyper-recruitment 50 10)
                                      (corrected-input 50))
             :color :blue
             :thick 2
             :alpha 0.5)
  #|
  ;; show "proper" hyper-recruitment correction
  (plt:fplot win '(0 100) (um:compose (hyper-recruitment 50 10)
                                      (hyper-corrected-input 50 10))
             :color :red  :thick 2)
  |#
  ;; show our pre-atten approximation for hyper-recruitment correction
  ;; starts to break down at very low levels
  (plt:fplot win '(0 100) (um:compose (hyper-recruitment 50 10)
                                      (um:compose
                                       (corrected-input 50)
                                       (atten 10)))
             :color :darkgreen  :thick 2)

  #|
  ;; what happens if we just correct as though there were no hyper-recruitment,
  ;; then attenuated the corrected signal before presentation to the ear?
  (plt:fplot win '(0 100) (um:compose (hyper-recruitment 50 10)
                                      (um:compose
                                       (atten 10)
                                       (corrected-input 50)
                                       ))
             :color :orange  :thick 2)
  |#
  )


;; show deficiency in hyper-recruitment correction approximation
(let ((win (plt:wset 2 :clear t)))
  (plt:fplot win '(0 100) (um:combine '-
                                      (um:compose
                                       (corrected-input 50)
                                       (atten 10))
                                      (hyper-corrected-input 50 10))
             :thick 2
             :title "HyperRecruitment Approx Deficiency"
             :xtitle "Input Level [dB]"
             :ytitle "Deficiency [dB]"
             ))

(let ((win (plt:wset 'xx :clear t)))
  (plt:fplot win '(30 100) #'identity :color :black :thick 1) ;; normal hearing
  (plt:fplot win '(30 100) (recruitment 50) :color :darkgreen :thick 2)
  (plt:fplot win '(30 100) (um:compose (recruitment 60) (gain 10)) :color :blue :thick 2)
  (plt:fplot win '(30 100) (um:compose (atten 10) (recruitment 60) (gain 10)) :color :orange :thick 2))

(let ((win (plt:wset 'x2 :clear t)))
  (plt:fplot win '(50 100) (lambda (p)
                            (- (funcall (um:compose (atten 10) (recruitment 60) (gain 10)) p)
                               (funcall (recruitment 50) p)))
             :color :darkgreen
             :thick 2))

;; show that HC(p+dp; p0) = HC(p; p0-dp)
(let ((win (plt:wset 'x2 :clear t))
      (dbel 50))
  ;; normal recruitment
  (plt:fplot win '(30 100) (correction-gain dbel)
             :color :darkgreen
             :thick 2
             :yrange '(-10 30)
             :title  (format nil "Correction Gains for ~D dB Elevation" dbel)
             :xtitle "Input Loudness [dB]"
             :ytitle "Correction Gain [dB]")
  ;; hyper-recruitment correction
  (plt:fplot win '(0 100) (um:compose (correction-gain (- dbel 10)) (atten 10))
             :color :orange
             :thick 2)
  (plt:fplot win '(0 100) (um:compose (correction-gain (+ dbel 10)) (gain 10))
             :color :magenta
             :thick 2))

;; show normal, hyper, and de-cruitment gain corrections
(let ((win (plt:wset 'xx :clear t))
      (dbel 50))
  ;; normal recruitment
  (plt:fplot win '(30 100) (correction-gain dbel)
             :color :darkgreen
             :thick 2
             :yrange '(-10 30)
             :legend "Normal Recruitment"
             :legend-x '(:frac 0.65)
             :title  (format nil "Correction Gains for ~D dB Elevation" dbel)
             :xtitle "Input Loudness [dB]"
             :ytitle "Correction Gain [dB]")
  ;; hyper-recruitment correction
  (plt:fplot win '(0 100) (um:compose (atten 10) (correction-gain dbel) (atten 10))
             :color :orange
             :thick 2
             :legend "Hyper-Recruitment")
  (plt:fplot win '(0 100) (um:compose (gain 10) (correction-gain dbel) (gain 10))
             :color :magenta
             :thick 2
             :legend "Decruitment"))

;; show error in HC(p+dp; p0) = HC(p; p0-dp)
(let ((win (plt:wset 'x2 :clear t)))
  (plt:fplot win '(50 100) (lambda (p)
                            (- (funcall (um:compose (correction-gain 50) (gain 10)) p)
                               (funcall (correction-gain 40) p)))
             :yrange '(-0.1 0.1)
             :color :darkgreen
             :thick 2))
|#

;; -------------------------------------------------------
;; (needs "Remez")

(defun nder (f x)
  ;; numerical derivative of f at x
  (funcall (deriv f 0.01d0) x))

(defun correction-slope-fn (dbelev)
  ;; the sensitivity of apparent loudness with
  ;; threshold elevation dbelev to errors in the correction level
  ;; at input level dbphons
  (lambda (dbphons)
    (/ (+ 1.0d0 (nder (correction-gain dbelev) dbphons)))))


(defun wt-unrestricted-rational-minimax (ord dbelev)
  (let* ((domain '(20 100)))
    (list :dbelev dbelev
          :domain domain
          :fit    (weighted-minimax-rational-approx
                   (correction-gain dbelev)
                   domain
                   (correction-slope-fn dbelev)
                   ord))
    ))

(defun get-wt-unrestricted-rational-minimax-fits ()
  (loop for dbelev from 5.0d0 to 90.0d0 by 5.0d0 collect
        (progn
          (print dbelev)
          (handler-case
              (wt-unrestricted-rational-minimax '(3 2) dbelev)
            (bad-error-curve ()
              (print "Retry with order (3,1)")
              (wt-unrestricted-rational-minimax '(3 1) dbelev)))
          )))

(defun show-peak-errors (fits)
  (let* ((dbelevs (mapcar (lambda (fit)
                            (getf fit :dbelev))
                          fits))
         (errs    (mapcar (lambda (fit)
                            (getf (getf fit :fit) :err))
                          fits)))
    (let ((win (plt:wset 1 :clear t)))
      (plt:plot win dbelevs errs
                :title  "MiniMax (3,2) Fits"
                :xtitle "Threshold Elevation [dB]"
                :ytitle "Max Abs Error [dB]"
                :symbol :circle
                :plot-joined t))
    ))

#|
;; keep retrying until all elevations succeed with (3,2) fits
(progn
  (setf fits (get-wt-unrestricted-rational-minimax-fits))
  (show-peak-errors fits)
  fits)
|#

;; -------------------------------------------------------------
;; Fletcher-Munson Corrections

(defun sq (x)
  (* x x))

(defun ath-raw (fkhz)
  (+ (/ 3.64d0 (expt fkhz 0.8d0))
                (* -6.5d0 (exp (* -0.6d0 (sq (- fkhz 3.3d0)))))
                (* 0.001d0 (expt fkhz 4.0d0))))

;; absolute threshold of hearing
(defun ath (fkhz)
  (min 120 (- (ath-raw (max 0.02 fkhz)) (* 1 (ath-raw 1.0d0)))))

#|
(defun inverse-hcmapping (dbhl fkhz)
  ;; map phons to dBSPL
  (let* ((th (ath fkhz)))
    (+ th (* dbhl (/ (max (- 120 th) 0) 120))
       )))
|#

(defun inverse-hcmapping (dbhl fkhz)
  ;; map phons to dBSPL
  (let* ((th  (ath fkhz))
         (thx (if (minusp th)
                          (* th (+ 1 (/ dbhl 120.0 2)))
                        (* th (- 1 (/ dbhl 120.0))))))
    (+ dbhl thx)))

(defun amusic (fkhz)
  ;; approximate spectral envelope for loud music
  (cond ((< fkhz 0.05)
         (- 70 (* 10 (/ (- (log 0.05 10) (log fkhz 10))
                        (- (log 0.05 10) (log 0.02 10))))))

        ((< fkhz 1)  70)
        (t  (- 70 (* 6 (/ (log fkhz) (log 2)))))
        ))

(defun phons-of-sones (s)
  (find-root #'(lambda (phons)
                 (- s (sones phons)))
             30))

(defun sones-contours ()
  (let ((domain  '(0.02 20))
        (win (plt:wset 'sones)))
    (plt:with-delayed-update (win)
      (plt:clear win)
      (plt:fplot win domain #'(lambda (fkhz)
                                (inverse-hcmapping 0 fkhz))
                 :xlog t
                 :color :gray75
                 :thick 2
                 :xrange '(0.019 20.1)
                 :yrange '(-10 110)
                 :title  "Sones Contours"
                 :xtitle "Frequency [kHz]"
                 :ytitle "Sound Intensity [dBSPL]"
                 :cright1 "Copyright (c) by SpectroDynamics, LLC")
      (let ((s0 (sones 0)))
        (loop for ix from 1 to 16 do
              (let ((s (* s0 (expt 2 ix))))
                (plt:fplot win domain
                           #'(lambda (fkhz)
                               (inverse-hcmapping (phons-of-sones s) fkhz))
                           :color :gray75)
                
                ))
        (plt:fplot win domain
                   #'(lambda (fkhz)
                       (inverse-hcmapping (phons-of-sones (* s0 512)) fkhz))
                   :color :darkgreen
                   :alpha 0.3)
        (loop for ix from 6 to 16 do
              (let* ((x 0.15)
                     (y (inverse-hcmapping (phons-of-sones (* s0 (expt 2 ix))) x)))
                (plt:draw-text win (format nil "~A" ;; (* 3 (- ix 9))
                                           (expt 2 (- ix 9)))
                               (list x y)
                               :anchor :ctr
                               :color :gray55
                               :font-size 10
                               :transparent nil))))
      (plt:fplot win '(0.02 20) #'amusic
                 :thick 2
                 :color :orange
                 :alpha (/ 96 256))
      )))

(defun phons-contours ()
  (let ((win (plt:wset 'phons)))
    (plt:with-delayed-update (win)
      (plt:axes win
                :xlog t
                :xrange '(0.019 20.1)
                :yrange '(-10 110)
                :title  "Phons Contours"
                :xtitle "Frequency [kHz]"
                :ytitle "Sound Intensity [dBSPL]"
                :watermarkfn nil)
      (loop for lvl from 0 to 100 by 10 do
            (plt:fplot win '(0.02 20) #'(lambda (fkhz)
                                          (inverse-hcmapping lvl fkhz))
                       :color :gray75)
            (plt:draw-text win (format nil "~A" lvl)
                           `(:data 0.19 ,(inverse-hcmapping lvl 0.2))
                           :anchor :se
                           :font-size 10
                           :color :gray50))
      (plt:fplot win '(0.02 20) #'amusic
                 :thick 2
                 :color :orange
                 :alpha (/ 96 256))
      )))

(defparameter *dm-hearing*
  ;; fkHz left right
  '((0.25 10 5)
    (0.5  15 15)
    (0.75 25 15)
    (1    30 15)
    (1.5  55 50)
    (2    60 55)
    (3    55 55)
    (4    60 55)
    (6    70 75)
    (8    70 65)))

#|
(defun hcmapping (dbspl fkhz)
  ;; map dbspl to phons
  (let ((th (ath fkhz)))
    (/ (* (- dbspl th) 120) (max (- 120 th) 1))))
|#

(defun hcmapping (dbspl fkhz)
  ;; map dbspl to phons
  (let* ((th  (ath fkhz)))
    (/ (max 0.0 (- dbspl th))
       (if (minusp th)
           (+ 1d0 (/ th 240d0))
         (- 1d0 (/ th 120d0)))
       )))

(defun dbel-to-phon (dbel fkhz)
  ;; convert elevation in dB above ATH into Phon
  (let ((th (ath fkhz)))
    (max 0 (min 120 (+ 120 (* (/ 120 (- 120 th)) (+ dbel th -120)))))))

(defun dbspl-to-phon (dbspl fkhz)
  (dbel-to-phon (- dbspl (ath fkhz)) fkhz))

(defun gdbphon-to-gdbspl (gdbphon fkhz)
  (let ((th (ath fkhz)))
    (* gdbphon (/ (- 120 th) 120))))

(defun dbspl-contours ()
  (let ((win (plt:wset 'dbspl)))
    (plt:with-delayed-update (win)
      (plt:axes win
                :xlog t
                :xrange '(0.019 20.1)
                :yrange '(0 110)
                :title  "dBSPL Contours"
                :xtitle "Frequency [kHz]"
                :ytitle "Sound Level [Phon]")
      (loop for lvl from 0 to 100 by 10 do
            (plt:fplot win '(0.02 20) (um:curry 'hcmapping lvl)
                       :color :gray75)
            (plt:draw-text win (format nil "~A" lvl)
                           `(:data 1.5 ,(hcmapping lvl 1.5))
                           :color :gray50
                           :font-size 10
                           :anchor :e))
      (plt:fplot win '(0.02 20) #'(lambda (fkhz)
                                    (hcmapping (amusic fkhz) fkhz))
                 :thick 2
                 :color :orange
                 :alpha (/ 96 256))
      )))

;; Bark frequency scale
;; Band centers:
(defparameter *bark-centers*
  '(50   150   250   350   450 
      570   700   840  1000  1170 
      1370  1600  1850  2150  2500 
      2900  3400  4000  4800  5800 
      7000  8500  10500  13500))

;; Band edges:
(defparameter *bark-edges*
  '(  0   100   200   300   400 
     510   630   770   920  1080 
    1270  1480  1720  2000  2320 
    2700  3150  3700  4400  5300 
    6400  7700  9500  12000  15500))

;; these can be extended by appending [20500, 27000] to accommodate sampling rates up to 54 kHz.
;; --------------------------------------
;; Allpass bilinear transformation vs sample rate Fsamp
;; Optimum allpass coeff = 0.8517*sqrt(atan(0.06583*Fs))-0.1916, for fs > 1 kHz

(defun bark-bw (fkhz)
  ;; bandwidth in Hz for Bark centered at fkhz
  (+ 94 (* 71 (expt fkhz 1.5))))

;; Critical band rate z (in bark)
(defun cbr (fkhz)
  ;; Traunmuller (1990)
  ;; this gives negative barks at low frequencies!
  ;; ... but it agrees within 0.05 bark with measured band edges
  ;; e.g., (cbr 0.4) = 4.01 bark
  (if (zerop fkhz)
      0
    (max 0 (- (/ 26.81d0 (+ 1d0 (/ 1.960d0 fkhz))) 0.53d0))))

(defun cbr-corr (fkhz)
  (let ((z (cbr fkhz)))
    (cond ((< z 2)     (+ z (* 0.15 (- 2 z))))
          ((> z 20.1)  (+ z (* 0.22 (- z 20.1))))
          (t           z))))

(defun inverse-cbr (zbark)
  ;; conversion back to khz
  (/ 1.960d0 (- (/ 26.81d0 (+ zbark 0.53d0)) 1d0)))
  
(defun inverse-cbr-corr (zbark &optional (fstart 1))
  (if (zerop zbark)
      0
    (find-root (lambda (f)
                 (- zbark (cbr-corr f)))
               fstart)))

;; Using cbr with newton root finding we have:
(defparameter *quarter-bark-freqs*
  ;; (zbark fkHz)
  '(
    ( 0.00    0.0d0)
    ( 0.25   0.03 #| 0.0587322331325569D0 |#) ;; adj by hand 

    (0.5 0.057592374143013D0)  ;; from inverse-cbr-corr
    (0.75 0.08063687783776609D0)
    (1.0 0.10421388231468207D0)
    (1.25 0.12834206045525162D0)
    (1.5 0.1530409685207375D0)
    (1.75 0.17833109901624567D0)
#|
    ( 0.00    0.0d0)
    ( 0.25   0.0587322331325569D0) 
    ( 0.5    0.07830876799049594D0) 
    ( 0.75   0.09826870539147667D0) 
    ( 1.0    0.11862342002759294D0) 
    ( 1.25   0.13938474102583573D0) 
    ( 1.5    0.16056497487605603D0) 
    ( 1.75   0.18217692975829525D0) 
|#
    ( 2.0    0.20423394137190606D0) 
    ( 2.25   0.2267499003765107D0) 
    ( 2.5    0.24973928156369857D0) 
    ( 2.75   0.273217174888302D0) 
    ( 3.0    0.29719931849908216D0) 
    ( 3.25   0.32170213392077546D0) 
    ( 3.5    0.34674276355277656D0) 
    ( 3.75   0.3723391106644038D0) 
    ( 4.0    0.39850988129730264D0) 
    ( 4.25   0.42527463345821504D0) 
    ( 4.5    0.45265381951175876D0) 
    ( 4.75   0.480668843486906D0) 
    ( 5.0    0.5093421151553674D0) 
    ( 5.25   0.5386971098578761D0) 
    ( 5.5    0.5687584326252008D0) 
    ( 5.75   0.599551886971761D0) 
    ( 6.0    0.6311045487714507D0) 
    ( 6.25   0.6634448456630049D0) 
    ( 6.5    0.6966026424757601D0) 
    ( 6.75   0.7306093332160245D0) 
    ( 7.0    0.7654979402098747D0) 
    ( 7.25   0.8013032210606076D0) 
    ( 7.5    0.8380617837260729D0) 
    ( 7.75   0.8758122134205719D0) 
    ( 8.0    0.9145952037887154D0) 
    ( 8.25   0.9544537068730857D0) 
    ( 8.5    0.995433090203232D0)  ;; 1 kHz
    ( 8.75   1.0375813094124307D0) 
    ( 9.0    1.0809490951066863D0) 
    ( 9.25   1.1255901569112797D0) 
    ( 9.5    1.1715614050465517D0) 
    ( 9.75   1.218923191843503D0) 
    (10.0    1.2677395824132148D0) 
    (10.25   1.3180786282734008D0) 
    (10.5    1.3700127009350925D0) 
    (10.75   1.4236188300246849D0) 
    (11.0    1.4789790863742583D0) 
    (11.25   1.5361810012833002D0) 
    (11.5    1.5953180283303912D0) 
    (11.75   1.6564900528592148D0) 
    (12.0    1.7198039549804398D0) 
    (12.25   1.7853742327626472D0) 
    (12.5    1.853323693251881D0) 
    (12.75   1.9237842200845197D0) 
    (13.0    1.9968976277707072D0) 
    (13.25   2.0728166174723724D0) 
    (13.5    2.151705832175385D0) 
    (13.75   2.233743060234433D0) 
    (14.0    2.319120566308852D0) 
    (14.25   2.4080465971576435D0) 
    (14.5    2.5007470775328766D0) 
    (14.75   2.5974675267010077D0) 
    (15.0    2.698475229820725D0) 
    (15.25   2.8040617046123626D0) 
    (15.5    2.914545511255776D0) 
    (15.75   3.030275462558315D0) 
    (16.0    3.1516343025319827D0) 
    (16.25   3.279042935105536D0) 
    (16.5    3.412965301409301D0) 
    (16.75   3.5539140247270105D0) 
    (17.0    3.7024569678713295D0) 
    (17.25   3.8592248797913165D0) 
    (17.5    4.024920348473617D0) 
    (17.75   4.200328328047534D0) 
    (18.0    4.386328587785254D0) 
    (18.25   4.583910425454156D0) 
    (18.5    4.794190324669622D0) 
    (18.75   5.018433032597315D0) 
    (19.0    5.258077025408785D0) 
    (19.25   5.514765398929299D0) 
    (19.5    5.790383593499363D0) 
    (19.75   6.087105784576593D0) 
    (20.0    6.407452353892727D0)

    (20.25 6.715445589658973D0) ;; from inverse-cbr-corr
    (20.5 7.019224712459529D0)
    (20.75 7.34505004453391D0)
    (21.0 7.6954118951150035D0)
    (21.25 8.07319031759314D0)
    (21.5 8.481734459995645D0)
    (21.75 8.924962124940873D0)
    (22.0 9.407485810181854D0)
    (22.25 9.934773827972009D0)
    (22.5 10.513358447361114D0)
    (22.75 11.151107889045664D0)
    (23.0 11.857586257171399D0)
    (23.25 12.644536465876959D0)
    (23.5 13.52653815613389D0)
    (23.75 14.521919334121922D0)
    (24.0 15.654043726055843D0)
    (24.25 16.953167826277568D0)
    (24.5 18.459185183283065D0)
    (24.75 20.225795252875916D0)
    (25.0 22.327041276910866D0))
#|
    (20.25   6.75436165691085D0) 
    (20.5    7.131280415033132D0) 
    (20.75   7.542278626419609D0) 
    (21.0    7.992197121690805D0) 
    (21.25   8.4868391215342D0) 
    (21.5    9.03322193313294D0) 
    (21.75   9.639911887396785D0) 
    (22.0   10.317476836315157D0) 
    (22.25  11.079106915364007D0) 
    (22.5   11.941481713835595D0) 
    (22.75  12.926005917092544D0) 
    (23.0   14.060610029025389D0) 
    (23.25  15.382442540748834D0) 
    (23.5   16.94201470466513D0) 
    (23.75  18.80980273763015D0) 
    (24.0   21.08719339285816D0) 
    (24.25  23.925517706931075D0) 
    (24.5   27.561124130986684D0) 
    (24.75  32.38483722074965D0) 
    (25.0   39.09281299829528D0) 
|#
    )

;; Critical bandwdith in kHz
(defun cbw (zbark)
  (/ 52.548d0 (+ 690.39 (* zbark (- zbark 52.56d0)))))

;; Mel scale - equal pitch intervals
;; 1000 mel = 1 kHz @ 40 dBHL
(defun mel (fkhz)
  (* 1127.01048d0 (log (+ 1.0d0 (/ fkhz 0.7d0)))))

;; !! Standard Musical Octave ??
(defun std-octave (fkhz)
  (log (/ fkhz 0.12709d0) 2))

(defun bark-ctr (fkhz)
  (+ (* 13.1478 (atan (* 0.821587 fkhz)))
     (* 3.87049 (atan (sq (/ fkhz 8.07315))))))

(defun bark-to-khz (b)
  (find-root #'(lambda (f)
                 (- b (bark-ctr f)))
             1))

;; --------------------------------------------

(defun erb (fkhz)
  ;; equiv rectangular bandwidth
  (+ 24.7 (* 108 fkhz)))

(defun erbs (fkhz)
  ;; ERB scale = number of ERB's below fkHz
  (* 21.4 (log (+ (* 4.37 fkhz) 1) 10)))

(defun erbs-to-khz (erbs)
  (/ (- (expt 10 (/ erbs 21.4)) 1) 4.37))

;; --------------------------------------------

(defun poly-eval (x &rest coffs)
  (if (endp coffs)
      0d0
    (+ (first coffs) (* x (apply #'poly-eval x (rest coffs))))
    ))

(defun rational-minimax-log10-bark-to-log10-khz (log10-zbark)
  ;; approx good to about 5%
  (let* ((rred (- (* 2d0 (/ (- log10-zbark (log (bark-ctr 0.02d0) 10d0))
                            (- (log (bark-ctr 20d0) 10d0) (log (bark-ctr 0.02d0) 10d0)))
                     )
                  1d0))
         (num (poly-eval rred
                         -0.6824869506755976D0 1.6308954905762994D0 -0.736856210803108D0))
         (den (poly-eval rred
                         1.0d0 -0.8268908540945528D0 -0.007640825732406986D0)))
    (/ num den)))

(defun rational-minimax-bark-to-khz (zbark)
  (expt 10.d0 (rational-minimax-log10-bark-to-log10-khz (log zbark 10d0))))

#|
(let ((win (plt:wset 'bark-ctr :clear t)))
  (plt:fplot win '(0.02 20) #'bark-ctr :xlog t :thick 2
             :title "Freq to Bark"
             :xtitle "Freq [kHz]"
             :ytitle "Bark Freq [Bk]"
             :legend "BARK-CTR")
  (plt:fplot win '(0.02 20) #'cbr :color :orange :thick 2
             :legend "CBR")
  (plt:plot  win (mapcar (um:rcurry #'/ 1000) *bark-centers*)
             (loop for ix from 0.5
                   for bc in *bark-centers*
                   collect ix)
             :symbol :circle
             :color :red
             :legend "Std Bark Centers"))
|#

(defun bark-ctr-x (fkhz)
  (cbr fkhz)) ;; or (bark-ctr fkhz)

(defun bark-to-khz-x (zbark)
  (inverse-cbr zbark)) ;; or (bark-to-khz zbark)

(defun phons-bark-contours ()
  (let ((domain '(0 24)))
    (let ((win (plt:wset 'bark-phons)))
      (plt:with-delayed-update (win)
        (plt:axes win
                  :xlog nil
                  :xrange '(0 24)
                  :yrange '(-10 110)
                  :title  "Phons Contours"
                  :xtitle "Frequency [Bark]"
                  :ytitle "Sound Intensity [dBSPL]"
                  :watermarkfn nil)
        (loop for f in '(0.25 0.5 0.75 1 1.5 2 3 4 6 8 12 16) do
              (let ((b (cbr f)))
                (plt:plot win (list b b) (list -20 120)
                          :color :violet
                          :alpha 0.5)
                (plt:draw-text win (format nil "~A" f)
                               `(:data ,b 110)
                               :color :gray50
                               :font-size 10
                               :anchor :s)))
        (loop for lvl from 0 to 100 by 10 do
              (plt:fplot win domain #'(lambda (b)
                                        (inverse-hcmapping lvl (bark-to-khz-x b)))
                         :color :gray75)
              (plt:draw-text win (format nil "~A" lvl)
                             `(:data ,(bark-ctr-x 0.19) ,(inverse-hcmapping lvl 0.2))
                             :color :gray50
                             :font-size 10
                             :anchor :se))
        (plt:fplot win domain
                   #'(lambda (b)
                       (amusic (bark-to-khz-x b)))
                   :thick 2
                   :color :orange
                   :alpha (/ 96 256))
        ))))

#|
;; ---------------------------------------
;; Audiology in dBHL

(progn
  (let* ((frqs (mapcar #'first *dm-hearing*))
         (lf   (mapcar #'second *dm-hearing*))
         (rt   (mapcar #'third  *dm-hearing*))
         (xfrqs (loop for f in frqs
                      for ix from 1
                      collect ix)))
    (plt:plot 'aud xfrqs lf
              :title  "Audiology"
              :xtitle "Frequency [kHz]"
              :ytitle "Threshold Elevation [dBHL]"
              ;; :xrange '(0.1 20)
              :yrange '(100 -10)
              ;; :xlog   t
              :color :blue
              :symbol :circle
              :plot-joined t
              :clear t
              :legend "Left"
              :x-values '("-" "0.25" "0.5" "0.75" "1" "1.5" "2" "3" "4" "6" "8")
              :axis-values t) 
    (plt:plot 'aud xfrqs rt
              :color :red
              :symbol :circle
              :plot-joined t
              :legend "Right")
    
    (when nil
    (let* ((frqs (mapcar (um:rcurry #'log 10) frqs))
           (rtlf (mapcar (lambda (lf rt)
                           (* 0.5 (+ lf rt)))
                         lf rt))
           (mnfrq (vm:mean xfrqs))
           (mnaud (vm:mean rtlf))
           (cfrqs (mapcar (lambda (f)
                            (- f mnfrq))
                          xfrqs))
           (caud  (mapcar (lambda (a)
                            (- a mnaud))
                          rtlf))
           (slope (/ (reduce #'+ (mapcar #'* caud cfrqs))
                     (reduce #'+ (mapcar #'* cfrqs cfrqs))))
           )
      (print slope)
      (plt:fplot 'aud '(0 10)
                 (lambda (f)
                   (+ mnaud (* (- f mnfrq) slope))))))))


;; do it again in Bark frequency space
(progn
  (let* ((frqs (mapcar #'cbr #| #'bark-ctr |# (mapcar #'first *dm-hearing*)))
         (lf   (mapcar #'second *dm-hearing*))
         (rt   (mapcar #'third  *dm-hearing*)))
    (plt:plot 'audb frqs lf
              :title  "Audiology"
              :xtitle "Frequency [Bark]"
              :ytitle "Threshold Elevation [dBHL]"
              :xrange '(1 23)
              :yrange '(100 -10)
              :color :blue
              :symbol :circle
              :plot-joined t
              :clear t
              :legend "Left")
    (plt:plot 'audb frqs rt
              :color :red
              :symbol :circle
              :plot-joined t
              :legend "Right")
    

    (let* ((rtlf (mapcar (lambda (lf rt)
                           (* 0.5 (+ lf rt)))
                         lf rt))
           (mnfrq (vm:mean frqs))
           (mnaud (vm:mean rtlf))
           (cfrqs (mapcar (lambda (f)
                            (- f mnfrq))
                          frqs))
           (caud  (mapcar (lambda (a)
                            (- a mnaud))
                          rtlf))
           (slope (/ (reduce #'+ (mapcar #'* caud cfrqs))
                     (reduce #'+ (mapcar #'* cfrqs cfrqs))))
           )
      (print slope)
      (plt:fplot 'audb '(1 23)
                 (lambda (f)
                   (+ mnaud (* (- f mnfrq) slope))))
      (plt:draw-text 'audb
                     (format nil "Slope = ~,3F" slope)
                     '(:data 15.5 18))
      (plt:draw-text 'audb
                     (format nil "Intcpt = ~,3F" (- mnfrq (/ mnaud slope)))
                     '(:data 15.5 24))
      )))

;; ------------------------------------------------------------------------------

(defvar *clm-hearing* '((0.25  20 20)
                      (0.5   20 35)
                      (0.75  35 35)
                      (1     40 40)
                      (1.5   50 50)
                      (2     50 60)
                      (3     65 65)
                      (4     75 65)
                      (6     70 70)
                      (8     80 80)))

(defvar *ms-hearing* '(;; (0.125 3 5)
                      (0.25  4 1)
                      (0.5   2 2)
                      (0.75  11 11)
                      (1     12 11)
                      (1.5   18 11)
                      (2     18 19)
                      (3     27 34)
                      (4     44 44)
                      (6     42 51)
                      (8     47 62)
                      ;; (12    65 68)
                      ))

;; --------------------------------------
;; Average slope among DM, CLM, MS = 3.33
;; --------------------------------------

;; Mark Schleune's audiology
(let ((hearing *ms-hearing*))
  (let* ((frqs (mapcar #'cbr #| #'bark-ctr |# (mapcar #'first hearing)))
         (lf   (mapcar #'second hearing))
         (rt   (mapcar #'third  hearing)))
    (plt:plot 'audb frqs lf
              :title  "Audiology"
              :xtitle "Frequency [Bark]"
              :ytitle "Threshold Elevation [dBHL]"
              :xrange '(1 23)
              :yrange '(100 -10)
              :color :blue
              :symbol :circle
              :plot-joined t
              :clear t
              :legend "Left")
    (plt:plot 'audb frqs rt
              :color :red
              :symbol :circle
              :plot-joined t
              :legend "Right")
    

    (let* ((rtlf (mapcar (lambda (lf rt)
                           (* 0.5 (+ lf rt)))
                         lf rt))

           #||#
           (pos   (position-if (um:rcurry #'> 20) rtlf))
           (rtlf  (subseq rtlf pos))
           (frqs  (subseq frqs pos))
           #||#
           
           (mnfrq (vm:mean frqs))
           (mnaud (vm:mean rtlf))
           (cfrqs (mapcar (lambda (f)
                            (- f mnfrq))
                          frqs))
           (caud  (mapcar (lambda (a)
                            (- a mnaud))
                          rtlf))
           (slope (/ (reduce #'+ (mapcar #'* caud cfrqs))
                     (reduce #'+ (mapcar #'* cfrqs cfrqs))))
           )
      (print slope)
      (plt:fplot 'audb '(1 23)
                 (lambda (f)
                   (+ mnaud (* (- f mnfrq) slope))))
      (plt:draw-text 'audb
                     (format nil "Slope = ~,3F" slope)
                     '(:data 15.5 18))
      (plt:draw-text 'audb
                     (format nil "Intcpt = ~,3F" (- mnfrq (/ mnaud slope)))
                     '(:data 15.5 24))
      )))

;; find intercept, forcing slope to 3.33
;; Mark Schleunes Audiology 01/13
(let ((hearing *ms-hearing*))
  (let* ((frqs (map 'vector #'cbr #| #'bark-ctr |# (mapcar #'first hearing)))
         (lf   (mapcar #'second hearing))
         (rt   (mapcar #'third  hearing)))
    (plt:plot-bars 'audb frqs
                   (list (mapcar (um:rcurry #'- 5) lf)
                         (mapcar (um:rcurry #'+ 5) lf))
                   :title  "Mark's Audiology"
                   :xtitle "Frequency [Bark]"
                   :ytitle "Threshold Elevation [dBHL]"
                   :xrange '(1 23)
                   :yrange '(100 -10)
                   :color :blue
                   :clear  t)
    (plt:plot 'audb frqs lf
              :color :blue
              :symbol :circle
              :plot-joined t
              :legend "Left")
    (plt:plot-bars 'audb frqs (list
                               (mapcar (um:rcurry #'- 5) rt)
                               (mapcar (um:rcurry #'+ 5) rt))
              :color :red)
    (plt:plot 'audb frqs rt
              :color :red
              :symbol :circle
              :plot-joined t
              :legend "Right")
    

    (let* ((rtlf (map 'vector (lambda (lf rt)
                                (* 0.5 (+ lf rt)))
                      lf rt))
           (slope #|3.575|# 3.33)
           (pos   (position-if (um:rcurry #'> 14) rtlf))
           (rtlf  (subseq rtlf pos))
           (fs    (subseq frqs pos)))
      (multiple-value-bind (b sigma) (linfit:regress-fixed-slope fs rtlf 1 slope)
        (format t "~%sigma = ~A" sigma)
        (let* ((best-intercept (- (/ b slope)))
               ;; marks preferred vTuning = 5.73
               ;; gives an intercept of 7.19
               (marks-intercept (- (* 3 5.73) 10)))
          (plt:fplot 'audb '(1 23)
                     (lambda (f)
                       (* slope (- f marks-intercept)))
                     :thick 2
                     :legend "Mark's Est")
          (plt:fplot 'audb '(1 23)
                     (lambda (f)
                       (* slope (- f best-intercept)))
                     :thick 2
                     :color :orange
                     :legend "Best Fit")
          #||#
          (plt:draw-text 'audb
                         (format nil "Slope = ~,3F" slope)
                     '(:data 16 18))
          (plt:draw-text 'audb
                         (format nil "Intcpt = ~,3F (~,3F)"
                                 marks-intercept
                                 best-intercept)
                         '(:data 16 24))
          #||#
          )))))

;; find intercept, forcing slope to 3.33
;; Clark McClain Audiology 02/13
(let ((hearing *clm-hearing*))
  (let* ((frqs (map 'vector #'cbr #| #'bark-ctr |# (mapcar #'first hearing)))
         (lf   (mapcar #'second hearing))
         (rt   (mapcar #'third  hearing)))
    (plt:plot-bars 'audb frqs
                   (list (mapcar (um:rcurry #'- 5) lf)
                         (mapcar (um:rcurry #'+ 5) lf))
                   :title  "Clark's Audiology"
                   :xtitle "Frequency [Bark]"
                   :ytitle "Threshold Elevation [dBHL]"
                   :xrange '(1 23)
                   :yrange '(100 -10)
                   :color :blue
                   :clear  t)
    (plt:plot 'audb frqs lf
              :color :blue
              :symbol :circle
              :plot-joined t
              :legend "Left")
    (plt:plot-bars 'audb frqs (list
                               (mapcar (um:rcurry #'- 5) rt)
                               (mapcar (um:rcurry #'+ 5) rt))
              :color :red)
    (plt:plot 'audb frqs rt
              :color :red
              :symbol :circle
              :plot-joined t
              :legend "Right")
    

    (let* ((rtlf (map 'vector (lambda (lf rt)
                                (* 0.5 (+ lf rt)))
                      lf rt))
           (slope #|3.575|# 3.33)
           (pos   (position-if (um:rcurry #'> 14) rtlf))
           (rtlf  (subseq rtlf pos))
           (fs    (subseq frqs pos)))
      (multiple-value-bind (b sigma) (linfit:regress-fixed-slope fs rtlf 1 slope)
        (format t "~%sigma = ~A" sigma)
        (let* ((best-intercept (- (/ b slope)))
               ;; marks preferred vTuning = 5.73
               ;; gives an intercept of 7.19
               (marks-intercept (- (* 3 5.73) 10)))
          #|
          (plt:fplot 'audb '(1 23)
                     (lambda (f)
                       (* slope (- f marks-intercept)))
                     :thick 2
                     :legend "Mark's Est")
          |#
          (plt:fplot 'audb '(1 23)
                     (lambda (f)
                       (* slope (- f best-intercept)))
                     :thick 2
                     :color :orange
                     :legend "Best Fit")
          #||#
          (plt:draw-text 'audb
                         (format nil "Slope = ~,3F" slope)
                     '(:data 16 18))
          (plt:draw-text 'audb
                         (format nil "Intcpt = ~,3F"
                                 best-intercept)
                         '(:data 16 24))
          #||#
          )))))

;; find intercept, forcing slope to 3.5
;; David McClain Audiology
(let* ((subject 'ms)
       (hearing (ecase subject
                  (dm  *dm-hearing*)
                  (ms  *ms-hearing*)
                  (clm *clm-hearing*)))
       (title (ecase subject
                (dm  "David's Audiology")
                (clm "Clark's Audiology")
                (ms  "Mark's Audiology")))
       (frqs (map 'vector #'cbr #| #'bark-ctr |# (mapcar #'first hearing)))
       (lf   (mapcar #'second hearing))
       (rt   (mapcar #'third  hearing)))
  (plt:plot-bars 'audb frqs
                 (list (mapcar (um:rcurry #'- 5) lf)
                       (mapcar (um:rcurry #'+ 5) lf))
                 :title  title
                 :xtitle "Frequency [Bark]"
                 :ytitle "Threshold Elevation [dBHL]"
                 :xrange '(1 23)
                 :yrange '(100 -10)
                 :color :blue
                 :clear  t)
  (plt:plot 'audb frqs lf
            :color :blue
            :symbol :circle
            :plot-joined t
            :legend "Left")
  (plt:plot-bars 'audb frqs (list
                             (mapcar (um:rcurry #'- 5) rt)
                             (mapcar (um:rcurry #'+ 5) rt))
                 :color :red)
  (plt:plot 'audb frqs rt
            :color :red
            :symbol :circle
            :plot-joined t
            :legend "Right")
  
  (let* ((rtlf (map 'vector (lambda (lf rt)
                              (* 0.5 (+ lf rt)))
                    lf rt))
         (slope #|3.575|# 3.5)
         (pos   (position-if (um:rcurry #'> 14) rtlf))
         (rtlf  (subseq rtlf pos))
         (fs    (subseq frqs pos)))
    (multiple-value-bind (b sigma) (linfit:regress-fixed-slope fs rtlf 1 slope)
      (format t "~%sigma = ~A" sigma)
      (let* ((best-intercept (- (/ b slope)))
             ;; DM/RAL 08/15: vTuning is now the elevation in dBHL at 4 kHz = 17.46 Bark
             )
        (plt:fplot 'audb '(1 23)
                   (lambda (f)
                     (* slope (- f best-intercept)))
                   :thick 2
                   :color :orange
                   :legend "Best Fit")
        #||#
        (plt:draw-text 'audb
                       (format nil "Slope = ~,3F" slope)
                       '(:data 16 18))
        (plt:draw-text 'audb
                       (format nil "vTun = ~,1F"
                               (* slope (- (cbr 4) best-intercept)))
                       '(:data 16 #|30|# 24))
        #||#
        ))))

;; ---------------------------------------
;; Hearing elevations plotted in dBSPL vs Log Hz

(let* ((frqs (mapcar #'first *dm-hearing*))
       (lf   (mapcar #'second *dm-hearing*))
       (rt   (mapcar #'third  *dm-hearing*))
       (lfa  (mapcar (lambda (f a)
                       (+ a (ath f)))
                     frqs lf))
       (rta  (mapcar (lambda (f a)
                       (+ a (ath f)))
                     frqs rt)))
  (phons-contours)
  (plt:plot 'phons frqs lfa
            :color :blue
            :symbol :circle
            :plot-joined t
            :legend "Left")
  (plt:plot 'phons frqs rta
            :color :red
            :symbol :circle
            :plot-joined t
            :legend "Right"))

(let* ((frqs   (mapcar #'first *dm-hearing*))
       (lfrqs  (mapcar (um:rcurry #'log 10) frqs))
       (rtlf   (mapcar (lambda (lf rt f)
                        (+ (ath f) (* 0.5 (+ lf rt))))
                      (mapcar #'second *dm-hearing*)
                      (mapcar #'third *dm-hearing*)
                      frqs))
       (mnfrq (vm:mean lfrqs))
       (mnaud (vm:mean rtlf))
       (cfrqs (mapcar (lambda (f)
                        (- f mnfrq))
                      lfrqs))
       (caud  (mapcar (lambda (a)
                        (- a mnaud))
                      rtlf))
       (slope (/ (reduce #'+ (mapcar #'* caud cfrqs))
                 (reduce #'+ (mapcar #'* cfrqs cfrqs))))
       )
  (print slope)
  (plt:fplot 'phons '(0.25 10)
             (lambda (f)
               (+ (ath f) mnaud (* (- (log f 10) mnfrq) slope)))))

;; ---------------------------------------
;; Hearing thresolds plotted in dB SPL vs Bark frequency

(let* ((frqs (mapcar #'first *dm-hearing*))
       (lf   (mapcar #'second *dm-hearing*))
       (rt   (mapcar #'third  *dm-hearing*))
       (lfa  (mapcar (lambda (f a)
                       (+ a (ath f)))
                     frqs lf))
       (rta  (mapcar (lambda (f a)
                       (+ a (ath f)))
                     frqs rt))
       (fbarks (mapcar #'cbr frqs)))
  (phons-bark-contours)
  (plt:plot 'bark-phons fbarks lfa
            :color :blue
            :symbol :circle
            :plot-joined t
            :legend "Left")
  (plt:plot 'bark-phons fbarks rta
            :color :red
            :symbol :circle
            :plot-joined t
            :legend "Right"))

(let* ((frqs     (mapcar #'first *dm-hearing*))
       (fbarks   (mapcar #'bark-ctr-x frqs))
       (lf       (mapcar #'second *dm-hearing*))
       (rt       (mapcar #'third  *dm-hearing*))
       (lfa      (mapcar (lambda (f a)
                           (+ a (ath f)))
                         frqs lf))
       (rta      (mapcar (lambda (f a)
                           (+ a (ath f)))
                         frqs rt))
       (rtlf     (mapcar (lambda (lf rt)
                           (* 0.5 (+ lf rt)))
                         lf rt))
       (mnfrq (vm:mean fbarks))
       (mnaud (vm:mean rtlf))
       (cfrqs (mapcar (lambda (f)
                        (- f mnfrq))
                      fbarks))
       (caud  (mapcar (lambda (a)
                        (- a mnaud))
                      rtlf))
       (slope (/ (reduce #'+ (mapcar #'* caud cfrqs))
                 (reduce #'+ (mapcar #'* cfrqs cfrqs))))
       )
  (print slope)
  (plt:fplot 'bark-phons '(2.5 21)
             (lambda (f)
               (+ (ath (inverse-cbr f)) mnaud (* (- f mnfrq) slope)))))
|#

#|
(loop for el in '(10 20 30 40 50 60 70 80 90 100) do
      (plt:fplot 'plt '(0.02 20)
                 (lambda (fkhz)
                   (dbel-to-phon el fkhz))
                 :clear (eql el 10)))
|#

(defun vtuning-elevation (vtuning zbark)
  (min 80 (* vtuning (max 0 (- zbark 2.5)))))

(defun cresc-gains (vtuning)
  (phons-bark-contours)
  (plt:fplot 'bark-phons '(0 24)
             (lambda (zbark)
               (inverse-hcmapping
                (vtuning-elevation vtuning zbark)
                (inverse-cbr zbark)))
             :color :violet
             :thick 2)
  (loop for ampl in '(0 -6 -12 -18 -24) do
        (progn
          (plt:fplot 'bark-phons (list (cbr 0.02) (cbr 18))
                     (lambda (zbark)
                       (+ ampl (amusic (bark-to-khz-x zbark))))
                     :color :orange
                     :thick 2
                     :alpha 0.5)
          (plt:fplot 'bark-phons (list (cbr 0.02) (cbr 18))
                     (lambda (zbark)
                       (let* ((fkhz  (inverse-cbr zbark))
                              (dbel  (vtuning-elevation vtuning zbark))
                              (pdb   (hcmapping (+ ampl (amusic fkhz)) fkhz))
                              (gdb   (gdbphon-to-gdbspl (min 24 (dbcorr dbel pdb)) fkhz)))
                         (+ (amusic fkhz) ampl gdb))))
          (plt:fplot 'bark-phons (list (cbr 0.02) (cbr 18))
                     (lambda (zbark)
                       (let* ((fkhz  (inverse-cbr zbark))
                              (dbel  (vtuning-elevation vtuning zbark))
                              (pdb   (hcmapping (+ ampl (amusic fkhz)) fkhz))
                              (gdb   (gdbphon-to-gdbspl (min 24 (dbcorr dbel pdb)) fkhz)))
                         gdb))
                     :color :blue
                     ))))
                              
                   
                 
        
    
(defun dbspl-bark-contours ()
  (let ((win (plt:wset 'dbspl-bark)))
    (plt:with-delayed-update (win)
      (plt:axes win
                :xlog nil
                :xrange '(0 24)
                :yrange '(0 115)
                :title  "dBSPL Contours"
                :xtitle "Frequency [Bark]"
                :ytitle "Sound Level [Phon]")
      (loop for lvl from 0 to 100 by 10 do
            (plt:paramplot win '(0.02 18.7)
                           #'bark-ctr-x
                           #'(lambda (fkhz)
                               (hcmapping lvl fkhz))
                           :color :gray75)
            (plt:draw-text win (format nil "~A" lvl)
                           `(:data 9.9 ,(hcmapping lvl (bark-to-khz-x 9.9)))
                           :color :gray50
                           :font-size 10
                           :anchor :e))
      (plt:paramplot win '(0.02 18.7)
                     #'bark-ctr-x
                     #'(lambda (fkhz)
                         (hcmapping (amusic fkhz) fkhz))
                     :thick 2
                     :color :orange
                     :alpha (/ 96 256))
      )))

(defun cresc-gains-phons (vtuning)
  (dbspl-bark-contours)
  (plt:fplot 'dbspl-bark '(0 24)
             (lambda (zbark)
               (vtuning-elevation vtuning zbark))
             :color :violet
             :thick 2)
  (loop for ampl in '(12 6 0 -6 -12 -18 -24) do
        (progn
          (plt:fplot 'dbspl-bark (list (cbr 0.02) (cbr 18))
                     (lambda (zbark)
                       (let ((fkhz (inverse-cbr zbark)))
                         (hcmapping (+ ampl (amusic (bark-to-khz-x zbark))) fkhz)))
                     :color :orange
                     :thick 2
                     :alpha 0.5)
          (plt:fplot 'dbspl-bark (list (cbr 0.02) (cbr 18))
                     (lambda (zbark)
                       (let* ((fkhz  (inverse-cbr zbark))
                              (dbel  (vtuning-elevation vtuning zbark))
                              (pdb   (hcmapping (+ ampl (amusic fkhz)) fkhz))
                              (gdb   (min 24 (dbcorr dbel pdb))))
                         (+ pdb gdb))))
          (plt:fplot 'dbspl-bark (list (cbr 0.02) (cbr 18))
                     (lambda (zbark)
                       (let* ((fkhz  (inverse-cbr zbark))
                              (dbel  (vtuning-elevation vtuning zbark))
                              (pdb   (hcmapping (+ ampl (amusic fkhz)) fkhz))
                              (gdb   (min 24 (dbcorr dbel pdb))))
                         gdb))
                     :color :blue
                     ))))
                              
                   
                 
(defparameter *dm-hearing-fbark-phons*
  (mapcar #'(lambda (entry)
              (destructuring-bind (fkhz ldb rdb) entry
                (list (bark-ctr-x fkhz)
                      (dbel-to-phon ldb fkhz)
                      (dbel-to-phon rdb fkhz))))
          *dm-hearing*))


  
#|
;; show apparent phons with a hearing impairment of 3.37 phons/bark
(defun app-phons (dbspl fkhz)
  (let* ((sl 3.37) ;; slope of phons/bark
         (phons (hcmapping dbspl fkhz))
         (dphel (* sl (bark-ctr-x fkhz))))
    (if (< phons dphel)
        -1
      (find-root
       #'(lambda (ph)
           (- (- (sones ph) (sones 0))
              (- (sones phons) (sones dphel))))
       phons))
    ))

(let ((win (plt:wset 0)))
  (plt:with-delayed-update (win)
    (plt:axes win
     :xrange '(0 25)
     :yrange '(0 110)
     :xtitle "Frequency [Bark]"
     :ytitle "Apparent Sound Level [Phon]")
    (loop for lvl from 0 to 100 by 10 do
          (plt:paramplot win
                         '(0.02 20)
                         #'bark-ctr-x
                         #'(lambda (fkhz)
                             (app-phons lvl fkhz))
                         :color :gray75)
          (plt:draw-text win
                         (format nil "~A" lvl)
                         `(:data 3 ,(app-phons lvl (bark-to-khz-x 3)))
                         :color :gray50
                         :font-size 10
                         :anchor :e))
    (plt:paramplot win
                   '(0.02 20)
                   #'bark-ctr-x
                   #'(lambda (fkhz)
                       (app-phons (amusic fkhz) fkhz))
                   :thick 2
                   :color :orange
                   :alpha (/ 96 256))
    ))

(let ((win (plt:wset 0)))
  (plt:with-delayed-update (win)
    (plt:axes win
     :xlog t
     :xrange '(0.019 20.1)
     :yrange '(0 110)
     :xtitle "Frequency [kHz]"
     :ytitle "Apparent Sound Level [Phon]")
    (loop for lvl from 0 to 100 by 10 do
          (plt:paramplot win
                         '(0.02 20)
                         #'identity
                         #'(lambda (fkhz)
                             (app-phons lvl fkhz))
                         :color :gray75)
          (plt:draw-text win
                         (format nil "~A" lvl)
                         `(:data 0.25 ,(app-phons lvl 0.25))
                         :color :gray50
                         :font-size 10
                         :anchor :e))
    (plt:paramplot win
                   '(0.02 20)
                   #'identity
                   #'(lambda (fkhz)
                       (app-phons (amusic fkhz) fkhz))
                   :thick 2
                   :color :orange
                   :alpha (/ 96 256))
    ))

;; Show the correction required for 60 dB threshold elevation
(let ((domain '(20 100))
      (win    (plt:wset 0)))
  (plt:with-delayed-update (win)
    (plt:clear win)
    (plt:fplot win
               domain #'identity
               :title "Hearing Recruitment"
               :xtitle "Sound Level [Phon]"
               :ytitle "Perceived Sound Level [Phon]"
               :yrange '(20 100)
               :thick  2
               :color :darkgreen)
    (plt:paramplot win
                   domain (corrected-input 60)
                   #'identity
                   :color :red
                   :thick 2)
    (plt:paramplot win
                   domain (corrected-input 40)
                   #'identity
                   :color :red
                   :thick 2)
    (plt:paramplot win
                   domain #'identity (corrected-input 60)
                   :color :magenta
                   :thick 2
                   :alpha (/ 64 256))
    (plt:paramplot win
                   domain #'identity (corrected-input 40)
                   :color :magenta
                   :thick 2
                   :alpha (/ 64 256))
    (plt:plot win
              '(45 45) '(-10 45)
              :color :blue
              :thick 2
              :alpha 0.5)
    (plt:plot win
              `(45 ,(funcall (corrected-input 60) 45)) '(45 45)
              :symbol :circle
              :plot-joined t
              :color :blue
              :thick 2
              :alpha 0.5)
    ))
                 
|#


#|
;; headphone excesses interpolated to Bark scale
;; Sennheiser HD-600 have excess of 3 dB at 250 Hz,
;; declining to -6 dB at 8 kHz. The decline is essential a
;; straight line in a log-frequency vs dB excess plot.
;; Excesses are measured relative to 0 dB at 1 kHz.
(let* ((freq '(0.25 0.5 0.75 1 1.5 2 3 4 6 8 12 16))
       (yfit (lambda (fkhz)
               (* 0.5
                  (round
                   (- 3 (/ (* 9 (- (bark-ctr-x fkhz) (bark-ctr-x 0.25)))
                           (- (bark-ctr-x 8) (bark-ctr-x 0.25)))
                      )
                   0.5)))
             ))
  (let ((*print-length* nil))
    (print (mapcar #'list freq (mapcar yfit freq)))))

;; Interestingly, using Bark interpolation produces the
;; estimate of 0 dB excess at 1 kHz, as it should be.
;; But using log scaling, instead, produces a -0.5 dB offset at 1 kHz...
;; Log scaling also overestimates the excess at higher frequencies by about 0.5-1 dB.

==>
((0.25 3.0)
 (0.5 2.0)
 (0.75 1.0)
 (1 0.0)
 (1.5 -1.5)
 (2 -2.0)
 (3 -3.5)
 (4 -4.0)
 (6 -5.0)
 (8 -6.0)
 (12 -7.0)
 (16 -7.5))

|#

#|
;; 4th order polynomial fits over the domain 30 dBSPL to 100 dBSPL for use by Chameleon DSP
;;
(defun wt-unrestricted-minimax (ord dbelev)
  (let* ((domain '(30 100)))
    (list :dbelev dbelev
          :domain domain
          :fit    (weighted-minimax-rational-approx
                   (correction-gain dbelev)
                   domain
                   (correction-slope-fn dbelev)
                   ord))
    ))

(defun get-wt-unrestricted-minimax-fits ()
  (loop for dbelev from 5.0d0 to 80.0d0 by 5.0d0 collect
        (progn
          (print dbelev)
          (handler-case
              (wt-unrestricted-minimax '(4 0) dbelev)
            (bad-error-curve ()
              (print "Retry with order (3,1)")
              (wt-unrestricted-minimax '(3 0) dbelev)))
          )))

(defun show-peak-errors (fits)
  (let* ((dbelevs (mapcar (lambda (fit)
                            (getf fit :dbelev))
                          fits))
         (errs    (mapcar (lambda (fit)
                            (getf (getf fit :fit) :err))
                          fits)))
    (let ((win (plt:wset 1 :clear t)))
      (plt:plot win dbelevs errs
                :title  "MiniMax 4th Order Fits"
                :xtitle "Threshold Elevation [dB]"
                :ytitle "Max Abs Error [dB]"
                :symbol :circle
                :plot-joined t))
    ))

(defun padd (lst1 lst2)
  (cond ((endp lst1) lst2)
        ((endp lst2) lst1)
        (t (cons (+ (first lst1) (first lst2))
                 (padd (rest lst1) (rest lst2))))
        ))

(defun pscale (sf lst)
  (mapcar (um:curry #'* sf) lst))

(defun pmul (lst1 lst2 &optional ans)
  ;; lst is a list of polynomial coefficients with ascending powers of x
  ;; e.g., (a b c) := a + b*x + c*x^2
  (cond ((endp lst1) ans)
        (t
         (let ((sum (padd (pscale (first lst1) lst2)
                          ans)))
           (cons (first sum)
                 (pmul (rest lst1) lst2 (rest sum)))
           ))
        ))



  
(defun print-c-fits (fits)
  (format t "~&Max Coff: ~A"
          (loop for fit in fits maximize
                (destructuring-bind (dbel-key dbel domain-key domain
                                              fit-key fit-params) fit
                  (destructuring-bind (nodes-key nodes coffs-key coffs &rest _) fit-params
                    (let* ((x0 (list 1d0))
                           (x1 (list 1d0 (/ 192.66d0 35.0d0)))
                           (x2 (pmul x1 x1))
                           (x3 (pmul x1 x2))
                           (x4 (pmul x1 x3))
                           (coffs (coerce coffs 'list))
                           (adj-coffs (reverse
                                       (reduce #'padd
                                               (mapcar (lambda (a xx)
                                                         (pscale (/ a 192.66d0 16) xx))
                                                       coffs (list x0 x1 x2 x3 x4))
                                               ))))
                      (format t "~&static float gfit~2,'0D[] = {-70, ~&~A~{, ~A~}};" (truncate dbel)
                              (first adj-coffs) (rest adj-coffs))
                      (reduce #'max (mapcar #'abs adj-coffs))
                      ))))))

;; keep retrying until all elevations succeed with (4,0) fits
(progn
  (setf fits (get-wt-unrestricted-minimax-fits))
  (show-peak-errors fits)
  (print-c-fits fits)
  nil)

(let ((win (plt:wset 'xx :clear t))
      (coffs #(-0.36651501757180893D0 -0.39580634215694116D0 -0.054896623048211204D0 -0.005051502782584217D0 7.440496514480444D-5)))
  (plt:fplot win '(30 100) (lambda (p)
                             (let ((x (/ (- p 100) 192.66d0)))
                               (* 192.66 16
                                  (+ (aref coffs 4)
                                     (* x (+ (aref coffs 3)
                                             (* x (+ (aref coffs 2)
                                                     (* x (+ (aref coffs 1)
                                                             (* x (aref coffs 0))))))))))
                               ))
             :thick 2
             :color :darkgreen))
                               
|#

#|
;; equal loudness filter
;; E(w) = w^4 (w^2 + 56.8e6)/((w^2 + 6.3e6)(w^2 + 0.38e9)(w^6 + 9.58e26))

(defun eqld-filt (fkhz)
  (let ((w (* 2 pi fkhz 1000)))
    (/ (* (expt w 4) (+ (expt w 2) 56.8e6))
       (* (+ (expt w 2) 6.3e6)
          (+ (expt w 2) 0.38e9)
          (+ (expt w 6) 9.58e26))
       )))

(plt:fplot 'eqld '(0.02 20) (um:compose #'db20 #'eqld-filt) :clear t :xlog t)
|#

#|
(setf faud-khz '(0.25 0.5 0.75 1 1.5 2 3 4 6 8 12 16 24))
(setf faud-bark (mapcar #'cbr faud-khz))
(setf fletch
      ;; (fbark, delta dB)
      '((0.0 36.607040586379284)
        (0.25 25.814769124709002)
        (0.5 19.583829825772057)
        (0.75 15.751949903789566)
        (1.0 13.218653604658158)
        (1.25 11.277773368732724)
        (1.5 9.801565253057348)
        (1.75 8.60981654133721)
        (2.0 7.640869063618908)
        (2.25 6.827122301769009)
        (2.5 6.13842123655003)
        (2.75 5.543539746801542)
        (3.0 5.026229587199509)
        (3.25 4.570006107720943)
        (3.5 4.165365617661266)
        (3.75 3.8099820040941808)
        (4.0 3.476225862325933)
        (4.25 3.132186368628126)
        (4.5 2.8081075573208665)
        (4.75 2.5228007499265175)
        (5.0 2.263677085859044)
        (5.25 2.030652712411154)
        (5.5 1.8086699922363607)
        (5.75 1.579706663572935)
        (6.0 1.3605851357049374)
        (6.25 1.1623414334382431)
        (6.5 0.9781637061923902)
        (6.75 0.8085381198676203)
        (7.0 0.6440728236675568)
        (7.25 0.4729006679199159)
        (7.5 0.30551700965779105)
        (7.75 0.14931750937196853)
        (8.0 0.0)
        (8.25 -0.14176545192580647)
        (8.5 -0.282878348047646)
        (8.75 -0.43097275588496586)
        (9.0 -0.5841109021291544)
        (9.25 -0.7447144598926485)
        (9.5 -0.9083913315623278)
        (9.75 -1.0677337712276938)
        (10.0 -1.231435106050327)
        (10.25 -1.4074538197702537)
        (10.5 -1.5950029204739167)
        (10.75 -1.7975342028364452)
        (11.0 -2.011036166528797)
        (11.25 -2.2271178253774053)
        (11.5 -2.4550644221857)
        (11.75 -2.701883055810086)
        (12.0 -2.9731665146687902)
        (12.25 -3.2843149216076455)
        (12.5 -3.6204081077064387)
        (12.75 -3.9571235079715166)
        (13.0 -4.314276365717376)
        (13.25 -4.713060489051813)
        (13.5 -5.136576775386715)
        (13.75 -5.566844493358531)
        (14.0 -6.008521126685437)
        (14.25 -6.466547026042377)
        (14.5 -6.908790891606429)
        (14.75 -7.294538042294093)
        (15.0 -7.650313300924353)
        (15.25 -7.987252787948506)
        (15.5 -8.22982331107552)
        (15.75 -8.341386768632091)
        (16.0 -8.329078298628687)
        (16.25 -8.156111768150055)
        (16.5 -7.808639602632553)
        (16.75 -7.351922989820849)
        (17.0 -6.756611521473388)
        (17.25 -5.932490757680554)
        (17.5 -5.026619609240782)
        (17.75 -4.226697604716605)
        (18.0 -3.4854992195670462)
        (18.25 -2.775356969968726)
        (18.5 -2.21098962422397)
        (18.75 -1.815539109425103)
        (19.0 -1.498294321222133)
        (19.25 -1.1937805508549388)
        (19.5 -0.8872677584433148)
        (19.75 -0.5716314158771496)
        (20.0 -0.20242491073685187)
        (20.25 0.2681145733449739)
        (20.5 0.8572431234816378)
        (20.75 1.5816609829554125)
        (21.0 2.507996710880434)
        (21.25 3.7779295063320433)
        (21.5 5.377062440067519)
        (21.75 7.103787922940413)
        (22.0 9.340813210043805)
        (22.25 12.963265260812605)
        (22.5 17.865538195724838)
        (22.75 23.32041455234321)
        (23.0 30.29976449959615)
        (23.25 40.500393500289576)
        (23.5 54.75728587300147)
        (23.75 59.9)
        (24.0 59.9)
        (24.25 59.9)
        (24.5 59.9)
        (24.75 59.9)
        (25.0 59.9)))

(plt:plot 'fletch (mapcar #'first fletch)
          (mapcar #'second fletch)
          :clear  t
          :yrange '(-10 130))

(setf fmc (mapcar (lambda (pair)
                    (list (first pair)
                          (/ 120.0 (- 120.0 (second pair)))))
                  fletch))
(plt:plot 'fletch (mapcar #'first fmc)
          (mapcar #'second fmc)
          :clear t
          )
(plt:plot 'invfletch (mapcar #'first fmc)
          (mapcar #'/ (mapcar #'second fmc))
          :clear t
          )

(defun interpolate (x xs ys)
  (cond ((<= x (aref xs 0))
         (aref ys 0))
        ((>= x (aref xs (1- (length xs))))
         (aref ys (1- (length xs))))
        (t
         (let* ((pos (position-if (um:curry #'< x) xs))
                (dx  (- x (aref xs (1- pos))))
                (delx (- (aref xs pos) (aref xs (1- pos))))
                (dely (- (aref ys pos) (aref ys (1- pos)))))
           (+ (aref ys (1- pos))
              (/ (* dx dely) delx))
           ))
        ))

(defun convert-to-phon (dbelev)
  (mapcar (lambda (pair)
            (let* ((ath (* 120.0 (- 1.0 (/ (second pair))))))
              (list (first pair)
                    (+ 120.0
                       (* (second pair)
                          (+ dbelev ath -120.0)))
                    )))
          fmc))

(setf dummy-aud-0 (convert-to-phon 0))

(progn
  (plt:plot 'fletch (mapcar #'first dummy-aud-0)
            (mapcar #'second dummy-aud-0)
            :clear t
            :yrange '(-10 120)
            )
  (loop for ix from 10 to 90 by 5 do
        (let ((aud (convert-to-phon ix)))
          (plt:plot 'fletch (mapcar #'first aud)
                    (mapcar #'second aud))))
  (loop for f in faud-bark do
        (plt:plot 'fletch (list f f) (list 0 120) :color :orange)))

(setf derl (loop for ix from 0.25 to 25 by 0.25 collect (list ix (min 2 (nder (um:rcurry #'hcmapping (bark-to-khz ix)) 77)))))

(dolist (der derl) (print der))


  
|#

(defun masking (zbark)
  (+ 15.81d0
     (* 7.5d0 (+ zbark 0.474d0))
     (* -17.5d0 (sqrt (+ 1d0 (sq (+ zbark 0.474d0)))))))

#|
(plt:fplot 'plt '(-5 10) #'masking :clear t
           :title "Noise Masking"
           :xtitle "Frequency [Bark]"
           :ytitle "Masking [dB]")

|#

(defun nn-to-khz (nn)
  (* 0.440 (expt 2 (/ (- nn 69) 12))))


#|
(progn
  (format t "~&FkHz  ZBark  dCorrdB")
  (dolist (f '(0.5 1 2 4 8))
    (format t "~&~3,1F   ~4,1F   ~5,1F"
            ;; using an assumed slope of 3.575 dB / ZBark
            f
            (cbr f)
            (* 3.575 (- (cbr f) (cbr 4))))))

;; HC threshold elevations relative to 4 kHz correction
;; based on VTuning with 3.575 dB/Bark slope
FkHz  ZBark  dCorrdB
0.5    4.9   -44.8
1.0    8.5   -31.9
2.0   13.0   -15.9
4.0   17.5     0.0  ;; 4 kHz ref correction
8.0   21.0    12.7

|#
