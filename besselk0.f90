function besselk0(x)

USE ddprecision,only : WP
REAL(WP) besselk0,x
integer i

!Returns the modied Bessel function K0(x) for positive real x.
REAL(WP) besseli0
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y !Accumulate polynomials in double precision.

SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
p1 = -0.57721566d0 
p2 = 0.42278420d0
p3 = 0.23069756d0
p4 = 0.3488590d-1 
p5 = 0.262698d-2
p6 = 0.10750d-3
p7 = 0.74d-5

q1 = 1.25331414d0
q2 = -0.7832358d-1
q3 = 0.2189568d-1
q4 = -0.1062446d-1
q5 = 0.587872d-2
q6 = -0.251540d-2
q7 = 0.53208d-3

if (x.le.2.0) then !Polynomial fit.
	y=x*x/4.0
	besselk0=(-log(x/2.0)*besseli0(x))+(p1+y*(p2+y*(p3+ &
	y*(p4+y*(p5+y*(p6+y*p7))))))
else
	y=(2.0/x)
	besselk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+ &
	y*(q4+y*(q5+y*(q6+y*q7))))))
endif

return

end function besselk0
