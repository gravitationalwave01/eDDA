function besselk1(x)
USE ddprecision, only : WP

REAL(WP) besselk1,x

!Returns the modified Bessel function K1(x) for positive real x.
REAL(WP) besseli1
integer i

DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1, q2,q3,q4,q5,q6,q7,y !Accumulate polynomials in double precision.
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
p1 = 1.0d0
p2 = 0.15443144d0
p3 = -0.67278579d0
p4 = -0.18156897d0
p5 = -0.1919402d-1
p6 = -0.110404d-2
p7 = -0.4686d-4

q1 = 1.25331414d0
q2 = 0.23498619d0
q3 = -0.3655620d-1
q4 = 0.1504268d-1
q5 = -0.780353d-2
q6 = 0.325614d-2
q7 = -0.68245d-3

if (x.le.2.0) then !Polynomial t.
	y=x*x/4.0
	besselk1=(log(x/2.0)*besseli1(x))+(1.0/x)*(p1+y*(p2+ &
	y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
	y=2.0/x
	besselk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+ & 
	y*(q4+y*(q5+y*(q6+y*q7))))))
endif

return
End function besselk1
