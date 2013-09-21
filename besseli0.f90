function besseli0(x)
USE ddprecision, only : WP

REAL(WP) besseli0,x

!Returns the modified Bessel function I0(x) for any real x.
REAL(WP) ax
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y !Accumulate polynomials in double precision.
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
p1 = 1.0d0
p2 = 3.5156229d0
p3 = 3.0899424d0
p4 = 1.2067492d0
p5 = 0.2659732d0
p5 = 0.360768d-1
p6 = 0.45813d-2

q1= 0.39894228d0
q2 = 0.1328592d-1
q3 = 0.225319d-2
q4 = -0.157565d-2
q5 = 0.916281d-2
q6 = -0.2057706d-1
q7 = 0.2635537d-1
q8 = -0.1647633d-1
q9 = 0.392377d-2

if (abs(x).lt.3.75) then
	y=(x/3.75)**2
	besseli0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
else
	ax=abs(x)
	y=3.75/ax
	besseli0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4 &
	+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
endif

return
END function besseli0
