function besseli1(x)
USE ddprecision, only : WP

REAL(WP) besseli1,x

!Returns the modified Bessel function I1(x) for any real x.
REAL(WP) ax
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y !Accumulate polynomials in double precision.
SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
p1 = 0.5d0 
p2 = 0.87890594d0 
p3 = 0.51498869d0
p4 = 0.15084934d0 
p5 = 0.2658733d-1 
p6 = 0.301532d-2
p7 = 0.32411d-3

q1 = 0.39894228d0 
q2 = -0.3988024d-1
q3 = -0.362018d-2 
q4 = 0.163801d-2
q5 = -0.1031555d-1
q6 = 0.2282967d-1
q7 = -0.2895312d-1
q8 = 0.1787654d-1
q9 = -0.420059d-2

if (abs(x).lt.3.75) then !Polynomial fit.
	y=(x/3.75)**2
	besseli1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
	ax=abs(x)
	y=3.75/ax
	besseli1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+ &
	y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))

        if(x.lt.0.) besseli1=-besseli1
endif

return
END function besseli1
