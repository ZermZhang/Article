ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	coef   function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function coef(x,y)
      real*8 eps,pp,factor,x,y,pi
	real*8 eps_1,eps_2
	common/epsilon/eps

      pp=1.8d0
      pi=dacos(-1.0d0)
      eps_1=eps
      eps_2=eps

          coef=(2.d0+pp*dsin(2.d0*pi*x/eps_1))/
     +           (2.d0+pp*dcos(2.d0*pi*y/eps_1))
     +    +(2.d0+dsin(2.d0*pi*y/eps_2))/
     +           (2.d0+pp*dsin(2.d0*pi*x/eps_2))
      return
      end
