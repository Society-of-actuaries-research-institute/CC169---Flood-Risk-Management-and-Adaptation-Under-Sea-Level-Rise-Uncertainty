
SUBROUTINE crank(mu, sigma, r0, g, nt, ns, Smin, Smax,T,x0,a0,b0,u0,x1,a1,b1,u1,delta,kappa, scale,xi,f)
USE value_fn
USE constants_NSWC

IMPLICIT NONE
INTEGER(4), INTENT(in):: nt,ns
REAL(dp), INTENT(in) :: mu, sigma, r0, g, Smin, Smax, T
real(dp), intent(in) :: u0,u1, x0,x1  
real(dp), intent(in) ::a0,b0, a1,b1 !a0 is before project, a1 is after project
real(dp), intent(in) :: delta,kappa, scale,xi

REAL(dp), INTENT(out), DIMENSION(nt,ns-1) :: f
REAL(dp):: a, b, c, a11, a12, an1,an2, ah, bh, ch
REAL(dp), DIMENSION(ns-1)::Rv, u, y, l
INTEGER :: i, j, i1
REAL(dp)  ::dt, dS, r

!==============================
  r = r0 - g
  dt = T/nt
  dS = (Smax-Smin)/ns
      a = dt*(mu*dS -sigma*sigma)/(4*dS*dS)
      b = 1 + dt*(r+sigma*sigma/(dS*dS))/2
      c = -dt*(mu*dS+sigma*sigma)/(4*dS*dS)
      ah = -a
      bh = 1 - dt*(r+sigma*sigma/(dS*dS))/2
      ch = -c
!boundary conditions
      a11=b 
      a12=c
      an1=a-c
      an2=b+2*c
  !-----------
  DO j = 1, (ns-1)
    f(nt, j) = 0 !Final period
  END DO

  DO i1=1, nt
    i=nt-i1
    !Calculate rhs vector
    DO j=1, (ns-1)
      IF (j==1) THEN 
        Rv(1) = bh*f(i+1, 1)+ch*f(i+1, 2)    
      ELSE IF (j==(ns-1)) THEN 
        Rv(ns-1) = ah*f(i+1, ns-2)+bh*f(i+1,ns-1)+ch*(2*f(i+1,ns-1)-f(i+1, ns-2))
      ELSE                
        Rv(j) = ah*f(i+1, j-1)+ bh*f(i+1,j)+ch*f(i+1,j+1)
      END IF
    END DO
    DO j=1, (ns-1) 
      Rv(j)=Rv(j)+ (1+delta)*kappa*(Dfn(x0,a0,b0,u0,Smin + j*dS,scale,xi)-Dfn(x1,a1,b1,u1,Smin + j*dS,scale,xi))*dt
    END DO
! Solve for option value f
    call minv(a11,a12, a,b,c, an1,an2, ns, Rv, f(i,:))
  END DO
END SUBROUTINE crank


