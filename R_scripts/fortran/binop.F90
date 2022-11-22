SUBROUTINE binop(mu, sigma, gamma, S0, r, nt, T, N, V,E, f, B) 
USE value_fn
USE constants_NSWC

IMPLICIT NONE
INTEGER(4), INTENT(in):: nt, N
REAL(dp), INTENT(in) :: mu, sigma, gamma, r, T, S0, E  !E is exercise price
REAL(dp), INTENT(in), DIMENSION(N,2) ::V
REAL(dp), INTENT(out) :: f
INTEGER :: i, j, i1, m,d
REAL(dp)  ::dt, dS, S, dx, Vj
real(dp), INTENT(OUT), DIMENSION(nt+1,2) ::B 	!B stores investment boundary
REAL(dp), DIMENSION(:),  ALLOCATABLE:: F1, F2
REAL(dp)  ::qu, EF 
  dx=V(2,1)-V(1,1)
  dt = T/nt
  dS = sigma*sqrt(dt)
  qu = 0.5*(1+mu*sqrt(dt)/sigma)
  
  !Final period
    ALLOCATE (F1(1:(nt+2)))
    ALLOCATE (F2(1:(nt+2)))
    F1(:) = 0
    F2(:) = 0
  !---------------------
  
  !== Other periods ====
  DO i1= 1, (nt+1)
    i = nt+2-i1
    B(i,1)=(i-1)*dt
    IF (MOD(i, 2) == 0) THEN   !i is even
        DEALLOCATE (F2)
        ALLOCATE (F2(1:i))
	    d=0
        DO j = 1, i
            S = S0 + (j-1)*dS-((i-1)-(j-1))*dS
            EF = exp(-r*dt)*(qu*F1(j+1)+(1-qu)*F1(j))   
            IF (S>V(N,1)) THEN
                Vj = V(N,2)
            ELSE IF(S<V(1,1)) THEN
                Vj=V(1,2)
            ELSE
    	        m = floor((S-V(1,1))/dx)+1  
	            Vj = V(m,2)+(V(m+1,2)-V(m,2))*(S-V(m,1))/dx
            END IF
	        IF (Vj-E < EF) THEN
                F2(j) =  EF
	        ELSE
                F2(j) =  Vj-E
                d=d+1
	        END IF		
	        IF(d==1) THEN
		        B(i,2)=S
	        END IF
        END DO    

    ELSE        !i is odd
        DEALLOCATE (F1)
        ALLOCATE (F1(1:i))
	    d=0
        DO j = 1, i
            S = S0 + (j-1)*dS-((i-1)-(j-1))*dS
            EF = exp(-r*dt)*(qu*F2(j+1)+(1-qu)*F2(j))   
            IF (S>V(N,1)) THEN
                Vj = V(N,2)
            ELSE IF(S<V(1,1)) THEN
                Vj=V(1,2)
            ELSE
    	        m = floor((S-V(1,1))/dx)+1  
	            Vj = V(m,2)+(V(m+1,2)-V(m,2))*(S-V(m,1))/dx
            END IF
	        IF (Vj-E < EF) THEN
                F1(j) =  EF
	        ELSE
                F1(j) =  Vj-E
		        d=d+1                
	        END IF		
	        IF(d==1) THEN
		        B(i,2)=S
            END IF

        END DO    
    END IF    
  END DO
f = F1(1)
DEALLOCATE (F1)
DEALLOCATE (F2)
END SUBROUTINE binop


