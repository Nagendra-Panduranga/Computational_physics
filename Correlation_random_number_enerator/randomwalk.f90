!This program simulates the random walk Nw times with each random walk consisting of n steps. 
!The convention is +1 if btest(rand,b)is true and -1 is btest(rand,b) is -1.

!The 'rand0' variable acts as the seed for the random number generator.

! xstep is the value of each intermediate step decided upon the value of each bit of the random number.

!P(X) represents the actual probablity distribution.
!Cn(b,X) represents the number of times the final position X is reached for the bit 'b'
!X(b) represents the final posistion for the bit 'b' after n steps
!delta(b) represents the total deviation from the actual distribution
!deviation_i = \Delta(x) for a given 'b' 

PROGRAM randomwalk

  IMPLICIT NONE

  INTEGER(KIND(1_8)),parameter::a = 2862933555777941757_8, c = 1013904243_8 
  INTEGER(KIND(1_8)):: n, Nw, rand0, rand, b, i, j, xstep
  INTEGER(KIND(1_8)), ALLOCATABLE:: Cn(:,:), X(:) 
  REAL(KIND(1.0D0))::deviation_i, delta_sqr
  REAL(KIND(1.0D0)),ALLOCATABLE::P(:), delta(:)
  CHARACTER(8)::filename ! This character variable is to generate the name "pnn.dat"
  
  OPEN(UNIT = 20, FILE = "read.in")
  READ(20,*)n, Nw, rand0
  CLOSE(20)
  PRINT *,"n is",n,"Nw is " ,Nw, "seed is " ,rand0 ! debug line
  
  ALLOCATE (Cn(0:63,-n:n))
  ALLOCATE (X(0:63))
  ALLOCATE (P(-n:n))  
  ALLOCATE (delta(0:63))
  
  Cn = 0_8 ! Initializes the array to zero
  
  rand = rand0 !Initializes the random number

!********************************************************* This part of the program deos the random walks and simulates the counts for each bit**************
  
5 DO j = 1, Nw
     
  X = 0_8 ! Initializes the array to zero
  
10 DO i = 1, n
     
     rand = a * rand + c
     
 20     DO b = 0, 63
        
        SELECT CASE(btest(rand,b))
        CASE (.TRUE.)
           xstep = 1_8
        CASE (.FALSE.)
           xstep = -1_8
        END SELECT
        
        X(b) = X(b) + xstep
        
     END DO

  END DO

  DO b = 0, 63
     
     Cn(b,X(b)) = Cn(b, X(b)) + 1_8   

  END DO
  
END DO

!**************************************************This part calculates the actual distribution ***************************************

CALL dist_act(n,P) 

!**********************************************This part of the progran calculates the deviation***************************************

DO b = 0, 63
   
   delta_sqr = 0.0D0
   
   DO i = -n, n,2
      deviation_i = (dble(Cn(b,i))/dble(Nw)) - P(i)
      delta_sqr = delta_sqr + (deviation_i*deviation_i)
   END DO
   
   delta(b)=sqrt(delta_sqr)
 
END DO

!***********************************************This part writes the deviation to the file "d.dat"*************************************

OPEN(UNIT = 30, FILE = "d.dat")
12 FORMAT(I4,2X,F12.5)

DO b = 0, 63
   WRITE(30,12)b,(delta(b)*sqrt(dble(Nw)))
END DO

CLOSE(30)

!*************************************This part writes the distributon to the file "pnn.dat"***********************************************

DO b = 0, 63
   
   filename = "p"//achar(48 + ((b - (mod(b,10_8)))/10_8))//achar(48 + mod(b,10_8))//".dat"

   OPEN(UNIT = 40, FILE = filename)
   DO i = -n, n, 2
15    FORMAT(I4, 2X, F12.9, 2X, F12.9)
      IF (Cn(b,i).NE.0) THEN
      WRITE(40,15)i,dble(Cn(b,i))/dble(Nw),P(i)
   END IF
   END DO
   CLOSE(40)

END DO

END PROGRAM randomwalk

!-----------------------------------------------Actual Distribution -------------------------------------------------------------------------

SUBROUTINE dist_act(n,P)
  
  IMPLICIT NONE

  INTEGER(KIND(1_8)),INTENT(IN)::n
  INTEGER(KIND(1_8))::i, Factorial
  REAL(KIND(1.0D0)),PARAMETER::Pi = 3.14159265358979323846D0
  REAL(KIND(1_8)), DIMENSION(-n:n),INTENT(OUT)::P
  REAL(KIND(1.0D0))::n_r, i_r, P_ln1, P_ln2, P_ln3
  
  n_r = dble(n)
  P = 0.0D0
  DO i = -n, n, 2
     IF(ABS(i)==n) THEN 
        P(i) = Exp(-n_r * log(2.0D0))
     ELSE    
     i_r = dble(i)
     P_ln1 = (n_r + 0.50D0) * log(n_r) - (n_r * log(2.0D0)) - (0.50D0 * log(2.0D0*Pi)) ! The expressions were too long to defined in the same line. So there are 3 parts of same expression, P_ln1,...
     P_ln2 = - (0.50D0 * (n_r + i_r + 1.0D0)* log(0.50D0 * (n_r + i_r)))- (0.50D0 * (n_r - i_r +1.0D0 )* log(0.50D0 * (n_r - i_r)))
     P_ln3 =  (((1.0D0/n_r) - (1.0D0/(n_r-i_r)) - (1.0D0/(n_r + i_r))  )/12.0D0)     
     P(i) = Exp(P_ln1 + P_ln2 + P_ln3) 
     END IF
  END DO
  
END SUBROUTINE dist_act


