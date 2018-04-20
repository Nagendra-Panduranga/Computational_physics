
!------------------!
 module systemparam
 implicit none

 real(KIND(1.D0)), parameter :: pi = 3.141592653589793d0 
 real(KIND(1.D0)), parameter :: gme = 3.987d14  ! Earth's mass times constant of gravity
 real(KIND(1.D0)), parameter :: Ts = 86164.0D0 ! Siderial day
 real(KIND(1.D0)), SAVE:: alpha
 real(KIND(1.D0)), SAVE :: dt,dt2  ! time step, half of the time step

 end module systemparam
!----------------------!

PROGRAM geosat

 use systemparam

 IMPLICIT NONE

 REAL(KIND(1.D0)):: a, Tsat, rsat
 REAL(KIND(1.D0)):: x, y, z, vx, vy, vz, t
 REAL(KIND(1.D0)):: r, phi, theta, phi_prev, deltaR, deltaPhi
 INTEGER(KIND(1_8)):: Nt, tmax, Nw
 INTEGER(KIND(1_8)):: i, nr

Nt = 86164! PRINT *, "The number of divisions of siderial day: Nt"; READ *, Nt
tmax = 50 !PRINT *, "Number of days of monitoring the satellite"; READ *, tmax
Nw = tmax ! PRINT *,"Steps at which the data has to be written"; READ *, Nw
 
 alpha = 25.0D0 * (pi/180.D0)

 a = 1.d0 ! This defines the paramter 'a'

 Tsat = a * Ts

 rsat = (gme * Tsat * Tsat/ (4.d0 * pi * pi ))**(1.d0/3.d0)

 dt = Ts/dble(Nt) ; dt2 = dt / 2.d0
 
! ******************Initialization *******************************
 
 x = rsat
 y = 0.d0
 z = 0.d0
 
 vx =0.d0
 vy = sqrt(gme/rsat)   
 vz = 0.d0 
 nr = 0
 phi_prev= 0.0d0

 !****************************************************************
 OPEN(UNIT = 10, FILE = "sat.dat")
 
 DO i = 0 , Nt * tmax
    
    t = dble(i) * dt
    
    CALL polarposition(x, y, z, r, phi, theta)
    
    IF ((phi_prev - phi)/(2.0d0 * pi) .GE. 0.99d0) THEN
       nr = nr + 1
    END IF
    phi_prev = phi
    
    phi = (2.d0 * pi * dble(nr)) + phi
       
    IF(mod(i, Nw) == 0) THEN
       deltaPHI = phi - (2.d0 * pi * t /Ts)
       deltaR = (r - rsat)/ 1.d3
       WRITE(10,15)t, phi, r, deltaPhi, deltaR, theta
       
15     FORMAT(f10.2,2X,f10.5,2X,f17.5,2X, f12.5, 2X,f 12.5, 2X, f7.5)
    END IF
    
    CALL rkstep(t, x, y, z, vx, vy, vz)
    
 END DO
 
 CLOSE(10)
 
END PROGRAM geosat

!-----------------------------------!
 subroutine rkstep(t0,x0,y0,z0,vx0,vy0,vz0)
!---------------------------------------!
!Integrates the equations of motion one !
!time step, using the Runge-Kutta method!
!---------------------------------------!
 use systemparam
 implicit none

 REAL(KIND(1.D0)):: x0,y0,z0,x1,y1,z1,vx0,vy0,vz0,vx1,vy1,vz1,t0,th,t1,ax,ay,az
 REAL(KIND(1.D0)) :: kx1,kx2,kx3,kx4,ky1,ky2,ky3,ky4,kz1,kz2,kz3,kz4
 REAL(KIND(1.D0)) :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4,lz1,lz2,lz3,lz4

 t1=t0+dt; th=t0+dt2   
 call accel(x0,y0,z0,vx0,vy0,vz0,t0,ax,ay,az) 
 kx1=dt2*ax 
 ky1=dt2*ay
 kz1=dt2*az
 lx1=dt2*vx0 
 ly1=dt2*vy0
 lz1=dt2*vz0
 call accel(x0+lx1,y0+ly1,z0+lz1,vx0+kx1,vy0+ky1,vz0+kz1,th,ax,ay,az)
 kx2=dt2*ax
 ky2=dt2*ay
 kz2=dt2*az
 lx2=dt2*(vx0+kx1)
 ly2=dt2*(vy0+ky1)
 lz2=dt2*(vz0+kz1)
 call accel(x0+lx2,y0+ly2,z0+lz2,vx0+kx2,vy0+ky2,vz0+kz2,th,ax,ay,az)
 kx3=dt*ax
 ky3=dt*ay
 kz3=dt*az
 lx3=dt*(vx0+kx2)
 ly3=dt*(vy0+ky2)
 lz3=dt*(vz0+kz2)
 call accel(x0+lx3,y0+ly3,z0+lz3,vx0+kx3,vy0+ky3,vz0+kz3,t1,ax,ay,az)
 kx4=dt2*ax
 ky4=dt2*ay
 kz4=dt2*az
 lx4=dt2*(vx0+kx3)
 ly4=dt2*(vy0+ky3)
 lz4=dt2*(vz0+kz3)
 x1=x0+(lx1+2.d0*lx2+lx3+lx4)/3.d0
 y1=y0+(ly1+2.d0*ly2+ly3+ly4)/3.d0
 z1=z0+(lz1+2.d0*lz2+lz3+lz4)/3.d0
 vx1=vx0+(kx1+2.d0*kx2+kx3+kx4)/3.d0
 vy1=vy0+(ky1+2.d0*ky2+ky3+ky4)/3.d0
 vz1=vz0+(kz1+2.d0*kz2+kz3+kz4)/3.d0

 x0=x1; y0=y1; z0=z1
 vx0=vx1; vy0=vy1;vz0=vz1

 end subroutine rkstep
!---------------------!


!-----------------------------------!
 subroutine accel(x,y,z,vx,vy,vz,t,ax,ay,az)
!------------------------------------------------------!
!calculates the x, y and z components of the acceleration,!
!due to gravitation of earth and moon!
!------------------------------------------------------!
 use systemparam
 implicit none

 REAL(KIND(1.D0)) :: x, y, z, vx, vy, vz, t, ax, ay, az, r, r3
 REAL(KIND(1.D0)) :: xm, ym, zm ! the co-ordinates of the moon
 REAL(KIND(1.D0)) :: xsm, ysm, zsm, rsm !  satellite with respect to moon
 REAL(KIND(1.D0)) :: gmm, Tm, Rm ! mass, time period and radius of orbit of moon

 r=sqrt(x**2 + y**2 + z**2)       

 !*** evaluates the acceleration due to gravitation
 r3 = 1.d0/r**3
 ax = -gme * x *r3           
 ay = -gme*y*r3  
 az = -gme*z*r3  
 
 gmm = 1.2300d-2 * gme
 Tm = 27.25d0 * Ts
 Rm = 384400.d3

! Evaluates the position of moon

 xm = Rm * Cos(alpha) * Cos(2.d0 * pi * t / Tm)
 ym = Rm * Sin(2.d0 * pi * t / Tm)
 zm = Rm * Sin(alpha) * Cos(2.d0 * pi * t / Tm)
 
!*** evaluates the acceleration due to gravitation of moon
 
 xsm = x - xm
 ysm = y - ym
 zsm = z - zm
 rsm = sqrt(xsm**2 + ysm**2 + zsm**2)       
 r3 = 1.d0/rsm**3

 ax = ax - gmm*xsm*r3           
 ay = ay - gmm*ysm*r3           
 az = az - gmm*zsm*r3           

10  end subroutine accel
!--------------------!

!---------------------------------!
 subroutine polarposition(x,y,z,r,phi,theta)
!------------------------------------------------------!
!converts the x, y & z coordinates to a distance r, phi, theta  !
!------------------------------------------------------!

 USE systemparam

 implicit none
 
 real(KIND(1.D0)) :: x,y,z,r,rho,theta, phi

 r = sqrt(x**2 + y**2 + z**2)

 theta = asin(z/r)

 rho = sqrt(x**2 + y**2)! r * Cos(theta) 

if (y >= 0.d0) then
   phi = acos(x/rho)
 else
   phi= 2.d0*pi - acos(x/rho)
end if



 end subroutine polarposition
!----------------------------!
