!------------------------------------------------------------------------
!! subroutine computeStressStrainVolumePlaneStrain computes the stress
!! field associated with deforming dipping strain volumes considering the
!! following geometry using the analytical solution of
!!
!!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
!!   Associated with Distributed Anelastic Deformation in a Half Space,
!!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
!!
!! and considering the following geometry:
!!
!!              surface
!!      -------------+-------------- E (x2)
!!                   |
!!                   | dip /
!!                   |----/  . w
!!                   |   /     . i 
!!                   |  /        . d           
!!                   | /           . t     
!!                   |/              . h   
!!           q2,q3 ->@                 .
!!                  /|                   . 
!!                 / :                  /  
!!                /  |                 /  s
!!               /   :                /  s
!!              /    |               /  e
!!             /     :              /  n
!!               .   |             /  k
!!                 . :            /  c
!!                   .           /  i
!!                   : .        /  h
!!                   |   .     /  t
!!                   :     .  /  
!!                   |       .    
!!                   q3 (x3)
!!
!! INPUT:
!! r2, x3             east coordinates and depth of the observation point,
!! q2, q3             east and depth coordinates of the shear zone,
!! T, W               thickness and width of the shear zone,
!! dip (radian)       dip angle from horizontal of the shear zone,
!! epsijp             source strain component 22, 23 and 33 in the shear zone
!!                    in the system of reference tied to the shear zone,
!! G, lambda          shear modulus and Lame parameter.
!!
!! OUTPUT:
!! s22, s23, s33      stress components
!!
!! \author Sylvain Barbot (06/09/17) - original Fortran form
!-----------------------------------------------------------------------
SUBROUTINE computeStressStrainVolumePlaneStrain(r2,x3,q2,q3,T,W,dip, &
                              eps22p,eps23p,eps33p,G,lambda,s22,s23,s33)

  IMPLICIT NONE

  REAL*8, INTENT(IN) :: r2,x3,q2,q3,T,W,dip
  REAL*8, INTENT(IN) :: eps22p,eps23p,eps33p
  REAL*8, INTENT(IN) :: G,lambda
  REAL*8, INTENT(OUT) :: s22,s23,s33
    
  REAL*8 :: nu
  REAL*8 :: x2
  REAL*8 :: eps22,eps23,eps33,epskk
  REAL*8 :: e22,e23,e33,u23,u32,ekk

  REAL*8, EXTERNAL :: S,omega

  REAL*8, PARAMETER :: PI = 3.141592653589793115997963468544185161_8
  REAL*8, PARAMETER :: one = 1.00000000000000000000000000000000000_8
  REAL*8, PARAMETER :: zero = 0.00000000000000000000000000000000000_8

  ! Poisson s ratio
  nu=lambda/2._8/(lambda+G)

  ! check valid parameters
  IF ((-1._8 .GT. nu) .OR. (0.5_8 .LT. nu)) THEN
     WRITE (0,'("error: -1<=nu<=0.5, nu=",ES9.2E2," given.")') nu
     STOP 1
  END IF

  IF (0 .GT. x3) THEN
     WRITE (0,'("error: observation depth (x3) must be positive")')
     STOP 1
  END IF

  IF (0 .GT. q3) THEN
     WRITE (0,'("error: source depth (q3) must be positive")')
     STOP 1
  END IF

  ! convert strain to reference coordinate system R=[[SIN(dip), COS(dip)];[-COS(dip), SIN(dip)]]
  eps22=+SIN(dip)*(+eps22p*SIN(dip)+eps23p*COS(dip))+COS(dip)*(+eps23p*SIN(dip)+eps33p*COS(dip))
  eps23=+SIN(dip)*(-eps22p*COS(dip)+eps23p*SIN(dip))+COS(dip)*(-eps23p*COS(dip)+eps33p*SIN(dip))
  eps33=-COS(dip)*(-eps22p*COS(dip)+eps23p*SIN(dip))+SIN(dip)*(-eps23p*COS(dip)+eps33p*SIN(dip))

  ! isotropic strain
  epskk=eps22+eps33

  ! translation
  x2=r2-q2

  ! displacement gradients
  e22=IU22(T/2,W)-IU22(-T/2,W)+IU22(-T/2,zero)-IU22(T/2,zero)
  u23=IU23(T/2,W)-IU23(-T/2,W)+IU23(-T/2,zero)-IU23(T/2,zero)
  u32=IU32(T/2,W)-IU32(-T/2,W)+IU32(-T/2,zero)-IU32(T/2,zero)
  e33=IU33(T/2,W)-IU33(-T/2,W)+IU33(-T/2,zero)-IU33(T/2,zero)

  e22=e22        -eps22*Omega((x2*SIN(dip)-(x3-q3)*COS(dip))/T)*S((x2*COS(dip)+(x3-q3)*SIN(dip))/W)
  e23=(u23+u32)/2-eps23*Omega((x2*SIN(dip)-(x3-q3)*COS(dip))/T)*S((x2*COS(dip)+(x3-q3)*SIN(dip))/W)
  e33=e33        -eps33*Omega((x2*SIN(dip)-(x3-q3)*COS(dip))/T)*S((x2*COS(dip)+(x3-q3)*SIN(dip))/W)

  ! stress components
  s22=lambda*(e22+e33)+2*G*e22
  s23=2*G*e23
  s33=lambda*(e22+e33)+2*G*e33

CONTAINS

  !---------------------------------------------------------------
  !> function IU22
  !! computes the derivative of the indefinite integral U2
  !! with respect to x2p
  !! function IU22 is the integrand for displacement gradient u2,2
  !---------------------------------------------------------------
  REAL*8 FUNCTION IU22(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    IU22=1._8/(8._8*pi*G*(1._8-nu))*( &
     SIN(dip)*((lambda*epskk+2*G*eps22)*I223d2(y2p,y3p) &
                            +2*G*eps23*(I222d2(y2p,y3p)+I323d2(y2p,y3p)) &
              +(lambda*epskk+2*G*eps33)*I322d2(y2p,y3p) ) &
    +COS(dip)*((lambda*epskk+2*G*eps22)*I222d2(y2p,y3p) &
                            +2*G*eps23*(I322d2(y2p,y3p)-I223d2(y2p,y3p)) &
              -(lambda*epskk+2*G*eps33)*I323d2(y2p,y3p) ) &
        )
  END FUNCTION IU22
    
  !---------------------------------------------------------------
  !> function IU23
  !! computes the derivative of the indefinite integral U2
  !! with respect to x3p
  !! function IU23 is the integrand for displacement gradient u2,3
  !---------------------------------------------------------------
  REAL*8 FUNCTION IU23(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    IU23=1._8/(8._8*pi*G*(1._8-nu))*( &
     SIN(dip)*((lambda*epskk+2*G*eps22)*I223d3(y2p,y3p) &
                            +2*G*eps23*(I222d3(y2p,y3p)+I323d3(y2p,y3p)) &
              +(lambda*epskk+2*G*eps33)*I322d3(y2p,y3p) ) &
    +COS(dip)*((lambda*epskk+2*G*eps22)*I222d3(y2p,y3p) &
                            +2*G*eps23*(I322d3(y2p,y3p)-I223d3(y2p,y3p)) &
              -(lambda*epskk+2*G*eps33)*I323d3(y2p,y3p) ) &
        )

  END FUNCTION IU23
    
  !---------------------------------------------------------------
  !> function IU32
  !! computes the derivative of the indefinite integral U3
  !! with respect to x2p
  !! function IU32 is the integrand for displacement gradient u3,2
  !---------------------------------------------------------------
  REAL*8 FUNCTION IU32(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    IU32=1._8/(8._8*pi*G*(1._8-nu))*( &
     SIN(dip)*((lambda*epskk+2*G*eps22)*I233d2(y2p,y3p) &
                            +2*G*eps23*(I232d2(y2p,y3p)+I333d2(y2p,y3p)) &
              +(lambda*epskk+2*G*eps33)*I332d2(y2p,y3p) ) &
    +COS(dip)*((lambda*epskk+2*G*eps22)*I232d2(y2p,y3p) &
                            +2*G*eps23*(I332d2(y2p,y3p)-I233d2(y2p,y3p)) &
              -(lambda*epskk+2*G*eps33)*I333d2(y2p,y3p) ) &
        )
    
  END FUNCTION IU32
    
  !---------------------------------------------------------------
  !> function IU33
  !! computes the derivative of the indefinite integral U3
  !! with respect to x3p
  !! function IU33 is the integrand for displacement gradient u3,3
  !---------------------------------------------------------------
  REAL*8 FUNCTION IU33(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    IU33=1._8/(8._8*pi*G*(1._8-nu))*( &
     SIN(dip)*((lambda*epskk+2*G*eps22)*I233d3(y2p,y3p) &
                            +2*G*eps23*(I232d3(y2p,y3p)+I333d3(y2p,y3p)) &
              +(lambda*epskk+2*G*eps33)*I332d3(y2p,y3p) ) &
    +COS(dip)*((lambda*epskk+2*G*eps22)*I232d3(y2p,y3p) &
                            +2*G*eps23*(I332d3(y2p,y3p)-I233d3(y2p,y3p)) &
              -(lambda*epskk+2*G*eps33)*I333d3(y2p,y3p) ) &
        )
    
  END FUNCTION IU33
    
  !---------------------------------------------------------------
  !> function I223d2
  !! computes the derivative of the I integral I223,2
  !---------------------------------------------------------------
  REAL*8 FUNCTION I223d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I223d2=atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1._8)*x2*SIN(dip))**(-1)*( &
        y3p+(-1._8)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*(3._8+(-4._8)*nu+ &
        COS(2._8*dip))*SIN(dip)+(1._8/2._8)*atan2((-1._8)*(y2p+(-1._8)*(q3+x3)*COS( &
        dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)* &
        SIN(dip)),one)*((-1)*(7._8+4._8*nu*((-5._8)+4._8*nu))*SIN(dip)+(-1)*(3._8+( &
        -4._8)*nu)*SIN(3._8*dip))+(1._8/2._8)*((-2._8)*(3._8+(-4._8)*nu+COS(2._8*dip))*(( &
        -1)*y2p+(q3+(-1)*x3)*COS(dip)+x2*SIN(dip))*((-1)*q3+x3+y2p* &
        COS(dip)+(-1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+(x2+(-1)*y3p* &
        COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8)*y3p+8._8* &
        nu*y3p+(5._8+(-8._8)*nu)*x2*COS(dip)+x2*COS(3*dip)+(-2._8)*(q3+(-1) &
        *x3)*(4._8+(-4._8)*nu+COS(2*dip))*SIN(dip)+2*y2p*SIN(2*dip))+( &
        q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-10._8)*y2p+24._8*nu* &
        y2p+(-16._8)*nu**2*y2p+((13._8+(-28._8)*nu+16._8*nu**2)*q3+(11._8+(-28._8)*nu+ &
        16._8*nu**2)*x3)*COS(dip)+2*((-3._8)+4._8*nu)*y2p*COS(2._8*dip)+3._8* &
        q3*COS(3._8*dip)+(-4._8)*nu*q3*COS(3._8*dip)+5._8*x3*COS(3._8*dip)+(-4._8) &
        *nu*x3*COS(3*dip)+7._8*x2*SIN(dip)+(-20._8)*nu*x2*SIN(dip)+16._8* &
        nu**2*x2*SIN(dip)+3._8*x2*SIN(3._8*dip)+(-4._8)*nu*x2*SIN(3._8*dip))+ &
        (x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3* &
        x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+ &
        2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*(5._8+4._8* &
        nu*((-3._8)+2._8*nu))*y3p+(7._8+4._8*nu*((-5._8)+4._8*nu))*x2*COS(dip)+(3._8+( &
        -4._8)*nu)*x2*COS(3._8*dip)+(-1)*(13._8*q3+11._8*x3+(-28._8)*nu*(q3+x3)+ &
        16._8*nu**2*(q3+x3))*SIN(dip)+2*(3._8+(-4._8)*nu)*y2p*SIN(2*dip)+(( &
        -3._8)*q3+(-5._8)*x3+4._8*nu*(q3+x3))*SIN(3._8*dip))+(-2._8)*x3*(q3**2+ &
        x2**2+2._8*q3*x3+x3**2+y2p**2+y3p**2+(-2._8)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2._8*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        2._8*q3*COS(dip)+(-3._8)*y2p*COS(2*dip)+y2p*COS(4._8*dip)+(-2._8)*x2* &
        SIN(dip)+y3p*SIN(2*dip)+2._8*x2*SIN(3*dip)+(-1)*y3p*SIN(4._8* &
        dip))+(-4._8)*x3*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2._8)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2._8*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -2)*((-2._8)*q3*x2*COS(dip)+(3._8*x2*y2p+(q3+(-1)*x3)*y3p)*COS( &
        2._8*dip)+(-2)*y2p*y3p*COS(3._8*dip)+((-1)*x2*y2p+(q3+x3)*y3p)* &
        COS(4._8*dip)+((-1)*q3**2+x2**2+x3**2)*SIN(dip)+(3._8*q3*y2p+x3* &
        y2p+(-1)*x2*y3p)*SIN(2._8*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2) &
        *SIN(3._8*dip)+((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip)))+xLogy((1._8/4._8)*( &
        (5._8+(-8._8)*nu)*COS(dip)+COS(3._8*dip)),q3**2+x2**2+(-2)*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p*COS(dip)+( &
        -2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN( &
        dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((7._8+4._8*nu*((-5._8)+4._8*nu) &
        )*COS(dip)+(3._8+(-4._8)*nu)*COS(3._8*dip)),q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p*COS( &
        dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p* &
        SIN(dip)+2*x3*y3p*SIN(dip))
  
  END FUNCTION I223d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I223d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I223d3=(-1)*atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**(-1) &
        *(y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*COS(dip)*( &
        3._8+(-4._8)*nu+COS(2*dip))+(1._8/2._8)*atan2((-1)*(y2p+(-1)*(q3+x3)* &
        COS(dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+ &
        x3)*SIN(dip)),one)*((-1)*(11._8+(-28._8)*nu+16._8*nu**2)*COS(dip)+(-1) &
        *(5._8+(-4._8)*nu)*COS(3._8*dip))+(1._8/4._8)*((-4._8)*(3._8+(-4._8)*nu+COS(2._8*dip) &
        )*((-1)*y2p+(q3+(-1)*x3)*COS(dip)+x2*SIN(dip))*((-1)*x2+ &
        y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+2._8*((-1)*q3+ &
        x3+y2p*COS(dip)+(-1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS( &
        dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8) &
        *y3p+8._8*nu*y3p+(5._8+(-8._8)*nu)*x2*COS(dip)+x2*COS(3._8*dip)+(-2)* &
        (q3+(-1)*x3)*(4._8+(-4._8)*nu+COS(2*dip))*SIN(dip)+2._8*y2p*SIN(2._8* &
        dip))+(-2._8)*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2._8*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        (-10._8)*y2p+24._8*nu*y2p+(-16._8)*nu**2*y2p+((13._8+(-28._8)*nu+16._8*nu**2) &
        *q3+(11._8+(-28._8)*nu+16._8*nu**2)*x3)*COS(dip)+2._8*((-3._8)+4._8*nu)* &
        y2p*COS(2._8*dip)+3._8*q3*COS(3._8*dip)+(-4._8)*nu*q3*COS(3._8*dip)+5._8* &
        x3*COS(3._8*dip)+(-4._8)*nu*x3*COS(3._8*dip)+7._8*x2*SIN(dip)+(-20._8)* &
        nu*x2*SIN(dip)+16._8*nu**2*x2*SIN(dip)+3._8*x2*SIN(3._8*dip)+(-4._8)* &
        nu*x2*SIN(3._8*dip))+2._8*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip)) &
        *(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2._8)*(q3*y2p+x3* &
        y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
        **(-1)*((-2._8)*(5._8+4._8*nu*((-3._8)+2._8*nu))*y3p+(7._8+4._8*nu*((-5._8)+4._8*nu) &
        )*x2*COS(dip)+(3._8+(-4._8)*nu)*x2*COS(3._8*dip)+(-1)*(13._8*q3+11._8* &
        x3+(-28._8)*nu*(q3+x3)+16._8*nu**2*(q3+x3))*SIN(dip)+2._8*(3._8+(-4._8)* &
        nu)*y2p*SIN(2._8*dip)+((-3._8)*q3+(-5._8)*x3+4._8*nu*(q3+x3))*SIN(3._8* &
        dip))+(-8._8)*x3*SIN(dip)*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2._8*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(q3+(-2._8)*y2p*COS(dip)+2._8*( &
        q3+x3)*COS(2*dip)+(-1)*y2p*COS(3._8*dip)+y3p*SIN(3._8*dip))+(-8._8) &
        *x3*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2._8* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2._8)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2._8*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*((-2._8)*q3* &
        x2*COS(dip)+(3._8*x2*y2p+(q3+(-1)*x3)*y3p)*COS(2*dip)+(-2._8)* &
        y2p*y3p*COS(3._8*dip)+((-1)*x2*y2p+(q3+x3)*y3p)*COS(4._8*dip)+(( &
        -1)*q3**2+x2**2+x3**2)*SIN(dip)+(3._8*q3*y2p+x3*y2p+(-1)*x2* &
        y3p)*SIN(2._8*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2)*SIN(3._8*dip)+ &
        ((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip))+4._8*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2._8*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*q3*x2*COS( &
        dip)+(3._8*x2*y2p+(q3+(-1)*x3)*y3p)*COS(2._8*dip)+(-2)*y2p*y3p* &
        COS(3._8*dip)+((-1)*x2*y2p+(q3+x3)*y3p)*COS(4._8*dip)+((-1)* &
        q3**2+x2**2+x3**2)*SIN(dip)+(3._8*q3*y2p+x3*y2p+(-1)*x2*y3p)* &
        SIN(2._8*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2)*SIN(3._8*dip)+((q3+ &
        x3)*y2p+x2*y3p)*SIN(4._8*dip)))+xLogy((1._8/2._8)*(4._8+(-4._8)*nu+COS(2* &
        dip))*SIN(dip),q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2._8) &
        *q3*y2p*COS(dip)+2._8*x3*y2p*COS(dip)+(-2._8)*x2*y3p*COS(dip)+( &
        -2._8)*x2*y2p*SIN(dip)+2._8*q3*y3p*SIN(dip)+(-2._8)*x3*y3p*SIN( &
        dip))+xLogy((1._8/4._8)*((-1)*(11._8+(-28._8)*nu+16._8*nu**2)*SIN(dip)+((-5._8) &
        +4._8*nu)*SIN(3._8*dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+( &
        -2)*q3*y2p*COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS( &
        dip)+(-2._8)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2._8*x3*y3p*SIN( &
        dip))

  END FUNCTION I223d3

  !---------------------------------------------------------------
  REAL*8 FUNCTION I222d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I222d2=(-1)*atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*( &
        y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))**(-1),one)*COS(dip) &
        *((-3._8)+4._8*nu+COS(2*dip))+(1._8/2._8)*atan2((y2p+(-1)*(q3+x3)*COS( &
        dip)+(-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip))**(-1),one)*((-1)*((-7._8)+4._8*(5._8+(-4._8)*nu)*nu)*COS(dip)+(-1)* &
        (3._8+(-4._8)*nu)*COS(3._8*dip))+(1._8/2._8)*(2*((-3._8)+4._8*nu+COS(2*dip))*(( &
        -1)*y3p+x2*COS(dip)+((-1)*q3+x3)*SIN(dip))*(q3+(-1)*x3+(-1) &
        *y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+(x2+(-1)*y3p* &
        COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8)*y2p+8._8* &
        nu*y2p+(-2)*(q3+(-1)*x3)*COS(dip)*((-4._8)+4._8*nu+COS(2*dip))+( &
        -2)*x2*((-2)+4._8*nu+COS(2*dip))*SIN(dip)+2*y3p*SIN(2*dip))+( &
        x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3* &
        x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+ &
        2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*(5._8+4._8* &
        nu*((-3._8)+2._8*nu))*y2p+(13._8*q3+11*x3+(-28._8)*nu*(q3+x3)+16._8* &
        nu**2*(q3+x3))*COS(dip)+((-3._8)*q3+(-5._8)*x3+4._8*nu*(q3+x3))*COS( &
        3._8*dip)+(7._8+4._8*nu*((-5._8)+4._8*nu))*x2*SIN(dip)+2*(3._8+(-4._8)*nu)* &
        y3p*SIN(2*dip)+((-3._8)+4._8*nu)*x2*SIN(3._8*dip))+(q3+x3+(-1)*y2p* &
        COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(10._8*y3p+(-24._8)*nu*y3p+16._8* &
        nu**2*y3p+(-1)*(7._8+(-20._8)*nu+16._8*nu**2)*x2*COS(dip)+2*((-3._8)+ &
        4._8*nu)*y3p*COS(2*dip)+3._8*x2*COS(3*dip)+(-4._8)*nu*x2*COS(3._8* &
        dip)+13._8*q3*SIN(dip)+(-28._8)*nu*q3*SIN(dip)+16._8*nu**2*q3*SIN( &
        dip)+11._8*x3*SIN(dip)+(-28._8)*nu*x3*SIN(dip)+16._8*nu**2*x3*SIN( &
        dip)+(-3._8)*q3*SIN(3._8*dip)+4._8*nu*q3*SIN(3*dip)+(-5._8)*x3*SIN( &
        3._8*dip)+4._8*nu*x3*SIN(3._8*dip))+2*x3*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*x2*COS(dip)+ &
        3._8*y3p*COS(2*dip)+(-2)*x2*COS(3._8*dip)+y3p*COS(4._8*dip)+(-2)* &
        q3*SIN(dip)+y2p*SIN(2*dip)+y2p*SIN(4._8*dip))+(-4._8)*x3*(x2+(-1) &
        *y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*((-1)*((-1)*q3**2+ &
        x2**2+x3**2)*COS(dip)+((-1)*q3*y2p+x3*y2p+3*x2*y3p)*COS(2* &
        dip)+(-1)*(x2**2+(q3+x3)**2+2*y3p**2)*COS(3*dip)+((q3+x3)* &
        y2p+x2*y3p)*COS(4._8*dip)+(-2)*q3*x2*SIN(dip)+(x2*y2p+(3*q3+ &
        x3)*y3p)*SIN(2*dip)+(-2)*y2p*y3p*SIN(3*dip)+(x2*y2p+(-1)* &
        (q3+x3)*y3p)*SIN(4._8*dip)))+xLogy((-1._8/2._8)*((-2._8)+4._8*nu+COS(2*dip) &
        )*SIN(dip),q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)* &
        q3*y2p*COS(dip)+2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+( &
        -2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN( &
        dip))+xLogy((1._8/4._8)*((7._8+4._8*nu*((-5._8)+4._8*nu))*SIN(dip)+((-3._8)+4._8*nu) &
        *SIN(3._8*dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)* &
        q3*y2p*COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+ &
        (-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))
        
  END FUNCTION I222d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I222d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I222d3=(-1)*atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*( &
        y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))**(-1),one)*((-3._8)+ &
        4._8*nu+COS(2*dip))*SIN(dip)+(1._8/2._8)*atan2((y2p+(-1)*(q3+x3)*COS( &
        dip)+(-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip))**(-1),one)*((-11)*SIN(dip)+(-4._8)*nu*((-7._8)+4._8*nu)*SIN(dip)+ &
        (-1)*((-5._8)+4._8*nu)*SIN(3._8*dip))+(1._8/2._8)*(2*((-3._8)+4._8*nu+COS(2* &
        dip))*(y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))*((-1)* &
        x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+((-1)*q3+x3+ &
        y2p*COS(dip)+(-1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS( &
        dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8) &
        *y2p+8*nu*y2p+(-2)*(q3+(-1)*x3)*COS(dip)*((-4._8)+4._8*nu+COS( &
        2*dip))+(-2)*x2*((-2._8)+4._8*nu+COS(2*dip))*SIN(dip)+2*y3p*SIN( &
        2*dip))+(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+ &
        2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)* &
        COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)* &
        (5._8+4._8*nu*((-3._8)+2._8*nu))*y2p+(13._8*q3+11._8*x3+(-28._8)*nu*(q3+x3)+ &
        16._8*nu**2*(q3+x3))*COS(dip)+((-3._8)*q3+(-5._8)*x3+4._8*nu*(q3+x3))* &
        COS(3._8*dip)+(7._8+4._8*nu*((-5._8)+4._8*nu))*x2*SIN(dip)+2*(3._8+(-4._8)*nu) &
        *y3p*SIN(2*dip)+((-3._8)+4._8*nu)*x2*SIN(3._8*dip))+(x2+(-1)*y3p* &
        COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)* &
        x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-10._8)*y3p+24._8*nu*y3p+( &
        -16._8)*nu**2*y3p+(7._8+(-20._8)*nu+16._8*nu**2)*x2*COS(dip)+(-2)*((-3._8) &
        +4._8*nu)*y3p*COS(2*dip)+(-3._8)*x2*COS(3*dip)+4._8*nu*x2*COS(3._8* &
        dip)+(-13._8)*q3*SIN(dip)+28._8*nu*q3*SIN(dip)+(-16._8)*nu**2*q3* &
        SIN(dip)+(-11._8)*x3*SIN(dip)+28._8*nu*x3*SIN(dip)+(-16._8)*nu**2* &
        x3*SIN(dip)+3*q3*SIN(3*dip)+(-4._8)*nu*q3*SIN(3._8*dip)+5._8*x3* &
        SIN(3*dip)+(-4._8)*nu*x3*SIN(3._8*dip))+4._8*x3*COS(dip)*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        q3+(-2)*(q3+x3)*COS(2*dip)+y2p*COS(3._8*dip)+2*y3p*SIN(dip)+( &
        -1)*y3p*SIN(3*dip))+(-4._8)*x3*(q3+x3+(-1)*y2p*COS(dip)+y3p* &
        SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3* &
        y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)* &
        SIN(dip))**(-2)*((-1)*((-1)*q3**2+x2**2+x3**2)*COS(dip)+((-1) &
        *q3*y2p+x3*y2p+3*x2*y3p)*COS(2*dip)+(-1)*(x2**2+(q3+x3) &
        **2+2*y3p**2)*COS(3*dip)+((q3+x3)*y2p+x2*y3p)*COS(4._8*dip)+( &
        -2)*q3*x2*SIN(dip)+(x2*y2p+(3._8*q3+x3)*y3p)*SIN(2*dip)+(-2) &
        *y2p*y3p*SIN(3._8*dip)+(x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip)) &
        +2*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
        y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
        **(-1)*((-1)*((-1)*q3**2+x2**2+x3**2)*COS(dip)+((-1)*q3*y2p+ &
        x3*y2p+3*x2*y3p)*COS(2*dip)+(-1)*(x2**2+(q3+x3)**2+2* &
        y3p**2)*COS(3*dip)+((q3+x3)*y2p+x2*y3p)*COS(4*dip)+(-2)* &
        q3*x2*SIN(dip)+(x2*y2p+(3*q3+x3)*y3p)*SIN(2*dip)+(-2)* &
        y2p*y3p*SIN(3*dip)+(x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip)))+ &
        xLogy((1._8/2._8)*COS(dip)*((-4._8)+4._8*nu+COS(2*dip)),q3**2+x2**2+(-2)* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p* &
        COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3* &
        y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((11._8+(-28._8)* &
        nu+16._8*nu**2)*COS(dip)+((-5._8)+4._8*nu)*COS(3._8*dip)),q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p* &
        COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3* &
        y3p*SIN(dip)+2*x3*y3p*SIN(dip))

  END FUNCTION I222d3

  !---------------------------------------------------------------
  REAL*8 FUNCTION I233d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I233d2=4._8*((-1)+nu)*((-1)+2*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+( &
        -1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**( &
        -1),one)*COS(dip)+COS(2*dip)*(y2p+((-1)*q3+x3)*COS(dip)+(-1)* &
        x2*SIN(dip))*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+ &
        x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3* &
        y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN( &
        dip))**(-1)+4._8*((-1)+nu)*((-1)+2*nu)*(x2*COS(dip)+(-1)*(q3+ &
        x3)*SIN(dip))*((-1)*q3+(-1)*x3+y2p*COS(dip)+(-1)*y3p*SIN( &
        dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+ &
        x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN( &
        dip))**(-1)+4._8*((-1)+nu)*((-1)+2*nu)*y3p*(q3+x3+(-1)*y2p* &
        COS(dip)+y3p*SIN(dip))**(-1)*(1+((-1)*x2+y3p*COS(dip)+y2p* &
        SIN(dip))**2*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))**(-2))**( &
        -1)+atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**(-1)* &
        (y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*SIN(dip)* &
        SIN(2*dip)+((-3._8)+4._8*nu)*atan2((-1)*(y2p+(-1)*(q3+x3)*COS(dip) &
        +(-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip)),one)*SIN(dip)*SIN(2*dip)+(y2p+((-1)*q3+x3)*COS(dip)+(-1) &
        *x2*SIN(dip))*((-1)*q3+x3+y2p*COS(dip)+(-1)*y3p*SIN(dip))* &
        (q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1) &
        *x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p) &
        *SIN(dip))**(-1)*SIN(2*dip)+(-1)*SIN(dip)*(q3+x3+(-1)*y2p* &
        COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*((-3._8)+4._8*nu)*y2p*COS( &
        dip)+((-3._8)*q3+4._8*nu*q3+(-1)*x3+4._8*nu*x3)*COS(2*dip)+((-3._8)+ &
        4._8*nu)*(q3+(-1)*x3+x2*SIN(2*dip)))+(1._8/2._8)*(x2+(-1)*y3p*COS( &
        dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-8._8)*((-1._8)+nu)*((-1._8)+2*nu) &
        *y2p+(11*q3+x3+16*nu**2*(q3+x3)+(-4._8)*nu*(7._8*q3+3._8*x3))*COS( &
        dip)+2*((-3._8)+4._8*nu)*y2p*COS(2*dip)+(3._8*q3+x3+(-4._8)*nu*(q3+x3) &
        )*COS(3._8*dip)+(5._8+4._8*nu*((-5._8)+4._8*nu))*x2*SIN(dip)+(3._8+(-4._8)*nu) &
        *x2*SIN(3._8*dip))+2*x3*SIN(dip)*(q3**2+x2**2+2*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)* &
        x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(q3+(-2)*y2p*COS(dip)+ &
        y2p*COS(3._8*dip)+2*x2*SIN(2*dip)+(-1)*y3p*SIN(3._8*dip))+(-2)* &
        x3*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*(((-1)* &
        q3**2+x2**2+x3**2)*COS(dip)+(3._8*q3*y2p+x3*y2p+(-1)*x2*y3p)* &
        COS(2*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2)*COS(3._8*dip)+((q3+ &
        x3)*y2p+x2*y3p)*COS(4._8*dip)+2*q3*x2*SIN(dip)+((-3._8)*x2*y2p+ &
        ((-1)*q3+x3)*y3p)*SIN(2*dip)+2*y2p*y3p*SIN(3._8*dip)+(x2* &
        y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip))+xLogy((1._8/2._8)*COS(2*dip)* &
        SIN(dip),q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3* &
        y2p*COS(dip)+2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)* &
        x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+ &
        xLogy((1._8/4._8)*((5._8+4._8*nu*((-5._8)+4._8*nu))*SIN(dip)+(3._8+(-4._8)*nu)*SIN( &
        3._8*dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p* &
        COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)* &
        x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I233d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I233d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I233d3=(-4._8)*((-1._8)+nu)*((-1._8)+2*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+ &
        (-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**( &
        -1),one)*SIN(dip)+atan2((-1)*(y2p+(-1)*(q3+x3)*COS(dip)+(-1)* &
        x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip)),one) &
        *(3._8+(-4._8)*nu+((-1._8)+4._8*nu)*COS(2*dip))*SIN(dip)+COS(2*dip)*(( &
        -1)*y2p+(q3+(-1)*x3)*COS(dip)+x2*SIN(dip))*((-1)*q3+x3+y2p* &
        COS(dip)+(-1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+4._8*((-1._8)+nu)*( &
        (-1._8)+2._8*nu)*y3p*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+( &
        -4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*(x2*COS(dip)+(-1)*(q3+x3)*SIN( &
        dip))*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+(-1)*atan2( &
        (y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1) &
        *x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*COS(dip)*SIN(2*dip)+ &
        (y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*((-1)*x2+y3p* &
        COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*SIN(2*dip)+(-1)* &
        SIN(dip)*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*(( &
        -3._8)+4._8*nu)*y2p*COS(dip)+((-3)*q3+4._8*nu*q3+(-1)*x3+4._8*nu*x3) &
        *COS(2*dip)+((-3._8)+4._8*nu)*(q3+(-1)*x3+x2*SIN(2*dip)))+(1._8/2._8)* &
        (q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-8._8)*((-1._8)+nu)*(( &
        -1._8)+2._8*nu)*y2p+(11*q3+x3+16*nu**2*(q3+x3)+(-4._8)*nu*(7._8*q3+3._8* &
        x3))*COS(dip)+2*((-3._8)+4._8*nu)*y2p*COS(2*dip)+(3._8*q3+x3+(-4._8)* &
        nu*(q3+x3))*COS(3*dip)+(5._8+4._8*nu*((-5._8)+4._8*nu))*x2*SIN(dip)+( &
        3._8+(-4._8)*nu)*x2*SIN(3*dip))+2*x3*COS(dip)*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(q3+2*x3+( &
        -2)*(q3+x3)*COS(2*dip)+y2p*COS(3*dip)+2*y3p*SIN(dip)+(-1)* &
        y3p*SIN(3._8*dip))+(-2)*x3*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN( &
        dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+ &
        x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN( &
        dip))**(-2)*(((-1)*q3**2+x2**2+x3**2)*COS(dip)+(3._8*q3*y2p+x3* &
        y2p+(-1)*x2*y3p)*COS(2*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2) &
        *COS(3._8*dip)+((q3+x3)*y2p+x2*y3p)*COS(4._8*dip)+2*q3*x2*SIN( &
        dip)+((-3._8)*x2*y2p+((-1)*q3+x3)*y3p)*SIN(2._8*dip)+2*y2p*y3p* &
        SIN(3*dip)+(x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip))+(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        ((-1)*q3**2+x2**2+x3**2)*COS(dip)+(3*q3*y2p+x3*y2p+(-1)*x2* &
        y3p)*COS(2*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2)*COS(3._8*dip)+ &
        ((q3+x3)*y2p+x2*y3p)*COS(4._8*dip)+2*q3*x2*SIN(dip)+((-3._8)* &
        x2*y2p+((-1)*q3+x3)*y3p)*SIN(2*dip)+2*y2p*y3p*SIN(3._8*dip)+ &
        (x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip))+xLogy((-1._8/2._8)*COS(dip) &
        *COS(2*dip),q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)* &
        q3*y2p*COS(dip)+2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+( &
        -2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN( &
        dip))+xLogy((1._8/4._8)*((1._8+(-12._8)*nu+16._8*nu**2)*COS(dip)+(1._8+(-4._8)*nu) &
        *COS(3._8*dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)* &
        q3*y2p*COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+ &
        (-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))
        
  END FUNCTION I233d3

  !---------------------------------------------------------------
  REAL*8 FUNCTION I232d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I232d2=(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+ &
        (-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip)),one)*SIN(dip)+COS(2*dip)*(y3p+(-1)*x2*COS(dip)+(q3+(-1)* &
        x3)*SIN(dip))*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+ &
        x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3* &
        y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN( &
        dip))**(-1)+(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*((q3+x3)*COS(dip)+x2* &
        SIN(dip))*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+ &
        4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*y2p*(q3+x3+(-1)*y2p*COS(dip)+y3p* &
        SIN(dip))**(-1)*(1._8+((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))**2*( &
        q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))**(-2))**(-1)+(-1)*atan2( &
        (y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*(y3p+(-1)*x2* &
        COS(dip)+(q3+(-1)*x3)*SIN(dip))**(-1),one)*COS(dip)*SIN(2*dip)+ &
        ((-3._8)+4._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+(-1)*x2*SIN( &
        dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**(-1),one)*COS( &
        dip)*SIN(2*dip)+((-1)*y3p+x2*COS(dip)+((-1)*q3+x3)*SIN(dip)) &
        *(q3+(-1)*x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+( &
        -2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2* &
        y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**( &
        -1)*SIN(2*dip)+COS(dip)*((-1)*q3+(-1)*x3+y2p*COS(dip)+(-1)* &
        y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*( &
        q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p) &
        *SIN(dip))**(-1)*((((-3._8)+4._8*nu)*q3+((-1._8)+4._8*nu)*x3)*COS(2* &
        dip)+(-1)*((-3._8)+4._8*nu)*(q3+(-1)*x3+2*y3p*SIN(dip)+(-1)*x2* &
        SIN(2._8*dip)))+(1._8/2._8)*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip)) &
        *(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
        y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
        **(-1)*(8._8*((-1._8)+nu)*((-1._8)+2._8*nu)*y3p+((-5._8)+4._8*(5._8+(-4._8)*nu)* &
        nu)*x2*COS(dip)+2*((-3._8)+4._8*nu)*y3p*COS(2*dip)+(3._8+(-4._8)*nu)* &
        x2*COS(3._8*dip)+(11._8*q3+x3+16._8*nu**2*(q3+x3)+(-4._8)*nu*(7._8*q3+3._8* &
        x3))*SIN(dip)+((-3._8)*q3+(-1)*x3+4._8*nu*(q3+x3))*SIN(3._8*dip))+( &
        -2)*x3*COS(dip)*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2) &
        *(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)* &
        y3p)*SIN(dip))**(-1)*(q3+(-1)*y2p*COS(3._8*dip)+2*y3p*SIN(dip) &
        +(-2)*x2*SIN(2*dip)+y3p*SIN(3._8*dip))+(-2)*x3*(x2+(-1)*y3p* &
        COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)* &
        x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*((-2)*q3*x2*COS(dip)+( &
        x2*y2p+(3*q3+x3)*y3p)*COS(2*dip)+(-2)*y2p*y3p*COS(3._8*dip)+ &
        (x2*y2p+(-1)*(q3+x3)*y3p)*COS(4._8*dip)+((-1)*q3**2+x2**2+ &
        x3**2)*SIN(dip)+(q3*y2p+(-1)*x3*y2p+(-3)*x2*y3p)*SIN(2* &
        dip)+(x2**2+(q3+x3)**2+2*y3p**2)*SIN(3._8*dip)+(-1)*((q3+x3)* &
        y2p+x2*y3p)*SIN(4._8*dip))+xLogy((1._8/2._8)*COS(dip)*COS(2._8*dip), &
        q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS( &
        dip)+2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p* &
        SIN(dip)+2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8) &
        *(((-5._8)+4._8*(5._8+(-4._8)*nu)*nu)*COS(dip)+(3._8+(-4._8)*nu)*COS(3._8*dip)) &
        ,q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip) &
        +(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p* &
        SIN(dip)+2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I232d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I232d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    REAL*8, EXTERNAL :: xLogy

    I232d3=(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+ &
        (-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip)),one)*COS(dip)+(-1)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+(-1)* &
        x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**(-1),one) &
        *COS(dip)*(3._8+(-4._8)*nu+(1._8+(-4._8)*nu)*COS(2*dip))+COS(2*dip)*( &
        y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))*(q3+(-1)*x3+(-1) &
        *y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+4._8*((-1._8)+nu)*( &
        (-1._8)+2._8*nu)*y2p*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+( &
        -4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*((q3+x3)*COS(dip)+x2*SIN(dip))*(( &
        -1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+(-1)*atan2((y2p+(( &
        -1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip) &
        +(q3+(-1)*x3)*SIN(dip))**(-1),one)*SIN(dip)*SIN(2*dip)+(y3p+( &
        -1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))*((-1)*x2+y3p*COS( &
        dip)+y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*SIN(2*dip)+(-1)* &
        COS(dip)*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((((-3._8)+4._8* &
        nu)*q3+((-1._8)+4._8*nu)*x3)*COS(2._8*dip)+(-1)*((-3._8)+4._8*nu)*(q3+( &
        -1)*x3+2*y3p*SIN(dip)+(-1)*x2*SIN(2*dip)))+(1._8/2._8)*(q3+x3+( &
        -1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)* &
        x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(8._8*((-1._8)+nu)*((-1._8)+2._8* &
        nu)*y3p+((-5._8)+4._8*(5._8+(-4._8)*nu)*nu)*x2*COS(dip)+2*((-3._8)+4._8*nu) &
        *y3p*COS(2*dip)+(3._8+(-4._8)*nu)*x2*COS(3._8*dip)+(11._8*q3+x3+16._8* &
        nu**2*(q3+x3)+(-4._8)*nu*(7._8*q3+3*x3))*SIN(dip)+((-3._8)*q3+(-1)* &
        x3+4._8*nu*(q3+x3))*SIN(3._8*dip))+2*x3*SIN(dip)*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(q3+2*x3+( &
        -2)*y2p*COS(dip)+2*(q3+x3)*COS(2._8*dip)+(-1)*y2p*COS(3._8*dip)+ &
        y3p*SIN(3._8*dip))+(-2)*x3*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN( &
        dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+ &
        x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN( &
        dip))**(-2)*((-2)*q3*x2*COS(dip)+(x2*y2p+(3*q3+x3)*y3p)* &
        COS(2*dip)+(-2)*y2p*y3p*COS(3*dip)+(x2*y2p+(-1)*(q3+x3)* &
        y3p)*COS(4._8*dip)+((-1)*q3**2+x2**2+x3**2)*SIN(dip)+(q3*y2p+( &
        -1)*x3*y2p+(-3)*x2*y3p)*SIN(2*dip)+(x2**2+(q3+x3)**2+2* &
        y3p**2)*SIN(3*dip)+(-1)*((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip))+( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -1)*((-2)*q3*x2*COS(dip)+(x2*y2p+(3*q3+x3)*y3p)*COS(2* &
        dip)+(-2)*y2p*y3p*COS(3*dip)+(x2*y2p+(-1)*(q3+x3)*y3p)* &
        COS(4._8*dip)+((-1)*q3**2+x2**2+x3**2)*SIN(dip)+(q3*y2p+(-1)* &
        x3*y2p+(-3)*x2*y3p)*SIN(2*dip)+(x2**2+(q3+x3)**2+2*y3p**2)* &
        SIN(3._8*dip)+(-1)*((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip))+xLogy((1._8/2._8) &
        *COS(2._8*dip)*SIN(dip),q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p*COS(dip)+(-2)*x2* &
        y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+(-2)* &
        x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((1._8+(-12._8)*nu+16._8*nu**2)*SIN(dip) &
        +((-1._8)+4._8*nu)*SIN(3._8*dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2* &
        y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2*x3* &
        y3p*SIN(dip))
      
  END FUNCTION I232d3

  !---------------------------------------------------------------
  REAL*8 FUNCTION I323d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I323d2=(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+ &
        (-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**( &
        -1),one)*COS(dip)+COS(2*dip)*(y2p+((-1)*q3+x3)*COS(dip)+(-1)* &
        x2*SIN(dip))*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+ &
        x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3* &
        y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN( &
        dip))**(-1)+(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*(x2*COS(dip)+(-1)*( &
        q3+x3)*SIN(dip))*((-1)*q3+(-1)*x3+y2p*COS(dip)+(-1)*y3p* &
        SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3* &
        y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)* &
        SIN(dip))**(-1)+(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*y3p*(q3+x3+(-1)* &
        y2p*COS(dip)+y3p*SIN(dip))**(-1)*(1._8+((-1)*x2+y3p*COS(dip)+ &
        y2p*SIN(dip))**2*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))**(-2) &
        )**(-1)+atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**( &
        -1)*(y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*SIN(dip) &
        *SIN(2._8*dip)+((-3._8)+4._8*nu)*atan2((-1)*(y2p+(-1)*(q3+x3)*COS( &
        dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)* &
        SIN(dip)),one)*SIN(dip)*SIN(2*dip)+(y2p+((-1)*q3+x3)*COS(dip)+( &
        -1)*x2*SIN(dip))*((-1)*q3+x3+y2p*COS(dip)+(-1)*y3p*SIN(dip) &
        )*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+( &
        -1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3* &
        y3p)*SIN(dip))**(-1)*SIN(2._8*dip)+(-1)*SIN(dip)*(q3+x3+(-1)* &
        y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*((-3._8)+4._8*nu)*y2p*COS( &
        dip)+((-3._8)*q3+4._8*nu*q3+(-5._8)*x3+4._8*nu*x3)*COS(2*dip)+((-3._8)+ &
        4._8*nu)*(q3+(-1)*x3+x2*SIN(2._8*dip)))+(1._8/2._8)*(x2+(-1)*y3p*COS( &
        dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(8._8*((-1._8)+nu)*((-1._8)+2._8*nu)* &
        y2p+(-1)*((5._8+4._8*nu*((-5._8)+4._8*nu))*q3+(19._8+4._8*nu*((-9._8)+4._8*nu))* &
        x3)*COS(dip)+2*((-3._8)+4._8*nu)*y2p*COS(2._8*dip)+(3._8*q3+5._8*x3+(-4._8) &
        *nu*(q3+x3))*COS(3._8*dip)+((-11._8)+4._8*(7._8+(-4._8)*nu)*nu)*x2*SIN( &
        dip)+(3._8+(-4._8)*nu)*x2*SIN(3._8*dip))+(-2)*x3*SIN(dip)*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        q3+(-2)*y2p*COS(dip)+y2p*COS(3._8*dip)+2*x2*SIN(2*dip)+(-1)* &
        y3p*SIN(3._8*dip))+(-2)*x3*(x2+(-1)*y3p*COS(dip)+(-1)*y2p* &
        SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3* &
        y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)* &
        SIN(dip))**(-2)*((q3**2+(-1)*x2**2+(-1)*x3**2)*COS(dip)+((-1) &
        *(3._8*q3+x3)*y2p+x2*y3p)*COS(2._8*dip)+(x2**2+(q3+x3)**2+2* &
        y2p**2)*COS(3*dip)+(-1)*((q3+x3)*y2p+x2*y3p)*COS(4._8*dip)+( &
        -2)*q3*x2*SIN(dip)+(3*x2*y2p+(q3+(-1)*x3)*y3p)*SIN(2._8*dip) &
        +(-2)*y2p*y3p*SIN(3._8*dip)+((-1)*x2*y2p+(q3+x3)*y3p)*SIN(4._8* &
        dip))+xLogy((1._8/2._8)*COS(2._8*dip)*SIN(dip),q3**2+x2**2+(-2)*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p*COS(dip)+( &
        -2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN( &
        dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*(((-11._8)+4._8*(7._8+(-4._8)*nu) &
        *nu)*SIN(dip)+(3._8+(-4._8)*nu)*SIN(3._8*dip)),q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p*COS( &
        dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p* &
        SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I323d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I323d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I323d3=4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+( &
        -1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**( &
        -1),one)*SIN(dip)+atan2((-1)*(y2p+(-1)*(q3+x3)*COS(dip)+(-1)* &
        x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip)),one) &
        *(3._8+(-4._8)*nu+((-5._8)+4._8*nu)*COS(2._8*dip))*SIN(dip)+COS(2._8*dip)*(( &
        -1)*y2p+(q3+(-1)*x3)*COS(dip)+x2*SIN(dip))*((-1)*q3+x3+y2p* &
        COS(dip)+(-1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+(-4._8)*((-1._8)+nu) &
        *((-1._8)+2._8*nu)*y3p*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x2*COS(dip)+(-1)*(q3+x3)*SIN( &
        dip))*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+(-1)*atan2( &
        (y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1) &
        *x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*COS(dip)*SIN(2*dip)+ &
        (y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*((-1)*x2+y3p* &
        COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*SIN(2*dip)+(-1)* &
        SIN(dip)*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*(( &
        -3._8)+4._8*nu)*y2p*COS(dip)+((-3._8)*q3+4._8*nu*q3+(-5._8)*x3+4._8*nu*x3) &
        *COS(2*dip)+((-3._8)+4._8*nu)*(q3+(-1)*x3+x2*SIN(2*dip)))+(1._8/2._8)* &
        (q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(8._8*((-1._8)+nu)*((-1._8) &
        +2._8*nu)*y2p+(-1)*((5._8+4._8*nu*((-5._8)+4._8*nu))*q3+(19._8+4._8*nu*((-9._8)+ &
        4._8*nu))*x3)*COS(dip)+2*((-3._8)+4._8*nu)*y2p*COS(2*dip)+(3._8*q3+ &
        5._8*x3+(-4._8)*nu*(q3+x3))*COS(3*dip)+((-11._8)+4._8*(7._8+(-4._8)*nu)*nu) &
        *x2*SIN(dip)+(3._8+(-4._8)*nu)*x2*SIN(3._8*dip))+(-2)*x3*COS(dip)* &
        (q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -1)*(q3+2*x3+(-2)*(q3+x3)*COS(2*dip)+y2p*COS(3._8*dip)+2* &
        y3p*SIN(dip)+(-1)*y3p*SIN(3*dip))+(-2)*x3*(q3+x3+(-1)*y2p* &
        COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*((q3**2+(-1)*x2**2+(-1)* &
        x3**2)*COS(dip)+((-1)*(3*q3+x3)*y2p+x2*y3p)*COS(2*dip)+( &
        x2**2+(q3+x3)**2+2*y2p**2)*COS(3*dip)+(-1)*((q3+x3)*y2p+x2* &
        y3p)*COS(4._8*dip)+(-2)*q3*x2*SIN(dip)+(3*x2*y2p+(q3+(-1)*x3) &
        *y3p)*SIN(2._8*dip)+(-2)*y2p*y3p*SIN(3*dip)+((-1)*x2*y2p+( &
        q3+x3)*y3p)*SIN(4._8*dip))+(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((q3**2+(-1)*x2**2+(-1)* &
        x3**2)*COS(dip)+((-1)*(3*q3+x3)*y2p+x2*y3p)*COS(2*dip)+( &
        x2**2+(q3+x3)**2+2*y2p**2)*COS(3*dip)+(-1)*((q3+x3)*y2p+x2* &
        y3p)*COS(4._8*dip)+(-2)*q3*x2*SIN(dip)+(3*x2*y2p+(q3+(-1)*x3) &
        *y3p)*SIN(2*dip)+(-2)*y2p*y3p*SIN(3*dip)+((-1)*x2*y2p+( &
        q3+x3)*y3p)*SIN(4._8*dip))+xLogy((-1._8/2._8)*COS(dip)*COS(2._8*dip), &
        q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS( &
        dip)+2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p* &
        SIN(dip)+2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8) &
        *((-1)*(19._8+4._8*nu*((-9._8)+4._8*nu))*COS(dip)+(5._8+(-4._8)*nu)*COS(3._8* &
        dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p* &
        COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)* &
        x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I323d3

  !---------------------------------------------------------------
  REAL*8 FUNCTION I322d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I322d2=4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+( &
        -1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip)),one)*SIN(dip)+2*((-3._8)+4._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS( &
        dip)+(-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip))**(-1),one)*COS(dip)**2*SIN(dip)+COS(2*dip)*(y3p+(-1)*x2* &
        COS(dip)+(q3+(-1)*x3)*SIN(dip))*((-1)*x2+y3p*COS(dip)+y2p* &
        SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*( &
        q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)*q3* &
        y3p+x3*y3p)*SIN(dip))**(-1)+4*((-1)+nu)*((-1._8)+2._8*nu)*((q3+x3) &
        *COS(dip)+x2*SIN(dip))*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip) &
        )*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
        y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
        **(-1)+(-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*y2p*(q3+x3+(-1)*y2p*COS( &
        dip)+y3p*SIN(dip))**(-1)*(1._8+((-1)*x2+y3p*COS(dip)+y2p*SIN( &
        dip))**2*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))**(-2))**(-1)+( &
        -1)*atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*(y3p+ &
        (-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))**(-1),one)*COS(dip)* &
        SIN(2*dip)+((-1)*y3p+x2*COS(dip)+((-1)*q3+x3)*SIN(dip))*(q3+ &
        (-1)*x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+(-2)* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)* &
        COS(dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)* &
        SIN(2._8*dip)+(-1)*COS(dip)*(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN( &
        dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+ &
        x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN( &
        dip))**(-1)*((((-3._8)+4._8*nu)*q3+((-5._8)+4._8*nu)*x3)*COS(2._8*dip)+( &
        -1)*((-3._8)+4._8*nu)*(q3+(-1)*x3+2*y3p*SIN(dip)+(-1)*x2*SIN(2._8* &
        dip)))+(1._8/2._8)*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -1)*((-8._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*y3p+(11._8+4._8*nu*((-7._8)+4._8*nu)) &
        *x2*COS(dip)+2*((-3._8)+4._8*nu)*y3p*COS(2*dip)+(3._8+(-4._8)*nu)* &
        x2*COS(3._8*dip)+(-1)*(5._8*q3+19._8*x3+16._8*nu**2*(q3+x3)+(-4._8)*nu*( &
        5._8*q3+9._8*x3))*SIN(dip)+((-3._8)*q3+(-5._8)*x3+4._8*nu*(q3+x3))*SIN( &
        3._8*dip))+2*x3*COS(dip)*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*(q3+(-1)*y2p*COS(3._8*dip)+2* &
        y3p*SIN(dip)+(-2)*x2*SIN(2*dip)+y3p*SIN(3*dip))+(-2)*x3*( &
        x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+2*q3* &
        x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+ &
        2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*(2*q3*x2*COS( &
        dip)+(-1)*(x2*y2p+(3*q3+x3)*y3p)*COS(2._8*dip)+2*y2p*y3p* &
        COS(3._8*dip)+((-1)*x2*y2p+(q3+x3)*y3p)*COS(4._8*dip)+(-1)*((-1) &
        *q3**2+x2**2+x3**2)*SIN(dip)+((-1)*q3*y2p+x3*y2p+3*x2*y3p) &
        *SIN(2._8*dip)+(-1)*(x2**2+(q3+x3)**2+2*y3p**2)*SIN(3._8*dip)+(( &
        q3+x3)*y2p+x2*y3p)*SIN(4._8*dip))+xLogy((1._8/2._8)*COS(dip)*COS(2._8* &
        dip),q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p* &
        COS(dip)+2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2* &
        y2p*SIN(dip)+2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy( &
        (1._8/4._8)*((11._8+4._8*nu*((-7._8)+4._8*nu))*COS(dip)+(3._8+(-4._8)*nu)*COS(3._8* &
        dip)),q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p* &
        COS(dip)+(-2)*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)* &
        x2*y2p*SIN(dip)+2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I322d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I322d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I322d3=4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+( &
        -1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN( &
        dip)),one)*COS(dip)+atan2((y2p+(-1)*(q3+x3)*COS(dip)+(-1)*x2* &
        SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**(-1),one)* &
        COS(dip)*((-3._8)+4._8*nu+((-5._8)+4._8*nu)*COS(2*dip))+COS(2*dip)*( &
        y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))*(q3+(-1)*x3+(-1) &
        *y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
        *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+(-4._8)*((-1._8)+nu) &
        *((-1._8)+2._8*nu)*y2p*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((q3+x3)*COS(dip)+x2*SIN(dip))* &
        ((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)+(-1)*atan2((y2p+(( &
        -1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip) &
        +(q3+(-1)*x3)*SIN(dip))**(-1),one)*SIN(dip)*SIN(2*dip)+(y3p+( &
        -1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))*((-1)*x2+y3p*COS( &
        dip)+y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*SIN(2*dip)+(-1)* &
        COS(dip)*((-1)*x2+y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS( &
        dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((((-3._8)+4._8* &
        nu)*q3+((-5._8)+4._8*nu)*x3)*COS(2._8*dip)+(-1)*((-3._8)+4._8*nu)*(q3+( &
        -1)*x3+2*y3p*SIN(dip)+(-1)*x2*SIN(2*dip)))+(1._8/2._8)*(q3+x3+( &
        -1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)* &
        x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-8._8)*((-1._8)+nu)*((-1._8)+ &
        2._8*nu)*y3p+(11._8+4._8*nu*((-7._8)+4._8*nu))*x2*COS(dip)+2*((-3._8)+4._8*nu) &
        *y3p*COS(2*dip)+(3._8+(-4._8)*nu)*x2*COS(3*dip)+(-1)*(5._8*q3+19._8* &
        x3+16._8*nu**2*(q3+x3)+(-4._8)*nu*(5._8*q3+9._8*x3))*SIN(dip)+((-3._8)* &
        q3+(-5._8)*x3+4._8*nu*(q3+x3))*SIN(3._8*dip))+(-2._8)*x3*SIN(dip)*( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -1)*(q3+2*x3+(-2)*y2p*COS(dip)+2*(q3+x3)*COS(2*dip)+(-1)* &
        y2p*COS(3._8*dip)+y3p*SIN(3._8*dip))+(-2)*x3*(q3+x3+(-1)*y2p* &
        COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*(2*q3*x2*COS(dip)+(-1)*( &
        x2*y2p+(3._8*q3+x3)*y3p)*COS(2*dip)+2*y2p*y3p*COS(3._8*dip)+(( &
        -1)*x2*y2p+(q3+x3)*y3p)*COS(4._8*dip)+(-1)*((-1)*q3**2+x2**2+ &
        x3**2)*SIN(dip)+((-1)*q3*y2p+x3*y2p+3*x2*y3p)*SIN(2*dip)+( &
        -1)*(x2**2+(q3+x3)**2+2*y3p**2)*SIN(3._8*dip)+((q3+x3)*y2p+x2* &
        y3p)*SIN(4._8*dip))+(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2) &
        *(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)* &
        y3p)*SIN(dip))**(-1)*(2*q3*x2*COS(dip)+(-1)*(x2*y2p+(3*q3+ &
        x3)*y3p)*COS(2._8*dip)+2*y2p*y3p*COS(3._8*dip)+((-1)*x2*y2p+( &
        q3+x3)*y3p)*COS(4._8*dip)+(-1)*((-1)*q3**2+x2**2+x3**2)*SIN( &
        dip)+((-1)*q3*y2p+x3*y2p+3*x2*y3p)*SIN(2*dip)+(-1)*(x2**2+ &
        (q3+x3)**2+2*y3p**2)*SIN(3*dip)+((q3+x3)*y2p+x2*y3p)*SIN(4._8* &
        dip))+xLogy((1._8/2._8)*COS(2*dip)*SIN(dip),q3**2+x2**2+(-2)*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p*COS(dip)+( &
        -2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p*SIN( &
        dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((-1)*(19._8+(-36._8)*nu+ &
        16._8*nu**2)*SIN(dip)+((-5._8)+4._8*nu)*SIN(3._8*dip)),q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p* &
        COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3* &
        y3p*SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I322d3

  !---------------------------------------------------------------
  REAL*8 FUNCTION I333d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I333d2=(-1)*atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**(-1) &
        *(y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*((-3._8)+4._8* &
        nu+COS(2._8*dip))*SIN(dip)+(1._8/2._8)*atan2((-1)*(y2p+(-1)*(q3+x3)* &
        COS(dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+ &
        x3)*SIN(dip)),one)*(((-13._8)+4._8*(7._8+(-4._8)*nu)*nu)*SIN(dip)+(3._8+(-4._8) &
        *nu)*SIN(3*dip))+(1._8/2._8)*(2*((-3._8)+4._8*nu+COS(2._8*dip))*((-1)* &
        y2p+(q3+(-1)*x3)*COS(dip)+x2*SIN(dip))*((-1)*q3+x3+y2p*COS( &
        dip)+(-1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+(x2+(-1)*y3p*COS( &
        dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8)*y3p+8*nu* &
        y3p+(7._8+(-8._8)*nu)*x2*COS(dip)+(-1)*x2*COS(3._8*dip)+2*(q3+(-1)* &
        x3)*((-2._8)+4._8*nu+COS(2._8*dip))*SIN(dip)+(-2._8)*y2p*SIN(2*dip))+( &
        q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-10._8)*y2p+24._8*nu* &
        y2p+(-16._8)*nu**2*y2p+((7._8+(-20._8)*nu+16._8*nu**2)*q3+(5._8+(-20._8)*nu+ &
        16._8*nu**2)*x3)*COS(dip)+(-2)*((-3._8)+4._8*nu)*y2p*COS(2*dip)+( &
        -3._8)*q3*COS(3._8*dip)+4._8*nu*q3*COS(3._8*dip)+(-1)*x3*COS(3._8*dip)+ &
        4._8*nu*x3*COS(3*dip)+13._8*x2*SIN(dip)+(-28._8)*nu*x2*SIN(dip)+ &
        16._8*nu**2*x2*SIN(dip)+(-3._8)*x2*SIN(3._8*dip)+4._8*nu*x2*SIN(3._8* &
        dip))+(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+ &
        2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)* &
        COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)* &
        (5._8+4._8*nu*((-3)+2*nu))*y3p+(13._8+4._8*nu*((-7._8)+4._8*nu))*x2*COS( &
        dip)+((-3._8)+4._8*nu)*x2*COS(3._8*dip)+(-1)*(7._8*q3+5._8*x3+(-20._8)*nu*( &
        q3+x3)+16._8*nu**2*(q3+x3))*SIN(dip)+2*((-3._8)+4._8*nu)*y2p*SIN(2* &
        dip)+(3._8*q3+x3+(-4._8)*nu*(q3+x3))*SIN(3*dip))+(-2)*x3*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        2*q3*COS(dip)+(-3._8)*y2p*COS(2*dip)+y2p*COS(4._8*dip)+(-2)*x2* &
        SIN(dip)+y3p*SIN(2*dip)+2*x2*SIN(3._8*dip)+(-1)*y3p*SIN(4._8* &
        dip))+(-4._8)*x3*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*( &
        q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+ &
        x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**( &
        -2)*((-2)*q3*x2*COS(dip)+(3*x2*y2p+(q3+(-1)*x3)*y3p)*COS( &
        2*dip)+(-2)*y2p*y3p*COS(3*dip)+((-1)*x2*y2p+(q3+x3)*y3p)* &
        COS(4._8*dip)+((-1)*q3**2+x2**2+x3**2)*SIN(dip)+(3._8*q3*y2p+x3* &
        y2p+(-1)*x2*y3p)*SIN(2._8*dip)+(-1)*(x2**2+(q3+x3)**2+2*y2p**2) &
        *SIN(3._8*dip)+((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip)))+xLogy((1._8/4._8)*( &
        (7._8+(-8._8)*nu)*COS(dip)+(-1)*COS(3._8*dip)),q3**2+x2**2+(-2)*q3* &
        x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p*COS( &
        dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p* &
        SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((13._8+4._8*nu*((-7._8)+ &
        4._8*nu))*COS(dip)+((-3._8)+4._8*nu)*COS(3._8*dip)),q3**2+x2**2+2*q3* &
        x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p*COS( &
        dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p* &
        SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I333d2

  !---------------------------------------------------------------
  REAL*8 FUNCTION I333d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I333d3=atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))**(-1)*( &
        y3p+(-1)*x2*COS(dip)+(q3+(-1)*x3)*SIN(dip)),one)*COS(dip)*(( &
        -3._8)+4._8*nu+COS(2*dip))+(1._8/2._8)*atan2((-1)*(y2p+(-1)*(q3+x3)*COS( &
        dip)+(-1)*x2*SIN(dip))**(-1)*(y3p+(-1)*x2*COS(dip)+(q3+x3)* &
        SIN(dip)),one)*((-1)*(5._8+(-20._8)*nu+16._8*nu**2)*COS(dip)+(1._8+(-4._8)* &
        nu)*COS(3._8*dip))+(1._8/2._8)*(2*((-3._8)+4._8*nu+COS(2*dip))*((-1)*y2p+ &
        (q3+(-1)*x3)*COS(dip)+x2*SIN(dip))*((-1)*x2+y3p*COS(dip)+ &
        y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2) &
        *(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1)* &
        q3*y3p+x3*y3p)*SIN(dip))**(-1)+((-1)*q3+x3+y2p*COS(dip)+(-1) &
        *y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+( &
        -2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1) &
        *q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8)*y3p+8._8*nu*y3p+(7._8+(-8._8) &
        *nu)*x2*COS(dip)+(-1)*x2*COS(3._8*dip)+2*(q3+(-1)*x3)*((-2._8)+ &
        4._8*nu+COS(2._8*dip))*SIN(dip)+(-2)*y2p*SIN(2*dip))+((-1)*x2+ &
        y3p*COS(dip)+y2p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-10._8)*y2p+24._8*nu*y2p+(-16._8) &
        *nu**2*y2p+((7._8+(-20._8)*nu+16._8*nu**2)*q3+(5._8+(-20._8)*nu+16._8*nu**2) &
        *x3)*COS(dip)+(-2)*((-3._8)+4._8*nu)*y2p*COS(2._8*dip)+(-3._8)*q3* &
        COS(3._8*dip)+4._8*nu*q3*COS(3._8*dip)+(-1)*x3*COS(3*dip)+4._8*nu* &
        x3*COS(3._8*dip)+13._8*x2*SIN(dip)+(-28._8)*nu*x2*SIN(dip)+16._8* &
        nu**2*x2*SIN(dip)+(-3._8)*x2*SIN(3._8*dip)+4._8*nu*x2*SIN(3._8*dip))+ &
        (q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+ &
        x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*( &
        (-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*(5._8+4._8*nu*(( &
        -3._8)+2._8*nu))*y3p+(13._8+4._8*nu*((-7._8)+4._8*nu))*x2*COS(dip)+((-3._8)+4._8* &
        nu)*x2*COS(3*dip)+(-1)*(7._8*q3+5*x3+(-20._8)*nu*(q3+x3)+16._8* &
        nu**2*(q3+x3))*SIN(dip)+2*((-3._8)+4._8*nu)*y2p*SIN(2*dip)+(3._8* &
        q3+x3+(-4._8)*nu*(q3+x3))*SIN(3._8*dip))+(-4._8)*x3*SIN(dip)*(q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
        y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
        q3+(-2)*y2p*COS(dip)+2*(q3+x3)*COS(2*dip)+(-1)*y2p*COS(3._8* &
        dip)+y3p*SIN(3*dip))+(-4._8)*x3*(q3+x3+(-1)*y2p*COS(dip)+y3p* &
        SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3* &
        y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)* &
        SIN(dip))**(-2)*((-2)*q3*x2*COS(dip)+(3._8*x2*y2p+(q3+(-1)*x3) &
        *y3p)*COS(2*dip)+(-2)*y2p*y3p*COS(3*dip)+((-1)*x2*y2p+( &
        q3+x3)*y3p)*COS(4._8*dip)+((-1)*q3**2+x2**2+x3**2)*SIN(dip)+(3._8* &
        q3*y2p+x3*y2p+(-1)*x2*y3p)*SIN(2*dip)+(-1)*(x2**2+(q3+x3) &
        **2+2*y2p**2)*SIN(3*dip)+((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip))+ &
        2*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
        y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
        **(-1)*((-2)*q3*x2*COS(dip)+(3._8*x2*y2p+(q3+(-1)*x3)*y3p)* &
        COS(2._8*dip)+(-2)*y2p*y3p*COS(3._8*dip)+((-1)*x2*y2p+(q3+x3)* &
        y3p)*COS(4._8*dip)+((-1)*q3**2+x2**2+x3**2)*SIN(dip)+(3._8*q3*y2p+ &
        x3*y2p+(-1)*x2*y3p)*SIN(2*dip)+(-1)*(x2**2+(q3+x3)**2+2* &
        y2p**2)*SIN(3._8*dip)+((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip)))+xLogy(( &
        -1._8/2._8)*((-2._8)+4._8*nu+COS(2*dip))*SIN(dip),q3**2+x2**2+(-2)*q3* &
        x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p*COS( &
        dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3*y3p* &
        SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((-1)*(5._8+(-20._8)* &
        nu+16._8*nu**2)*SIN(dip)+(1._8+(-4._8)*nu)*SIN(3._8*dip)),q3**2+x2**2+2* &
        q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p* &
        COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3* &
        y3p*SIN(dip)+2*x3*y3p*SIN(dip))
      
  END FUNCTION I333d3
      
  !---------------------------------------------------------------
  REAL*8 FUNCTION I332d2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    IF (0 .EQ. MOD(dip,360._8)) THEN
       I332d2=(1._8/2._8)*(8._8*x3*(x2+(-1)*y3p)*(x2**2+q3*x3+x3**2+(-1)*x3*y2p+( &
              -2)*x2*y3p+y3p**2)*(q3**2+x2**2+x3**2+2*q3*(x3+(-1)*y2p)+( &
              -2)*x3*y2p+y2p**2+(-2)*x2*y3p+y3p**2)**(-2)+(-16._8)*((-1._8)+nu) &
              **2*(q3+x3+(-1)*y2p)*(x2+(-1)*y3p)*(q3**2+x2**2+x3**2+2*q3* &
              (x3+(-1)*y2p)+(-2)*x3*y2p+y2p**2+(-2)*x2*y3p+y3p**2)**(-1)+ &
              8._8*((-1._8)+nu)*(q3+(-1)*x3+(-1)*y2p)*(x2+(-1)*y3p)*(q3**2+ &
              x2**2+x3**2+2*x3*y2p+y2p**2+(-2)*q3*(x3+y2p)+(-2)*x2*y3p+ &
              y3p**2)**(-1)+(-2)*((-3._8)+4._8*nu)*(q3+(-1)*x3+(-1)*y2p)*(x2+( &
              -1)*y3p)*(q3**2+x2**2+x3**2+2*x3*y2p+y2p**2+(-2)*q3*(x3+y2p) &
              +(-2)*x2*y3p+y3p**2)**(-1)+(10._8*q3+6._8*x3+(-24._8)*nu*(q3+x3)+16._8* &
              nu**2*(q3+x3)+(-2)*(5._8+4._8*nu*((-3._8)+2._8*nu))*y2p)*(x2+(-1)*y3p) &
              *(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
              y2p+x2*y3p))**(-1)+2*x3*((-4._8)*x2+4._8*y3p)*(q3**2+x2**2+2*q3* &
              x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p))**(-1))+8._8* &
              ((-1._8)+nu)**2*atan2((q3+x3+(-1)*y2p)*(x2+(-1)*y3p)**(-1),one)+( &
              -4._8)*((-1._8)+nu)*atan2(((-1)*q3+x3+y2p)*((-1)*x2+y3p)**(-1),one)
    ELSE
       IF (180 .EQ. MOD(dip,360._8)) THEN
          I332d2=(1._8/2._8)*(x2+y3p)*((-8._8)*((-1._8)+nu)*(q3+(-1)*x3+y2p)*(q3**2+ &
              x2**2+(-2)*q3*x3+x3**2+2*q3*y2p+(-2)*x3*y2p+y2p**2+2*x2* &
              y3p+y3p**2)**(-1)+2*((-3._8)+4._8*nu)*(q3+(-1)*x3+y2p)*(q3**2+ &
              x2**2+(-2)*q3*x3+x3**2+2*q3*y2p+(-2)*x3*y2p+y2p**2+2*x2* &
              y3p+y3p**2)**(-1)+(-8._8)*x3*(x2**2+q3*x3+x3**2+x3*y2p+2*x2* &
              y3p+y3p**2)*(q3**2+x2**2+x3**2+2*x3*y2p+y2p**2+2*q3*(x3+y2p)+ &
              2*x2*y3p+y3p**2)**(-2)+8._8*x3*(q3**2+x2**2+x3**2+2*x3*y2p+ &
              y2p**2+2*q3*(x3+y2p)+2*x2*y3p+y3p**2)**(-1)+16._8*((-1._8)+nu)**2* &
              (q3+x3+y2p)*(q3**2+x2**2+x3**2+2*x3*y2p+y2p**2+2*q3*(x3+y2p)+ &
              2*x2*y3p+y3p**2)**(-1)+(-2)*(5._8*q3+3*x3+(-12._8)*nu*(q3+x3)+8._8* &
              nu**2*(q3+x3)+(5._8+4._8*nu*((-3._8)+2._8*nu))*y2p)*(q3**2+x2**2+x3**2+ &
              2*x3*y2p+y2p**2+2*q3*(x3+y2p)+2*x2*y3p+y3p**2)**(-1))+4._8*(( &
              -1._8)+nu)*atan2((q3+(-1)*x3+y2p)*(x2+y3p)**(-1),one)+(-8._8)*((-1._8)+ &
              nu)**2*atan2((q3+x3+y2p)*(x2+y3p)**(-1),one)
       ELSE
          I332d2= &
            atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*(y3p+(-1) &
              *x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))**(-1),one)*COS(dip)*(3._8+( &
              -4._8)*nu+COS(2._8*dip))+(1._8/2._8)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+( &
              -1)*x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**( &
              -1),one)*((-1)*((-13._8)+4._8*(7._8+(-4._8)*nu)*nu)*COS(dip)+(-1)*((-3._8)+ &
              4._8*nu)*COS(3._8*dip))+(1._8/4._8)*((-4._8)*(3._8+(-4._8)*nu+COS(2._8*dip))*((-1) &
              *y3p+x2*COS(dip)+((-1)*q3+x3)*SIN(dip))*(q3+(-1)*x3+(-1)* &
              y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+ &
              y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2) &
              *(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)+2*(x2+(-1)* &
              y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+ &
              x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS( &
              dip)+(-2)*(x2*y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8) &
              *y2p+8._8*nu*y2p+2*(q3+(-1)*x3)*COS(dip)*(2._8+(-4._8)*nu+COS(2._8* &
              dip))+2*x2*(4._8+(-4._8)*nu+COS(2._8*dip))*SIN(dip)+(-2._8)*y3p*SIN(2._8* &
              dip))+2*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip))*(q3**2+ &
              x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2* &
              y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*( &
              (-2)*(5._8+4._8*nu*((-3._8)+2._8*nu))*y2p+(7._8*q3+5._8*x3+(-20._8)*nu*(q3+x3) &
              +16._8*nu**2*(q3+x3))*COS(dip)+(3._8*q3+x3+(-4._8)*nu*(q3+x3))*COS( &
              3._8*dip)+(13._8+4._8*nu*((-7._8)+4._8*nu))*x2*SIN(dip)+2*((-3._8)+4._8*nu)* &
              y3p*SIN(2*dip)+(3._8+(-4._8)*nu)*x2*SIN(3._8*dip))+(-2)*(y3p+(-1)* &
              x2*COS(dip)+(q3+x3)*SIN(dip))**(-2)*(q3+x3+(-1)*y2p*COS(dip)+ &
              y3p*SIN(dip))*(1+((-1)*y2p+(q3+x3)*COS(dip)+x2*SIN(dip))**2* &
              (y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**(-2))**(-1)*((-2)*( &
              5._8+4._8*nu*((-3._8)+2._8*nu))*y3p+(13._8+(-28._8)*nu+16._8*nu**2)*x2*COS(dip) &
              +2*((-3._8)+4._8*nu)*y3p*COS(2*dip)+(-1)*((-3._8)+4._8*nu)*x2*COS(3._8* &
              dip)+(-4._8)*nu*((-5._8)+4._8*nu)*(q3+x3)*SIN(dip)+(-1)*(7._8*q3+5._8*x3) &
              *SIN(dip)+(-1)*(3._8*q3+x3+(-4._8)*nu*(q3+x3))*SIN(3._8*dip))+4._8* &
              x3*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
              y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
              **(-1)*((-2)*x2*COS(dip)+3._8*y3p*COS(2*dip)+(-2)*x2*COS(3._8* &
              dip)+y3p*COS(4._8*dip)+(-2)*q3*SIN(dip)+y2p*SIN(2*dip)+y2p* &
              SIN(4._8*dip))+(-8)*x3*(x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip) &
              )*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
              y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
              **(-2)*((-1)*((-1)*q3**2+x2**2+x3**2)*COS(dip)+((-1)*q3*y2p+ &
              x3*y2p+3*x2*y3p)*COS(2*dip)+(-1)*(x2**2+(q3+x3)**2+2* &
              y3p**2)*COS(3._8*dip)+((q3+x3)*y2p+x2*y3p)*COS(4._8*dip)+(-2)* &
              q3*x2*SIN(dip)+(x2*y2p+(3._8*q3+x3)*y3p)*SIN(2*dip)+(-2)* &
              y2p*y3p*SIN(3._8*dip)+(x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip)))+ &
              xLogy((1._8/2._8)*(4._8+(-4._8)*nu+COS(2*dip))*SIN(dip),q3**2+x2**2+(-2)* &
              q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2*x3*y2p* &
              COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3* &
              y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((13._8+4._8*nu*(( &
              -7._8)+4._8*nu))*SIN(dip)+(3._8+(-4._8)*nu)*SIN(3._8*dip)),q3**2+x2**2+2* &
              q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)*x3*y2p* &
              COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2*q3* &
              y3p*SIN(dip)+2*x3*y3p*SIN(dip))
       END IF
    END IF

  END FUNCTION I332d2
      
  !---------------------------------------------------------------
  REAL*8 FUNCTION I332d3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
      
    REAL*8, EXTERNAL :: xLogy

    I332d3=atan2((y2p+((-1)*q3+x3)*COS(dip)+(-1)*x2*SIN(dip))*(y3p+(-1) &
        *x2*COS(dip)+(q3+(-1)*x3)*SIN(dip))**(-1),one)*(3._8+(-4._8)*nu+COS( &
        2*dip))*SIN(dip)+(1._8/2._8)*atan2((y2p+(-1)*(q3+x3)*COS(dip)+(-1) &
        *x2*SIN(dip))*(y3p+(-1)*x2*COS(dip)+(q3+x3)*SIN(dip))**(-1), &
        1._8)*((-5._8)*SIN(dip)+(-4._8)*nu*((-5._8)+4._8*nu)*SIN(dip)+(-1)*(1._8+(-4._8) &
        *nu)*SIN(3._8*dip))+(1._8/2._8)*(2*(3._8+(-4._8)*nu+COS(2*dip))*((-1)* &
        y3p+x2*COS(dip)+((-1)*q3+x3)*SIN(dip))*((-1)*x2+y3p*COS(dip) &
        +y2p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+( &
        -2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2*y2p+(-1) &
        *q3*y3p+x3*y3p)*SIN(dip))**(-1)+((-1)*q3+x3+y2p*COS(dip)+( &
        -1)*y3p*SIN(dip))*(q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+(-1)*x3*y2p+x2*y3p)*COS(dip)+(-2)*(x2* &
        y2p+(-1)*q3*y3p+x3*y3p)*SIN(dip))**(-1)*((-6._8)*y2p+8._8*nu* &
        y2p+2*(q3+(-1)*x3)*COS(dip)*(2._8+(-4._8)*nu+COS(2._8*dip))+2*x2*( &
        4._8+(-4._8)*nu+COS(2*dip))*SIN(dip)+(-2)*y3p*SIN(2*dip))+(q3+x3+( &
        -1)*y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+ &
        y2p**2+y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)* &
        x2*y2p+(q3+x3)*y3p)*SIN(dip))**(-1)*((-2)*(5._8+4._8*nu*((-3._8)+2._8* &
        nu))*y2p+(7._8*q3+5._8*x3+(-20._8)*nu*(q3+x3)+16._8*nu**2*(q3+x3))* &
        COS(dip)+(3._8*q3+x3+(-4._8)*nu*(q3+x3))*COS(3._8*dip)+(13._8+4._8*nu*(( &
        -7._8)+4._8*nu))*x2*SIN(dip)+2*((-3._8)+4._8*nu)*y3p*SIN(2*dip)+(3._8+( &
        -4._8)*nu)*x2*SIN(3._8*dip))+((-1)*x2+y3p*COS(dip)+y2p*SIN(dip)) &
        *(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*(q3*y2p+x3* &
        y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p)*SIN(dip)) &
        **(-1)*(10._8*y3p+(-24._8)*nu*y3p+16._8*nu**2*y3p+(-1)*(13._8+(-28._8)* &
        nu+16._8*nu**2)*x2*COS(dip)+(-2)*((-3._8)+4._8*nu)*y3p*COS(2*dip)+( &
        -3._8)*x2*COS(3._8*dip)+4._8*nu*x2*COS(3._8*dip)+7._8*q3*SIN(dip)+(-20._8) &
        *nu*q3*SIN(dip)+16._8*nu**2*q3*SIN(dip)+5._8*x3*SIN(dip)+(-20._8)* &
        nu*x3*SIN(dip)+16._8*nu**2*x3*SIN(dip)+3._8*q3*SIN(3*dip)+(-4._8)* &
        nu*q3*SIN(3._8*dip)+x3*SIN(3._8*dip)+(-4._8)*nu*x3*SIN(3._8*dip))+4._8* &
        x3*COS(dip)*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*( &
        q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3)*y3p) &
        *SIN(dip))**(-1)*(q3+(-2)*(q3+x3)*COS(2*dip)+y2p*COS(3._8*dip) &
        +2*y3p*SIN(dip)+(-1)*y3p*SIN(3*dip))+(-4._8)*x3*(q3+x3+(-1)* &
        y2p*COS(dip)+y3p*SIN(dip))*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+ &
        y3p**2+(-2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2* &
        y2p+(q3+x3)*y3p)*SIN(dip))**(-2)*((-1)*((-1)*q3**2+x2**2+ &
        x3**2)*COS(dip)+((-1)*q3*y2p+x3*y2p+3*x2*y3p)*COS(2*dip)+( &
        -1)*(x2**2+(q3+x3)**2+2*y3p**2)*COS(3._8*dip)+((q3+x3)*y2p+x2* &
        y3p)*COS(4._8*dip)+(-2)*q3*x2*SIN(dip)+(x2*y2p+(3*q3+x3)*y3p) &
        *SIN(2*dip)+(-2)*y2p*y3p*SIN(3._8*dip)+(x2*y2p+(-1)*(q3+x3)* &
        y3p)*SIN(4._8*dip))+2*(q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+( &
        -2)*(q3*y2p+x3*y2p+x2*y3p)*COS(dip)+2*((-1)*x2*y2p+(q3+x3) &
        *y3p)*SIN(dip))**(-1)*((-1)*((-1)*q3**2+x2**2+x3**2)*COS( &
        dip)+((-1)*q3*y2p+x3*y2p+3*x2*y3p)*COS(2*dip)+(-1)*(x2**2+ &
        (q3+x3)**2+2*y3p**2)*COS(3*dip)+((q3+x3)*y2p+x2*y3p)*COS(4._8* &
        dip)+(-2)*q3*x2*SIN(dip)+(x2*y2p+(3._8*q3+x3)*y3p)*SIN(2*dip) &
        +(-2)*y2p*y3p*SIN(3._8*dip)+(x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8* &
        dip)))+xLogy((-1._8/2._8)*COS(dip)*(2._8+(-4._8)*nu+COS(2._8*dip)),q3**2+ &
        x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+2* &
        x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+ &
        2*q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip))+xLogy((1._8/4._8)*((5._8+( &
        -20._8)*nu+16._8*nu**2)*COS(dip)+(1._8+(-4._8)*nu)*COS(3._8*dip)),q3**2+ &
        x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2)* &
        x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+ &
        2*q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip))

  END FUNCTION I332d3
      
END








