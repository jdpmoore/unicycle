!------------------------------------------------------------------------
!! subroutine computeDisplacementStrainVolumePlaneStrain computes the
!! displacement field associated with deforming dipping strain volumes 
!! considering the following geometry using the analytical solution of
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
!! G, lambda          shear modulus and Lame parameter in the half space.
!!
!! OUTPUT:
!! u2                 displacement component in the east direction,
!! u3                 displacement component in the depth direction.
!!
!! \author Sylvain Barbot (06/15/17) - original Fortran form
!-----------------------------------------------------------------------
SUBROUTINE computeDisplacementStrainVolumePlaneStrain(r2,x3,q2,q3,T,W,dip, &
                              eps22p,eps23p,eps33p,G,lambda,u2,u3)

  IMPLICIT NONE

  REAL*8, INTENT(IN) :: r2,x3,q2,q3,T,W,dip
  REAL*8, INTENT(IN) :: eps22p,eps23p,eps33p
  REAL*8, INTENT(IN) :: G,lambda
  REAL*8, INTENT(OUT) :: u2,u3
    
  REAL*8 :: nu
  REAL*8 :: x2
  REAL*8 :: eps22,eps23,eps33,epskk

  REAL*8, EXTERNAL :: S,omega

  REAL*8, PARAMETER :: PI = 3.141592653589793115997963468544185161_8

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

  ! displacement components
  u2=IU2(T/2,W)-IU2(-T/2,W)+IU2(-T/2,0._8)-IU2(T/2,0._8)
  u3=IU3(T/2,W)-IU3(-T/2,W)+IU3(-T/2,0._8)-IU3(T/2,0._8)

CONTAINS

  !---------------------------------------------------------------
  !> function IU2
  !! is the integrand for displacement u2
  !---------------------------------------------------------------
  REAL*8 FUNCTION IU2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    IU2=1._8/(8._8*pi*G*(1._8-nu))*( &
       SIN(dip)*((lambda*epskk+2*G*eps22)*I223(y2p,y3p) &
                               +2*G*eps23*(I222(y2p,y3p)+I323(y2p,y3p)) &
                 +(lambda*epskk+2*G*eps33)*I322(y2p,y3p) ) &
      +COS(dip)*((lambda*epskk+2*G*eps22)*I222(y2p,y3p) &
                               +2*G*eps23*(I322(y2p,y3p)-I223(y2p,y3p)) &
                 -(lambda*epskk+2*G*eps33)*I323(y2p,y3p) ) &
        )

  END FUNCTION IU2
    
  !---------------------------------------------------------------
  !> function IU3
  !! is the integrand for displacement gradient u3
  !---------------------------------------------------------------
  REAL*8 FUNCTION IU3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    IU3=1._8/(8._8*pi*G*(1._8-nu))*( &
       SIN(dip)*((lambda*epskk+2*G*eps22)*I233(y2p,y3p) &
                               +2*G*eps23*(I232(y2p,y3p)+I333(y2p,y3p)) &
                 +(lambda*epskk+2*G*eps33)*I332(y2p,y3p) ) &
      +COS(dip)*((lambda*epskk+2*G*eps22)*I232(y2p,y3p) &
                               +2*G*eps23*(I332(y2p,y3p)-I233(y2p,y3p)) &
                 -(lambda*epskk+2*G*eps33)*I333(y2p,y3p) ) &
        )

  END FUNCTION IU3
    
  REAL*8 FUNCTION y2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    y2=y2p*SIN(dip)+y3p*COS(dip)

  END FUNCTION y2

  REAL*8 FUNCTION y3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    y3=-y2p*COS(dip)+y3p*SIN(dip)+q3

  END FUNCTION y3

  REAL*8 FUNCTION r1p(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    r1p=SQRT((x3-y3(y2p,y3p))**2+(x2-y2(y2p,y3p))**2)

  END FUNCTION r1p

  REAL*8 FUNCTION r2p(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    r2p=SQRT((x3+y3(y2p,y3p))**2+(x2-y2(y2p,y3p))**2)

  END FUNCTION r2p

  REAL*8 FUNCTION r12p(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    r12p= q3**2+x2**2+(-2)*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+ &
          2*x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2* &
          q3*y3p*SIN(dip)+(-2)*x3*y3p*SIN(dip)

  END FUNCTION r12p

  REAL*8 FUNCTION p2(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    p2 = (x3-q3)*COS(dip)-SIN(dip)*x2+y2p

  END FUNCTION p2

  REAL*8 FUNCTION p3(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    p3 =-(x3-q3)*SIN(dip)-COS(dip)*x2+y3p

  END FUNCTION p3

  REAL*8 FUNCTION p2p(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    p2p=-(q3+x3)*COS(dip)-SIN(dip)*x2+y2p

  END FUNCTION p2p

  REAL*8 FUNCTION p3p(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    p3p= (q3+x3)*SIN(dip)-COS(dip)*x2+y3p

  END FUNCTION p3p

  REAL*8 FUNCTION r22p(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p
    
    r22p= q3**2+x2**2+2*q3*x3+x3**2+y2p**2+y3p**2+(-2)*q3*y2p*COS(dip)+(-2) &
         *x3*y2p*COS(dip)+(-2)*x2*y3p*COS(dip)+(-2)*x2*y2p*SIN(dip)+2* &
          q3*y3p*SIN(dip)+2*x3*y3p*SIN(dip)

  END FUNCTION r22p

  !---------------------------------------------------------------
  !> function I223
  !! computes the I integral I223
  !---------------------------------------------------------------
  REAL*8 FUNCTION I223(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I223= &
    (-1)*p2(y2p,y3p)*atan2(p3(y2p,y3p)/p2(y2p,y3p),1._8)*(3._8-4._8*nu+COS(2._8*dip)) &
    +(1._8/4._8)*log(r12p(y2p,y3p))*( &
        -6._8*y3p+8._8*y3p*nu &
        +x2*(5._8-8._8*nu)*COS(dip) &
        +x2*COS(3._8*dip) &
        -2*(q3-x3)*(4._8-4._8*nu+COS(2._8*dip))*SIN(dip) &
        +2*y2p*SIN(2._8*dip) &
    ) &
    +(1._8/2._8)*atan2(-p3p(y2p,y3p)/p2p(y2p,y3p),1._8)*( &
        2*y2p*(5._8+4._8*nu*(-3._8+2._8*nu)) &
        -(13._8*q3+11._8*x3-28._8*(q3+x3)*nu+16._8*(q3+x3)*nu**2)*COS(dip) &
        -2*y2p*(-3._8+4._8*nu)*COS(2*dip) &
        -(3._8*q3+5._8*x3-4._8*(q3+x3)*nu)*COS(3._8*dip) &
        -x2*(7._8+4._8*nu*(-5._8+4._8*nu))*SIN(dip) &
        -x2*(3._8-4._8*nu)*SIN(3._8*dip) &
    ) &
    +(1._8/4._8)*log(r22p(y2p,y3p))*( &
        -2*y3p*(5._8+4._8*nu*(-3._8+2._8*nu)) &
        +x2*(7._8+4._8*nu*(-5._8+4._8*nu))*COS(dip) &
        +x2*(3._8-4._8*nu)*COS(3._8*dip) &
        -(13._8*q3+11._8*x3-28._8*(q3+x3)*nu+16._8*(q3+x3)*nu**2)*SIN(dip) &
        +2*y2p*(3._8-4._8*nu)*SIN(2._8*dip) &
        +(-3._8*q3-5._8*x3+4._8*(q3+x3)*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        -2*q3*x2*COS(dip) &
        +(3._8*x2*y2p+q3*y3p-x3*y3p)*COS(2._8*dip) &
        -2*y2p*y3p*COS(3._8*dip) &
        +(-x2*y2p+(q3+x3)*y3p)*COS(4._8*dip) &
        +(-q3**2+x2**2+x3**2)*SIN(dip) &
        +(3._8*q3*y2p+x3*y2p-x2*y3p)*SIN(2._8*dip) &
        -(x2**2+(q3+x3)**2+2*y2p**2)*SIN(3._8*dip) &
        +((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip) &
    )

  END FUNCTION I223

  !---------------------------------------------------------------
  !> function I222
  !! computes the I integral I222
  !---------------------------------------------------------------
  REAL*8 FUNCTION I222(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I222= &
    p3(y2p,y3p)*atan2(p2(y2p,y3p)/p3(y2p,y3p),1._8)*(-3._8+4._8*nu+COS(2._8*dip)) &
    +(1._8/4._8)*log(r12p(y2p,y3p))*( &
        -6._8*y2p &
        +8._8*y2p*nu &
        -2*(q3-x3)*COS(dip)*(-4._8+4._8*nu+COS(2*dip)) &
        -2*x2*(-2._8+4._8*nu+COS(2._8*dip))*SIN(dip) &
        +2*y3p*SIN(2._8*dip) &
	) &
	+(1._8/4._8)*log(r22p(y2p,y3p))*( &
        -2*y2p*(5._8+4._8*nu*((-3._8)+2._8*nu)) &
        +(13._8*q3+11._8*x3-28._8*(q3+x3)*nu+16._8*(q3+x3)*nu**2)*COS(dip) &
        +(-3._8*q3-5._8*x3+4._8*(q3+x3)*nu)*COS(3._8*dip) &
        +x2*(7._8+4._8*nu*(-5._8+4._8*nu))*SIN(dip) &
        +2*y3p*(3._8-4._8*nu)*SIN(2._8*dip) &
        +x2*(-3._8+4._8*nu)*SIN(3._8*dip) &
    ) &
    +(1._8/2._8)*atan2(p2p(y2p,y3p)/p3p(y2p,y3p),1._8)*( &
        -2*y3p*(5._8+4._8*nu*(-3._8+2._8*nu)) &
        -x2*(-7._8+4._8*(5._8-4._8*nu)*nu)*COS(dip) &
        -2*y3p*(-3._8+4._8*nu)*COS(2._8*dip) &
        -x2*(3._8-4._8*nu)*COS(3._8*dip) &
        -(13._8*q3+11._8*x3)*SIN(dip) &
        -4._8*(q3+x3)*nu*(-7._8+4._8*nu)*SIN(dip) &
        -(-3._8*q3-5._8*x3+4._8*(q3+x3)*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        -(-q3**2+x2**2+x3**2)*COS(dip) &
        +(-q3*y2p+x3*y2p+3*x2*y3p)*COS(2._8*dip) &
        -(x2**2+(q3+x3)**2+2*y3p**2)*COS(3._8*dip) &
        +((q3+x3)*y2p+x2*y3p)*COS(4._8*dip) &
        -2._8*q3*x2*SIN(dip) &
        +(x2*y2p+(3*q3+x3)*y3p)*SIN(2._8*dip) &
        -2._8*y2p*y3p*SIN(3._8*dip) &
        +(x2*y2p-(q3+x3)*y3p)*SIN(4._8*dip) &
	)

  END FUNCTION I222

  !---------------------------------------------------------------
  !> function I233
  !! computes the I integral I233
  !---------------------------------------------------------------
  REAL*8 FUNCTION I233(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I233= &
    4._8*y3p*(-1._8+nu)*(-1._8+2._8*nu)*atan2((x2-y2(y2p,y3p))/(x3+y3(y2p,y3p)),1._8) &
    +(-1._8/2._8)*p2(y2p,y3p)*COS(2._8*dip)*log(r12p(y2p,y3p)) &
    +4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2(p2p(y2p,y3p)*p3p(y2p,y3p)**(-1),1._8)*(x2*COS(dip)-(q3+x3)*SIN(dip)) &
    -p2(y2p,y3p)*atan2(p2(y2p,y3p)**(-1)*p3(y2p,y3p),1._8)*SIN(2._8*dip) &
    +atan2((-1)*p2p(y2p,y3p)**(-1)*p3p(y2p,y3p),1._8)*SIN(dip)*( &
        2*y2p*(3._8+(-4._8)*nu)*COS(dip) &
        +(-3._8*q3-x3+4._8*(q3+x3)*nu)*COS(2._8*dip) &
        +(-3._8+4._8*nu)*(q3-x3+x2*SIN(2._8*dip)) &
    ) &
    +0.25_8*log(r22p(y2p,y3p))*( &
        -8._8*y2p*(-1._8+nu)*(-1._8+2._8*nu) &
        +(11._8*q3+x3-4._8*(7._8*q3+3*x3)*nu+16._8*(q3+x3)*nu**2)*COS(dip) &
        +2*y2p*(-3+4*nu)*COS(2._8*dip) &
        +(3._8*q3+x3-4._8*(q3+x3)*nu)*COS(3._8*dip) &
        +x2*(5._8+4._8*nu*(-5._8+4._8*nu))*SIN(dip) &
        +x2*(3._8-4._8*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        (-q3**2+x2**2+x3**2)*COS(dip) &
        +(3._8*q3*y2p+x3*y2p-x2*y3p)*COS(2._8*dip) &
        +(-1)*(x2**2+(q3+x3)**2+2*y2p**2)*COS(3._8*dip) &
        +((q3+x3)*y2p+x2*y3p)*COS(4._8*dip) &
        +2._8*q3*x2*SIN(dip) &
        +(-3._8*x2*y2p-q3*y3p+x3*y3p)*SIN(2._8*dip) &
        +2._8*y2p*y3p*SIN(3._8*dip) &
        +(x2*y2p-(q3+x3)*y3p)*SIN(4._8*dip) &
	)

  END FUNCTION I233

  !---------------------------------------------------------------
  !> function I232
  !! computes the I integral I232
  !---------------------------------------------------------------
  REAL*8 FUNCTION I232(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I232= &
    4._8*y2p*(-1._8+nu)*(-1._8+2._8*nu)*atan2( &
        (x2-y3p*COS(dip)-y2p*SIN(dip)) &
        /(q3+x3-y2p*COS(dip)+y3p*SIN(dip)),1._8) &
    -(1._8/2._8)*p3(y2p,y3p)*COS(2*dip)*log(r12p(y2p,y3p)) &
    -4._8*(-1._8+nu)*(-1._8+2._8*nu)*atan2(p2p(y2p,y3p)**(-1)*p3p(y2p,y3p),1._8)*((q3+x3)*COS(dip)+x2*SIN(dip)) &
    +p3(y2p,y3p)*atan2(p2(y2p,y3p)*p3(y2p,y3p)**(-1),1._8)*SIN(2._8*dip) &
    -atan2(p2p(y2p,y3p)*p3p(y2p,y3p)**(-1),1._8)*COS(dip)*( &
        (3._8*q3+x3-4._8*(q3+x3)*nu)*COS(2._8*dip) &
        +(-3._8+4._8*nu)*(q3-x3+2*y3p*SIN(dip)-x2*SIN(2._8*dip)) &
    ) &
    +0.25_8*log(r22p(y2p,y3p))*( &
        8._8*y3p*(-1._8+nu)*(-1._8+2._8*nu) &
        +x2*(-5._8+4._8*(5._8-4._8*nu)*nu)*COS(dip) &
        +2._8*y3p*(-3._8+4._8*nu)*COS(2._8*dip) &
        +x2*(3._8-4._8*nu)*COS(3._8*dip) &
        +(11._8*q3+x3-4._8*(7._8*q3+3._8*x3)*nu+16._8*(q3+x3)*nu**2)*SIN(dip) &
        +(-3._8*q3-x3+4._8*(q3+x3)*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        -2._8*q3*x2*COS(dip) &
        +(x2*y2p+3*q3*y3p+x3*y3p)*COS(2._8*dip) &
        -2._8*y2p*y3p*COS(3._8*dip) &
        +(x2*y2p-(q3+x3)*y3p)*COS(4._8*dip) &
        +(-q3**2+x2**2+x3**2)*SIN(dip) &
        +(q3*y2p-x3*y2p-3*x2*y3p)*SIN(2._8*dip) &
        +(x2**2+(q3+x3)**2+2*y3p**2)*SIN(3._8*dip) &
        -((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip) &
    )

  END FUNCTION I232

  !---------------------------------------------------------------
  !> function I323
  !! computes the I integral I323
  !---------------------------------------------------------------
  REAL*8 FUNCTION I323(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I323= &
    -4._8*y3p*(-1._8+nu)*(-1._8+2._8*nu)*atan2( &
        (x2-y3p*COS(dip)-y2p*SIN(dip)) &
        /(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip)),1._8) &
    -(1._8/2._8)*p2(y2p,y3p)*COS(2._8*dip)*log(r12p(y2p,y3p)) &
    -4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2(p2p(y2p,y3p)*p3p(y2p,y3p)**(-1),1._8)*(x2*COS(dip)-(q3+x3)*SIN(dip)) &
    -p2(y2p,y3p)*atan2(p2(y2p,y3p)**(-1)*p3(y2p,y3p),1._8)*SIN(2._8*dip) &
    +atan2((-1)*p2p(y2p,y3p)**(-1)*p3p(y2p,y3p),1._8)*SIN(dip)*( &
        2._8*y2p*(3._8+(-4._8)*nu)*COS(dip) &
        +((-3._8)*q3+(-5._8)*x3 &
        +4._8*(q3+x3)*nu)*COS(2._8*dip) &
        +((-3._8)+4._8*nu)*(q3+(-1)*x3+x2*SIN(2._8*dip)) &
    ) &
    +0.25_8*log(r22p(y2p,y3p))*( &
        8._8*y2p*(-1._8+nu)*(-1._8+2._8*nu) &
        -(x3*(19._8+4._8*nu*(-9._8+4._8*nu)) &
        +q3*(5._8+4._8*nu*(-5._8+4._8*nu)))*COS(dip) &
        +2._8*y2p*(-3._8+4._8*nu)*COS(2._8*dip) &
        +(3._8*q3+5._8*x3-4._8*(q3+x3)*nu)*COS(3._8*dip) &
        +x2*(-11._8+4._8*(7._8-4._8*nu)*nu)*SIN(dip) &
        +x2*(3._8-4._8*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        (q3**2-x2**2-x3**2)*COS(dip) &
        +(-(3._8*q3+x3)*y2p+x2*y3p)*COS(2._8*dip) &
        +(x2**2+(q3+x3)**2+2*y2p**2)*COS(3._8*dip) &
        -((q3+x3)*y2p+x2*y3p)*COS(4._8*dip) &
        -2._8*q3*x2*SIN(dip) &
        +(3._8*x2*y2p+q3*y3p-x3*y3p)*SIN(2._8*dip) &
        -2._8*y2p*y3p*SIN(3._8*dip) &
        +(-x2*y2p+(q3+x3)*y3p)*SIN(4._8*dip) &
    )

  END FUNCTION I323

  !---------------------------------------------------------------
  !> function I322
  !! computes the I integral I322
  !---------------------------------------------------------------
  REAL*8 FUNCTION I322(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I322= &
    -4._8*y2p*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2( &
        (x2+(-1)*y3p*COS(dip)+(-1)*y2p*SIN(dip)) &
        *(q3+x3+(-1)*y2p*COS(dip)+y3p*SIN(dip))**(-1),1._8) &
    -(0.5_8)*p3(y2p,y3p)*COS(2*dip)*log(r12p(y2p,y3p)) &
    +4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan2( &
        p2p(y2p,y3p)**(-1)*p3p(y2p,y3p),1._8)*((q3+x3)*COS(dip)+x2*SIN(dip)) &
    +atan2(p2p(y2p,y3p)*p3p(y2p,y3p)**(-1),1._8)*COS(dip)*( &
        (-3._8*q3-5._8*x3+4._8*(q3+x3)*nu)*COS(2._8*dip) &
        -(-3._8+4._8*nu)*(q3-x3+2*(y3p-x2*COS(dip))*SIN(dip)) &
    ) &
    +p3(y2p,y3p)*atan2(p2(y2p,y3p)/p3(y2p,y3p),1._8)*SIN(2._8*dip) &
    +0.25_8*log(r22p(y2p,y3p))*( &
        -8._8*y3p*(-1._8+nu)*(-1._8+2._8*nu) &
        +x2*(11._8+4._8*nu*(-7._8+4._8*nu))*COS(dip) &
        +2._8*y3p*(-3._8+4._8*nu)*COS(2._8*dip) &
        +x2*(3._8-4._8*nu)*COS(3._8*dip) &
        -(5._8*q3+19._8*x3-4._8*(5._8*q3+9._8*x3)*nu+16._8*(q3+x3)*nu**2)*SIN(dip) &
        +(-3._8*q3-5._8*x3+4._8*(q3+x3)*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        2._8*q3*x2*COS(dip) &
        -(x2*y2p+3*q3*y3p+x3*y3p)*COS(2._8*dip) &
        +2._8*y2p*y3p*COS(3._8*dip) &
        +(-x2*y2p+(q3+x3)*y3p)*COS(4._8*dip) &
        -(-q3**2+x2**2+x3**2)*SIN(dip) &
        +(-q3*y2p+x3*y2p+3*x2*y3p)*SIN(2._8*dip) &
        -(x2**2+(q3+x3)**2+2*y3p**2)*SIN(3._8*dip) &
        +((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip) &
    )

  END FUNCTION I322

  !---------------------------------------------------------------
  !> function I333
  !! computes the I integral I333
  !---------------------------------------------------------------
  REAL*8 FUNCTION I333(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I333= &
    p2(y2p,y3p)*atan2(p3(y2p,y3p)/p2(y2p,y3p),1._8)*((-3._8)+4._8*nu+COS(2._8*dip)) &
    +(1._8/4._8)*log(r12p(y2p,y3p))*( &
        -6._8*y3p+8._8*y3p*nu &
        +x2*(7._8-8._8*nu)*COS(dip) &
        -x2*COS(3._8*dip) &
        +2._8*(q3-x3)*(-2._8+4._8*nu+COS(2._8*dip))*SIN(dip) &
        -2._8*y2p*SIN(2._8*dip) &
    ) &
    +0.5_8*atan2(-p3p(y2p,y3p)/p2p(y2p,y3p),1._8)*( &
        2._8*y2p*(5._8+4._8*nu*(-3._8+2._8*nu)) &
        -(7._8*q3+5._8*x3-20._8*(q3+x3)*nu+16._8*(q3+x3)*nu**2)*COS(dip) &
        +2._8*y2p*(-3._8+4._8*nu)*COS(2._8*dip) &
        +(3._8*q3+x3-4._8*(q3+x3)*nu)*COS(3._8*dip) &
        +x2*(-13._8+4._8*(7._8-4._8*nu)*nu)*SIN(dip) &
        +x2*(3._8-4._8*nu)*SIN(3._8*dip) &
    ) &
    +0.25_8*log(r22p(y2p,y3p))*( &
        -2._8*y3p*(5._8+4._8*nu*(-3._8+2._8*nu)) &
        +x2*(13._8+4._8*nu*(-7._8+4._8*nu))*COS(dip) &
        +x2*(-3._8+4._8*nu)*COS(3._8*dip) &
        -(7._8*q3+5._8*x3-20._8*(q3+x3)*nu+16._8*(q3+x3)*nu**2)*SIN(dip) &
        +2._8*y2p*(-3._8+4._8*nu)*SIN(2._8*dip) &
        +(3._8*q3+x3-4._8*(q3+x3)*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        -2._8*q3*x2*COS(dip) &
        +(3._8*x2*y2p+q3*y3p-x3*y3p)*COS(2._8*dip) &
        -2._8*y2p*y3p*COS(3._8*dip) &
        +(-x2*y2p+(q3+x3)*y3p)*COS(4._8*dip)+(-q3**2+x2**2+x3**2)*SIN(dip) &
        +(3._8*q3*y2p+x3*y2p-x2*y3p)*SIN(2._8*dip) &
        -(x2**2+(q3+x3)**2+2*y2p**2)*SIN(3._8*dip) &
        +((q3+x3)*y2p+x2*y3p)*SIN(4._8*dip) &
    )

  END FUNCTION I333

  !---------------------------------------------------------------
  !> function I332
  !! computes the I integral I332
  !---------------------------------------------------------------
  REAL*8 FUNCTION I332(y2p,y3p)
    REAL*8, INTENT(IN) :: y2p,y3p

    I332= &
    -p3(y2p,y3p)*atan2(p2(y2p,y3p)/p3(y2p,y3p),1._8)*(3._8-4._8*nu+COS(2._8*dip)) &
    +(1._8/4._8)*log(r12p(y2p,y3p))*( &
        -6._8*y2p &
        +8._8*y2p*nu &
        +2._8*(q3-x3)*COS(dip)*(2._8-4._8*nu+COS(2._8*dip)) &
        +2._8*x2*(4._8-4._8*nu+COS(2._8*dip))*SIN(dip) &
        -2._8*y3p*SIN(2._8*dip) &
    ) &
    +0.25_8*log(r22p(y2p,y3p))*( &
        (-2._8)*y2p*(5._8+4._8*nu*((-3._8)+2._8*nu)) &
        +(7._8*q3+5._8*x3+(-20._8)*(q3+x3)*nu+16._8*(q3+x3)*nu**2)*COS(dip) &
        +(3._8*q3+x3+(-4._8)*(q3+x3)*nu)*COS(3._8*dip) &
        +x2*(13._8+4._8*nu*((-7._8)+4._8*nu))*SIN(dip) &
        +2._8*y3p*((-3._8)+4._8*nu)*SIN(2._8*dip) &
        +x2*(3._8+(-4._8)*nu)*SIN(3._8*dip) &
    ) &
    +0.5_8*atan2(p2p(y2p,y3p)/p3p(y2p,y3p),1._8)*( &
        -2._8*y3p*(5._8+4._8*nu*(-3._8+2._8*nu)) &
        -x2*(-13._8+4._8*(7._8-4._8*nu)*nu)*COS(dip) &
        -2._8*y3p*(3._8-4._8*nu)*COS(2._8*dip) &
        -x2*((-3._8)+4._8*nu)*COS(3._8*dip) &
        -(7._8*q3+5._8*x3)*SIN(dip) &
        -4._8*(q3+x3)*nu*((-5._8)+4._8*nu)*SIN(dip) &
        -(3._8*q3+x3+(-4._8)*(q3+x3)*nu)*SIN(3._8*dip) &
    ) &
    +r2p(y2p,y3p)**(-2)*x3*( &
        (-1)*((-1)*q3**2+x2**2+x3**2)*COS(dip) &
        +((-1)*q3*y2p+x3*y2p+3*x2*y3p)*COS(2._8*dip) &
        -(x2**2+(q3+x3)**2+2*y3p**2)*COS(3._8*dip) &
        +((q3+x3)*y2p+x2*y3p)*COS(4._8*dip) &
        +(-2._8)*q3*x2*SIN(dip) &
        +(x2*y2p+(3*q3+x3)*y3p)*SIN(2._8*dip) &
        -2._8*y2p*y3p*SIN(3._8*dip) &
        +(x2*y2p+(-1)*(q3+x3)*y3p)*SIN(4._8*dip) &
	)

  END FUNCTION I332

END








