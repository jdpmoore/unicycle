ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Printing routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
        module PRINI_MOD
        implicit none

        integer, save :: IP = 0
        integer, save :: IQ = 0

        end module PRINI_MOD




        SUBROUTINE PRINI(IP1,IQ1)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: IP1
        integer, intent(in) :: IQ1

        IP=IP1
        IQ=IQ1

        RETURN
        END SUBROUTINE PRINI

C
C
C
C
C
        SUBROUTINE PRINC(MES)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)

        CALL  MESSPR(MES,IP,IQ)

        RETURN
        END SUBROUTINE PRINC
C
C
C
C
C
        SUBROUTINE PRIN(MES,A,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        real(4), intent(in) :: A(N)

        integer :: J

        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E12.5))
        RETURN
        END SUBROUTINE PRIN
C
C
C
C
C
        SUBROUTINE PRINS(MES,A)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        real(4), intent(in) :: A

        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0) WRITE(IP,1200)A
        IF(IQ.NE.0) WRITE(IQ,1200)A
 1200 FORMAT(6(2X,E12.5))
        RETURN
        END SUBROUTINE PRINS
C
C
C
C
        SUBROUTINE PRIN2(MES,A2,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        real(8), intent(in) :: A2(N)

        integer :: J

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E12.5))
        RETURN
        END SUBROUTINE PRIN2
C
C
C
C
        SUBROUTINE PRIN2S(MES,A2)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        real(8), intent(in) :: A2

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0) WRITE(IP,1400)A2
        IF(IQ.NE.0) WRITE(IQ,1400)A2
 1400 FORMAT(6(2X,E12.5))
        RETURN
        END SUBROUTINE PRIN2S
C
C
C
C
        SUBROUTINE PRINQ(MES,A4,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        real(8), intent(in) :: A4(N)

        integer :: J

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e12.5))
        RETURN
        END SUBROUTINE PRINQ
C
C
C
C
        SUBROUTINE PRINQS(MES,A4)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        real(8), intent(in) :: A4

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0) WRITE(IP,1500)A4
        IF(IQ.NE.0) WRITE(IQ,1500)A4
 1500 FORMAT(6(2X,e12.5))
        RETURN
        END SUBROUTINE PRINQS
C
C
C
C
        SUBROUTINE PRINF(MES,IA,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        integer, intent(in) :: IA(N)

        integer :: J

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
        END SUBROUTINE PRINF
C
C
C
C
        SUBROUTINE PRINFS(MES,IA)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        integer, intent(in) :: IA

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0) WRITE(IP,1600)IA
        IF(IQ.NE.0) WRITE(IQ,1600)IA
 1600 FORMAT(10(1X,I7))
        RETURN
        END SUBROUTINE PRINFS
C
C
C
C
        SUBROUTINE PRINF1(MES,IA1,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        integer(4), intent(in) :: IA1(N)

        integer :: J

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA1(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA1(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
        END SUBROUTINE PRINF1
C
C
C
C
        SUBROUTINE PRINF1S(MES,IA1)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        integer(4), intent(in) :: IA1

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0) WRITE(IP,1600)IA1
        IF(IQ.NE.0) WRITE(IQ,1600)IA1
 1600 FORMAT(10(1X,I7))
        RETURN
        END SUBROUTINE PRINF1S
C
C
C
C
        SUBROUTINE PRINF2(MES,IA2,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        integer(2), intent(in) :: IA2(N)

        integer :: J

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
        END SUBROUTINE PRINF2
C
C
C
C
        SUBROUTINE PRINF2S(MES,IA2)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        integer(2), intent(in) :: IA2

        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0) WRITE(IP,1600)IA2
        IF(IQ.NE.0) WRITE(IQ,1600)IA2
 1600 FORMAT(10(1X,I7))
        RETURN
        END SUBROUTINE PRINF2S
C
C
C
C
        SUBROUTINE PRINA(MES,AA,N)

        use PRINI_MOD
        implicit none

        integer, intent(in) :: N
        character(1), intent(in) :: MES(*)
        character(1), intent(in) :: AA(N)

        integer :: J

        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END SUBROUTINE PRINA
C
C
C
C
        SUBROUTINE PRINAS(MES,AA)

        use PRINI_MOD
        implicit none

        character(1), intent(in) :: MES(*)
        character(1), intent(in) :: AA

        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0) WRITE(IP,2000)AA
        IF(IQ.NE.0) WRITE(IQ,2000)AA
        RETURN
        END SUBROUTINE PRINAS
c
c
c
c
c
        SUBROUTINE MESSPR(MES,IP,IQ)

        implicit none

        character(1), intent(in) :: MES(*)
        integer, intent(in) :: IP
        integer, intent(in) :: IQ

        character(1), parameter :: AST = '*'

        integer :: I
        integer :: I1
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END SUBROUTINE MESSPR
c
c
c
c
c
        subroutine msgmerge(a,b,c)

        implicit none

        character(1), intent(in) :: a(*)
        character(1), intent(in) :: b(*)
        character(1), intent(inout) :: c(*)

        character(1), parameter :: ast = '*'

        integer :: i
        integer :: iadd
c
        do 1200 i=1,1000
c
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)       
        iadd=i
 1200 continue
c
 1400 continue
c
        do 1800 i=1,1000
c
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c
