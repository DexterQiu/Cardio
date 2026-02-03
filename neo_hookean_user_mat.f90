! ---------------------------------------------------------------------------
! ANSYS USERMAT: Large-deformation Neo-Hookean (incompressible / penalty)
! ---------------------------------------------------------------------------
! This implementation follows the common USERMAT argument list used by ANSYS
! Mechanical APDL when TB,USER is enabled with large-deformation kinematics.
! If your ANSYS version uses a different interface, adjust the subroutine
! signature accordingly.
!
! Material parameters (PROPS):
!   PROPS(1) = shear modulus, mu
!   PROPS(2) = augmented Lagrangian bulk parameter, kappa
!   PROPS(3) = optional ramp time for kappa (same time units as TIME)
!
! State variables (STATEV):
!   STATEV(1) = pressure-like Lagrange multiplier p
!
! Notes:
! - This is a fully incompressible (augmented Lagrangian) formulation.
!   The pressure-like variable p is stored in STATEV(1) and updated as
!   p_{n+1} = p_n + kappa * (J - 1).
! - Provide NSTATV >= 1 in your material definition.
!
! Voigt order (ANSYS):
!   STRESS(1:6) = [Sxx, Syy, Szz, Sxy, Syz, Sxz]
! ---------------------------------------------------------------------------
      SUBROUTINE USERMAT(
     &  STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     &  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC )

      IMPLICIT NONE

!     Arguments
      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS
      INTEGER NOEL,NPT,LAYER,KSPT,KSTEP,KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION SSE,SPD,SCD,RPL,DDSDDT(NTENS,NTENS)
      DOUBLE PRECISION DRPLDE(NTENS),DRPLDT
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION TIME(2),DTIME,TEMP,DTEMP,PREDEF,DPRED
      DOUBLE PRECISION PROPS(NPROPS),COORDS(3),DROT(3,3),PNEWDT,CELENT
      DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)

!     Locals
      DOUBLE PRECISION F(3,3),B(3,3),CMAT(3,3),CMATINV(3,3),I3(3,3)
      DOUBLE PRECISION S2PK(3,3)
      DOUBLE PRECISION MU,KAPPA,J,TRB,P,KRAMP
      DOUBLE PRECISION LAMBDA,I1BAR,PSI
      DOUBLE PRECISION CWORK(3,3),SPERT(3,3),SMINUS(3,3)
      DOUBLE PRECISION EPS,DC,DE
      DOUBLE PRECISION MAT_DET3
      EXTERNAL MAT_INV3
      INTEGER I,JIDX

!     Initialize identity
      CALL MAT_IDENTITY(I3)

!     Deformation gradient at end of step
      F = DFGRD1

!     Left Cauchy-Green tensor: B = F * F^T
      CALL MAT_MATMULT(F, TRANSPOSE(F), B)
!     Right Cauchy-Green tensor: C = F^T * F
      CALL MAT_MATMULT(TRANSPOSE(F), F, CMAT)

!     Determinant of F
      J = MAT_DET3(F)

!     Material properties
      MU = PROPS(1)
      IF (NPROPS .GE. 2) THEN
        KAPPA = PROPS(2)
      ELSE
        KAPPA = 1.0D6 * MU
      END IF
      KRAMP = KAPPA
      IF (NPROPS .GE. 3) THEN
        IF (PROPS(3) .GT. 0.0D0) THEN
          KRAMP = KAPPA * MIN(1.0D0, TIME(1) / PROPS(3))
        END IF
      END IF
      P = 0.0D0
      IF (NSTATV .GE. 1) THEN
        P = STATEV(1)
      END IF

!     2nd PK stress (Total Lagrangian, incompressible)
!     S = mu*J^(-2/3) * (I - (I1bar/3)*Cinv) + p * Cinv
      CALL NEOHOOK_S2PK(CMAT, MU, P, S2PK, I1BAR)

!     Strain energy density (Neo-Hookean, deviatoric only)
      PSI = 0.5D0*MU*(I1BAR - 3.0D0)
      SSE = PSI

!     Fill STRESS (Voigt) with 2nd PK stress (Total Lagrangian)
      STRESS(1) = S2PK(1,1)
      STRESS(2) = S2PK(2,2)
      STRESS(3) = S2PK(3,3)
      STRESS(4) = S2PK(1,2)
      STRESS(5) = S2PK(2,3)
      STRESS(6) = S2PK(1,3)

!     Consistent tangent via numerical differentiation in Green-Lagrange strain
      CALL ZERO_TANGENT(DDSDDE, NTENS)
      EPS = 1.0D-8
      DO JIDX = 1,6
        DC = EPS
        CWORK = CMAT
        IF (JIDX .EQ. 1) THEN
          CWORK(1,1) = CWORK(1,1) + DC
        ELSEIF (JIDX .EQ. 2) THEN
          CWORK(2,2) = CWORK(2,2) + DC
        ELSEIF (JIDX .EQ. 3) THEN
          CWORK(3,3) = CWORK(3,3) + DC
        ELSEIF (JIDX .EQ. 4) THEN
          CWORK(1,2) = CWORK(1,2) + DC
          CWORK(2,1) = CWORK(2,1) + DC
        ELSEIF (JIDX .EQ. 5) THEN
          CWORK(2,3) = CWORK(2,3) + DC
          CWORK(3,2) = CWORK(3,2) + DC
        ELSE
          CWORK(1,3) = CWORK(1,3) + DC
          CWORK(3,1) = CWORK(3,1) + DC
        END IF
        CALL NEOHOOK_S2PK(CWORK, MU, P, SPERT, I1BAR)

        CWORK = CMAT
        IF (JIDX .EQ. 1) THEN
          CWORK(1,1) = CWORK(1,1) - DC
        ELSEIF (JIDX .EQ. 2) THEN
          CWORK(2,2) = CWORK(2,2) - DC
        ELSEIF (JIDX .EQ. 3) THEN
          CWORK(3,3) = CWORK(3,3) - DC
        ELSEIF (JIDX .EQ. 4) THEN
          CWORK(1,2) = CWORK(1,2) - DC
          CWORK(2,1) = CWORK(2,1) - DC
        ELSEIF (JIDX .EQ. 5) THEN
          CWORK(2,3) = CWORK(2,3) - DC
          CWORK(3,2) = CWORK(3,2) - DC
        ELSE
          CWORK(1,3) = CWORK(1,3) - DC
          CWORK(3,1) = CWORK(3,1) - DC
        END IF
        CALL NEOHOOK_S2PK(CWORK, MU, P, SMINUS, I1BAR)

        IF (JIDX .GE. 4) THEN
          DE = DC
        ELSE
          DE = 0.5D0 * DC
        END IF
        DDSDDE(1,JIDX) = (SPERT(1,1) - SMINUS(1,1)) / (2.0D0 * DE)
        DDSDDE(2,JIDX) = (SPERT(2,2) - SMINUS(2,2)) / (2.0D0 * DE)
        DDSDDE(3,JIDX) = (SPERT(3,3) - SMINUS(3,3)) / (2.0D0 * DE)
        DDSDDE(4,JIDX) = (SPERT(1,2) - SMINUS(1,2)) / (2.0D0 * DE)
        DDSDDE(5,JIDX) = (SPERT(2,3) - SMINUS(2,3)) / (2.0D0 * DE)
        DDSDDE(6,JIDX) = (SPERT(1,3) - SMINUS(1,3)) / (2.0D0 * DE)
      END DO

!     Update pressure-like multiplier (augmented Lagrangian)
      IF (NSTATV .GE. 1) THEN
        STATEV(1) = P + KRAMP * (J - 1.0D0)
      END IF

      RETURN
      END SUBROUTINE USERMAT

! ---------------------------------------------------------------------------
      SUBROUTINE MAT_IDENTITY(I3)
      IMPLICIT NONE
      DOUBLE PRECISION I3(3,3)
      I3 = 0.0D0
      I3(1,1) = 1.0D0
      I3(2,2) = 1.0D0
      I3(3,3) = 1.0D0
      END SUBROUTINE MAT_IDENTITY

! ---------------------------------------------------------------------------
      SUBROUTINE MAT_MATMULT(A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),B(3,3),C(3,3)
      INTEGER I,J,K
      C = 0.0D0
      DO I=1,3
        DO J=1,3
          DO K=1,3
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
          END DO
        END DO
      END DO
      END SUBROUTINE MAT_MATMULT

! ---------------------------------------------------------------------------
      SUBROUTINE MAT_INV3(A, AINV)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),AINV(3,3)
      DOUBLE PRECISION DET
      DOUBLE PRECISION MAT_DET3

      DET = MAT_DET3(A)
      AINV(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2)) / DET
      AINV(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2)) / DET
      AINV(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2)) / DET
      AINV(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1)) / DET
      AINV(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1)) / DET
      AINV(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1)) / DET
      AINV(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1)) / DET
      AINV(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1)) / DET
      AINV(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1)) / DET
      END SUBROUTINE MAT_INV3

! ---------------------------------------------------------------------------
      SUBROUTINE NEOHOOK_S2PK(CMAT, MU, P, S2PK, I1BAR)
      IMPLICIT NONE
      DOUBLE PRECISION CMAT(3,3), S2PK(3,3)
      DOUBLE PRECISION MU, P, I1BAR
      DOUBLE PRECISION CMATINV(3,3), I3(3,3), J, TRC
      DOUBLE PRECISION MAT_DET3
      EXTERNAL MAT_INV3

      CALL MAT_IDENTITY(I3)
      J = DSQRT(MAT_DET3(CMAT))
      TRC = CMAT(1,1) + CMAT(2,2) + CMAT(3,3)
      I1BAR = TRC * J**(-2.0D0/3.0D0)
      CALL MAT_INV3(CMAT, CMATINV)
      S2PK = MU * J**(-2.0D0/3.0D0) * (I3 - (I1BAR/3.0D0) * CMATINV)
      S2PK = S2PK + P * CMATINV
      END SUBROUTINE NEOHOOK_S2PK

! ---------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION MAT_DET3(A)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3)
      MAT_DET3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     &          - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     &          + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      END FUNCTION MAT_DET3

! ---------------------------------------------------------------------------
      SUBROUTINE ZERO_TANGENT(DDSDDE, NTENS)
      IMPLICIT NONE
      INTEGER NTENS, I, J
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J) = 0.0D0
        END DO
      END DO
      END SUBROUTINE ZERO_TANGENT

! ---------------------------------------------------------------------------
      SUBROUTINE ISOTROPIC_TANGENT(DDSDDE, LAMBDA, MU)
      IMPLICIT NONE
      DOUBLE PRECISION DDSDDE(6,6), LAMBDA, MU

!     3D isotropic elasticity in Voigt order
      DDSDDE(1,1) = LAMBDA + 2.0D0*MU
      DDSDDE(2,2) = LAMBDA + 2.0D0*MU
      DDSDDE(3,3) = LAMBDA + 2.0D0*MU
      DDSDDE(1,2) = LAMBDA
      DDSDDE(1,3) = LAMBDA
      DDSDDE(2,1) = LAMBDA
      DDSDDE(2,3) = LAMBDA
      DDSDDE(3,1) = LAMBDA
      DDSDDE(3,2) = LAMBDA
      DDSDDE(4,4) = MU
      DDSDDE(5,5) = MU
      DDSDDE(6,6) = MU
      END SUBROUTINE ISOTROPIC_TANGENT
