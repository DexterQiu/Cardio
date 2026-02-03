! ---------------------------------------------------------------------------
! ANSYS USERMAT: Large-deformation Neo-Hookean (compressible penalty)
! ---------------------------------------------------------------------------
! This implementation follows the common USERMAT argument list used by ANSYS
! Mechanical APDL when TB,USER is enabled with large-deformation kinematics.
! If your ANSYS version uses a different interface, adjust the subroutine
! signature accordingly.
!
! Material parameters (PROPS):
!   PROPS(1) = shear modulus, mu (ANSYS Neo-Hookean "Initial Shear Modulus")
!   PROPS(2) = incompressibility parameter D1 (ANSYS "Incompressibility Parameter D1")
!   PROPS(3) = optional ramp time for penalty (use 0 to disable)
!
! Notes:
! - Compressible penalty formulation:
!   W = (mu/2) * (I1 - 3) + (1/D1) * (J - 1)^2
!
! Voigt order (ANSYS):
!   STRESS(1:6) = [Sxx, Syy, Szz, Sxy, Syz, Sxz] (Cauchy)
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
      DOUBLE PRECISION F(3,3),CMAT(3,3),I3(3,3)
      DOUBLE PRECISION S2PK(3,3),CAUCHY(3,3)
      DOUBLE PRECISION MU,D1,D1_EFF,J,TRC,PSI,RAMP
      DOUBLE PRECISION CWORK(3,3),SPERT(3,3),SMINUS(3,3)
      DOUBLE PRECISION EPS,DC,DE
      DOUBLE PRECISION MAT_DET3
      EXTERNAL MAT_INV3
      INTEGER JIDX

!     Initialize identity
      CALL MAT_IDENTITY(I3)

!     Deformation gradient at end of step
      F = DFGRD1

!     Right Cauchy-Green tensor: C = F^T * F
      CALL MAT_MATMULT(TRANSPOSE(F), F, CMAT)

!     Determinant of F
      J = MAT_DET3(F)

!     Material properties
      MU = PROPS(1)
      IF (NPROPS .GE. 2) THEN
        D1 = PROPS(2)
      ELSE
        D1 = 1.0D-6 / MAX(MU, 1.0D-12)
      END IF
      RAMP = 1.0D0
      IF (NPROPS .GE. 3) THEN
        IF (PROPS(3) .GT. 0.0D0) THEN
          RAMP = MIN(1.0D0, TIME(1) / PROPS(3))
        END IF
      END IF
      D1_EFF = D1 / MAX(RAMP, 1.0D-12)

!     2nd PK stress (Total Lagrangian, penalty)
!     S = mu * I + (2/D1_eff) * (J - 1) * J * Cinv
      CALL NEOHOOK_S2PK(CMAT, MU, D1_EFF, S2PK, TRC)

!     Strain energy density
      PSI = 0.5D0*MU*(TRC - 3.0D0) + (1.0D0/D1_EFF) * (J - 1.0D0)**2
      SSE = PSI

!     Cauchy stress: sigma = (1/J) * F * S * F^T
      CALL MAT_MATMULT(F, S2PK, CWORK)
      CALL MAT_MATMULT(CWORK, TRANSPOSE(F), CAUCHY)
      CAUCHY = CAUCHY / J

!     Fill STRESS (Voigt) with Cauchy stress
      STRESS(1) = CAUCHY(1,1)
      STRESS(2) = CAUCHY(2,2)
      STRESS(3) = CAUCHY(3,3)
      STRESS(4) = CAUCHY(1,2)
      STRESS(5) = CAUCHY(2,3)
      STRESS(6) = CAUCHY(1,3)

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
        CALL NEOHOOK_S2PK(CWORK, MU, D1_EFF, SPERT, TRC)

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
        CALL NEOHOOK_S2PK(CWORK, MU, D1_EFF, SMINUS, TRC)

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
      SUBROUTINE NEOHOOK_S2PK(CMAT, MU, D1, S2PK, TRC)
      IMPLICIT NONE
      DOUBLE PRECISION CMAT(3,3), S2PK(3,3)
      DOUBLE PRECISION MU, D1, TRC
      DOUBLE PRECISION CMATINV(3,3), I3(3,3), J
      DOUBLE PRECISION MAT_DET3
      EXTERNAL MAT_INV3

      CALL MAT_IDENTITY(I3)
      J = DSQRT(MAT_DET3(CMAT))
      TRC = CMAT(1,1) + CMAT(2,2) + CMAT(3,3)
      CALL MAT_INV3(CMAT, CMATINV)
      S2PK = MU * I3 + (2.0D0/D1) * (J - 1.0D0) * J * CMATINV
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
