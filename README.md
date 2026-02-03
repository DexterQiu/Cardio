# Cardio

## ANSYS Transient Structural command snippet (direct solver)

If you see a warning that an iterative solver was used, you can force a
direct solver in Mechanical APDL with `EQSLV,SPARSE` before `SOLVE`:

```apdl
/SOLU
ANTYPE,TRANS
TRNOPT,FULL
EQSLV,SPARSE
SOLVE
FINISH
```

## ANSYS Mechanical APDL command snippet (Transient Structural + USERMAT)

This snippet is suitable for a Mechanical "Commands (APDL)" object or a
standalone APDL input. It defines a SOLID186 element, sets a USERMAT
material with three properties (mu, D1, ramp time), and enforces a direct
solver.

```apdl
/prep7
! Define element type (example: SOLID186 for large deformation)
! NOTE: Large deformation is controlled by NLGEOM,ON (not SOLID186 KEYOPT(3)).
! If your Mechanical GUI already inserted a KEYOPT(3)=2 line, remove it.
et,1,186

! Define material with USERMAT
mp,ex,1,1.0       ! Dummy (not used by USERMAT)
mp,nuxy,1,0.3     ! Dummy (not used by USERMAT)

! USERMAT definition: TB,USER
! TB,USER,mat_id,1,number_of_properties
! For PROPS: mu, D1, ramp_time
TB,USER,1,1,3
TBDATA,1,1e6,1.0e-9,0.75

! Geometry / mesh assumed defined here
EQSLV,SPARSE

finish

/solu
trnopt,full       ! Full transient
nlgeom,on         ! Large deformation
kbc,1             ! Ramped loads
```
