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
