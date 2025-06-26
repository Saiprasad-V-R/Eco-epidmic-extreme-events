# Predator–Prey Statistical Analysis

This project contains a Fortran program that simulates a predator–prey system with sinusoidal forcing and computes two outputs:

1. A time series of all three state variables (`ts_0.8.dat`).  
2. A single “hs” statistic per run (`hs_0.8.dat`), combining peak mean and variability.

You can then import these data files into MATLAB’s **Distribution Fitter** app (or use built-in functions) to fit probability distributions and generate publication-quality plots.

---

## Prerequisites

- **Fortran compiler** (e.g. `gfortran`, `ifort`, `pgfortran`)  
- **MATLAB** (R2016b or later recommended, for the Distribution Fitter app)

---

## 1. Compile & Run

1. **Place** your Fortran source (`PredatorPrey.f90`) in a working directory.  
2. **Compile**:
   ```bash
   gfortran -O2 PredatorPrey.f90 -o PredatorPrey
