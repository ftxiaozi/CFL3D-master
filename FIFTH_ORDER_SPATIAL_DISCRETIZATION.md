# Implementing Fifth-Order Spatial Discretization in CFL3D-Style Code

This document describes how fifth-order spatial discretization is implemented in this project, and how to reproduce the same approach in another version of the codebase.

The implementation in this repository is a fifth-order WENO interface reconstruction for the inviscid flux terms. It is not a wholesale replacement of the solver. The existing residual assembly, flux splitting/Roe flux routines, boundary data layout, viscous flux routines, and time integration remain in place. The fifth-order change is localized to the construction of left and right states at each cell face.

## 1. Where the Current Implementation Lives

The relevant files are:

- `source/cfl3d/libs/xlim.F`
  - Contains `xlim`, the older limiter-based second/third-order reconstruction.
  - Contains `xlim2`, the WENO-5 reconstruction routine.

- `source/cfl3d/libs/ffluxr.F`
  - I-direction inviscid residual contribution.
  - Uses `xlim2` when `xkap >= 2.0`.

- `source/cfl3d/libs/gfluxr.F`
  - J-direction inviscid residual contribution.
  - Uses the same WENO-5 pattern as `ffluxr.F`.

- `source/cfl3d/libs/hfluxr.F`
  - K-direction inviscid residual contribution.
  - Uses the same WENO-5 pattern as `ffluxr.F`.

- `source/cfl3d/dist/global.F`
  - Reads `IFLIM` and `RKAP0` from the CFL3D input file.

- `source/cfl3d/dist/resid.F`
  - Copies `RKAP0` into directional `rkap`.
  - Passes `rkap(1:3)` into `ffluxr`, `gfluxr`, and `hfluxr`.

The user input file enables this mode with:

```text
IFDS(I) IFDS(J) IFDS(K) RKAP0(I) RKAP0(J) RKAP0(K)
      1       1       1   2.0000   2.0000   2.0000
```

In this project, `RKAP0 = 2.0` is the practical switch that selects the WENO-5 branch.

## 2. High-Level Algorithm

For each coordinate direction and for each conservative/primitive variable component `l = 1..5`, the code does the following:

1. Copy cell-centered states into two temporary arrays:
   - `t(:,20+l)` stores one side of the reconstructed interface state.
   - `t(:,25+l)` stores the other side of the reconstructed interface state.

2. Compute first differences between neighboring cell values.

3. Assemble four neighboring differences around each target face:
   - `del(-2)`
   - `del(-1)`
   - `del(+1)`
   - `del(+2)`

4. Call `xlim2` to compute WENO-5 interface correction values.

5. Add/subtract those corrections from the cell-centered state to get left and right face states.

6. Pass the reconstructed left/right states to the existing inviscid flux routine:
   - `fluxp`/`fluxm` for van Leer flux-vector splitting when `IFDS = 0`.
   - `fhat` for Roe flux-difference splitting when `IFDS = 1`.
   - `fmaps` for MAPS+ when selected.

7. Accumulate residuals by taking flux differences between adjacent faces.

The important design point is that WENO-5 only changes how interface states are reconstructed. The flux calculation itself is reused.

## 3. Switch Logic

The current code uses `xkap` to select between the legacy reconstruction and WENO-5.

In `ffluxr.F`, the structure is:

```fortran
      if (real(xkap).lt.+2.e0) then
         ...
         call xlim(...)
         ...
      else
         ...
         call xlim2(...)
         ...
      end if
```

The same structure appears in `gfluxr.F` and `hfluxr.F`.

To reproduce this in another version:

1. Preserve the existing reconstruction path for `xkap < 2.0`.
2. Add a new branch for `xkap >= 2.0`.
3. In that new branch, call a WENO-5 reconstruction routine equivalent to `xlim2`.
4. Ensure the input parser can set `RKAP0(I/J/K) = 2.0`.
5. Ensure `resid.F` or the equivalent residual driver passes the directional `rkap` into each directional inviscid flux routine.

## 4. Data Meaning in the WENO-5 Routine

The current WENO-5 routine is:

```fortran
      subroutine xlim2(n,x1,x2,x3,x4)
```

Its arguments have this meaning:

```text
Input:
  x1 = del(-2)
  x2 = del(-1)
  x3 = del(+1)
  x4 = del(+2)

Output:
  x2 = correction for right-interface construction
  x1 = correction for left-interface construction
```

The code reconstructs corrections rather than absolute states. The calling flux routine then applies those corrections as:

```fortran
right_state = cell_value + x2
left_state  = cell_value - x1
```

The exact array offsets differ by direction, but the concept is the same.

## 5. Interior I-Direction Template Assembly

The I-direction implementation in `ffluxr.F` is the clearest model.

First, the code computes nearest-neighbor differences:

```fortran
      do izz=1,nr-jv
         t(izz+jv,1) = t(izz+jv,25+l)-t(izz,25+l)
      end do
```

Then it builds the four shifted difference arrays:

```fortran
      do izz=1,nr
         t(izz,2) = t(izz+jv,1)
      end do

      do izz=1,nr2
         t(izz,3) = t(izz+2*jv,1)
      end do

      do izz=1,nr3
         t(izz,4) = t(izz+3*jv,1)
      end do
```

At the call:

```fortran
      call xlim2(nr3,t(1,1),t(1,2),t(1,3),t(1,4))
```

the four arguments correspond to:

```text
t(:,1) = del(-2)
t(:,2) = del(-1)
t(:,3) = del(+1)
t(:,4) = del(+2)
```

After `xlim2`, the reconstructed corrections are applied:

```fortran
      t(izz+2*jv,20+l) = t(izz+2*jv,20+l)+t(izz,2)
      t(izz+jv,  25+l) = t(izz+jv,  25+l)-t(izz,1)
```

This produces the two states at each interface.

## 6. J- and K-Direction Assembly

The J- and K-direction implementations are structurally identical, but the stride changes.

For J-direction in `gfluxr.F`, the neighbor stride is `1` in the packed temporary array:

```fortran
      t(izz+1,1) = q(izz+1,1,i,l)-q(izz,1,i,l)
      ...
      call xlim2(nr3,t(1,1),t(1,2),t(1,3),t(1,4))
```

For K-direction in `hfluxr.F`, the neighbor stride is `jdim`:

```fortran
      t(izz+jdim,1) = q(izz+jdim,1,i,l)-q(izz,1,i,l)
      ...
      call xlim2(nr3,t(1,1),t(1,2),t(1,3),t(1,4))
```

When porting, do not copy the I-direction indexing blindly. Identify the packed-array stride in the target direction and build the four shifted differences consistently.

## 7. WENO-5 Reconstruction Formula

The current `xlim2` implementation uses three candidate substencils and nonlinear WENO-Z style weights.

For each point:

```fortran
      t1 = x1(izz)
      t2 = x2(izz)
      t3 = x3(izz)
      t4 = x4(izz)
```

The three left-interface candidate corrections are:

```fortran
      ql0 =  5.0/6*t2 - 1.0/3*t1
      ql1 =  1.0/3*t3 + 1.0/6*t2
      ql2 = -1.0/6*t4 + 2.0/3*t3
```

The three right-interface candidate corrections are:

```fortran
      qr0 =  2.0/3*t2 - 1.0/6*t1
      qr1 =  1.0/6*t3 + 1.0/3*t2
      qr2 = -1.0/3*t4 + 5.0/6*t3
```

The smoothness indicators are:

```fortran
      beta0 = 13.0/12*(t2-t1)**2 + 0.25*(3*t2-t1)**2
      beta1 = 13.0/12*(t3-t2)**2 + 0.25*(t3 + t2)**2
      beta2 = 13.0/12*(t4-t3)**2 + 0.25*(t4-3*t3)**2
```

The WENO-Z-like global smoothness measure is:

```fortran
      tau = abs(beta0-beta2)
```

The code converts the smoothness indicators into weight multipliers:

```fortran
      b0 = 1 + tau/(eps + beta0)
      b1 = 1 + tau/(eps + beta1)
      b2 = 1 + tau/(eps + beta2)
```

with:

```fortran
      eps = 1.0e-02
```

The ideal linear weights for the left state are:

```fortran
      d0 = 0.1
      d1 = 0.6
      d2 = 0.3
```

The ideal linear weights for the opposite side are reversed:

```fortran
      d0 = 0.3
      d1 = 0.6
      d2 = 0.1
```

The nonlinear weights are:

```fortran
      alpha0 = d0*b0
      alpha1 = d1*b1
      alpha2 = d2*b2

      w0 = alpha0/(alpha0+alpha1+alpha2)
      w1 = alpha1/(alpha0+alpha1+alpha2)
      w2 = alpha2/(alpha0+alpha1+alpha2)
```

The final reconstructed corrections are:

```fortran
      x2(izz) = wl0*ql0 + wl1*ql1 + wl2*ql2
      x1(izz) = wr0*qr0 + wr1*qr1 + wr2*qr2
```

The key feature of WENO is that in smooth regions the nonlinear weights approach the ideal linear weights, giving fifth-order accuracy. Near discontinuities, the weights shift away from nonsmooth substencils, reducing oscillations.

## 8. Boundary Treatment

Fifth-order reconstruction needs enough neighboring data. Near physical boundaries, block boundaries, patched boundaries, and ghost-cell regions, the interior five-point stencil is not directly available.

The current code handles this by explicitly reconstructing boundary-adjacent faces with boundary arrays:

- `qi0` for I-direction boundary/ghost data.
- `qj0` for J-direction boundary/ghost data.
- `qk0` for K-direction boundary/ghost data.

For example, near the left I-boundary, `ffluxr.F` builds:

```fortran
      t(jl,1) = boundary_or_ghost_difference_1
      t(jl,2) = boundary_or_ghost_difference_2
      t(jl,3) = interior_difference_1
      t(jl,4) = interior_difference_2
      call xlim2(...)
```

It does this separately for:

- First inner cell near the left boundary.
- First ghost/interface cell near the left boundary.
- First inner cell near the right boundary.
- First ghost/interface cell near the right boundary.

When porting, implement the interior WENO path first, but do not consider the feature complete until all boundary-adjacent interfaces have matching one-sided/ghost-supported WENO calls. Without these, the scheme may silently fall back to lower order near boundaries or access invalid data.

## 9. Positivity and Safety Checks

After reconstruction, `ffluxr.F` checks density and pressure when `ichk == 1`:

```fortran
      if (real(t(kc,21)).lt.real(epsz) .or.
     .    real(t(kc,25)).lt.real(epsz) .or.
     .    real(t(kc,21)).gt.real(epss) .or.
     .    real(t(kc,25)).gt.real(epss)) then
          stop
      end if
```

The older `xlim` path also has optional calls to `prolim` and `prolim2` for density/pressure protection. The `xlim2` path does not call `prolim` directly in this file, so the reconstructed states rely on the surrounding checks and the robustness of the WENO weights.

When porting to another version, keep or add:

- Density positivity checks.
- Pressure positivity checks.
- Optional clipping or fallback to lower-order reconstruction if positivity fails.
- Diagnostics that print block, index, and reconstructed values.

For production robustness, a common strategy is:

1. Try WENO-5.
2. If density or pressure becomes nonphysical at a face, retry with a lower-order limited reconstruction.
3. If that still fails, fall back to first order for that face.

The current repository primarily checks and stops; a more robust port can add local fallback.

## 10. Interaction With Flux Splitting

The fifth-order reconstruction is independent of the flux function. It supplies interface states.

After left and right states are reconstructed:

- If `IFDS = 0`, the code computes van Leer flux-vector splitting:

```fortran
      call fluxp(..., left_state,  ...)
      call fluxm(..., right_state, ...)
```

- If `IFDS = 1`, the code computes Roe flux-difference splitting:

```fortran
      call fhat(..., flux, right_state, left_state, ...)
```

Therefore, when implementing the fifth-order scheme in another version, avoid mixing the reconstruction change with a flux-function rewrite. First preserve the existing flux routine and feed it higher-order left/right states.

## 11. Residual Assembly

After fluxes are computed at interfaces, the residual update remains conservative:

```fortran
      res(cell,l) = res(cell,l) + flux(right_face,l) - flux(left_face,l)
```

The fifth-order method must preserve this structure. Do not reconstruct residuals directly. Reconstruct states, compute face fluxes, then difference the face fluxes.

## 12. Implementation Checklist for Another Version

Use this checklist when porting the method.

1. Add a WENO-5 reconstruction routine equivalent to `xlim2`.

2. Add a scheme-selection switch.
   - Match this project by using `RKAP0 = 2.0`.
   - Or introduce a clearer named input such as `IORDER = 5`, while still mapping it internally to the WENO path.

3. In each inviscid directional flux routine:
   - Keep the original lower-order branch.
   - Add a WENO-5 branch.
   - Compute first differences along the active direction.
   - Build four shifted difference arrays.
   - Call the WENO routine.
   - Apply returned corrections to form left/right interface states.

4. Implement all three directions:
   - I direction: equivalent to `ffluxr.F`.
   - J direction: equivalent to `gfluxr.F`.
   - K direction: equivalent to `hfluxr.F`.

5. Implement boundary-adjacent reconstruction.
   - Use existing ghost/boundary arrays.
   - Reconstruct first inner and ghost-adjacent faces.
   - Avoid out-of-bounds stencils.

6. Preserve existing flux computation.
   - `fluxp/fluxm`, `fhat`, or the equivalent flux routine should receive reconstructed states.
   - Do not alter residual differencing.

7. Preserve parallel/block behavior.
   - Ensure ghost/halo exchange provides enough values for fifth-order stencils.
   - Ensure block-split cases reconstruct consistently across block boundaries.

8. Add positivity checks or fallback.
   - Density and pressure are the critical quantities.
   - At minimum, print enough information to locate failure.

9. Verify with progressively harder tests.
   - 1D smooth advection or vortex-like smooth case: check fifth-order convergence.
   - 1D shock/contact case: check non-oscillatory behavior.
   - 2D/3D smooth case: check directional indexing.
   - Multi-block case: check block-interface behavior.
   - Production case: compare residual history, forces, and pressure distributions.

## 13. Accuracy Expectations

This implementation gives fifth-order spatial accuracy for inviscid interface reconstruction in smooth regions. The full simulation may not show fifth-order convergence if any of the following dominate:

- Lower-order boundary conditions.
- Lower-order viscous discretization.
- Lower-order time integration.
- Shock waves or discontinuities.
- Turbulence-model source terms.
- Grid singularities, stretching, skewness, or block-interface interpolation.
- Positivity fallback to lower-order reconstruction.

For accuracy verification, use a smooth problem, fine enough grids, fixed CFL/time-discretization effects, and compare grid refinement rates.

## 14. Minimal Pseudocode

The core reconstruction can be expressed as:

```text
for each direction dir:
  for each variable l = 1..5:
    initialize left_state  = cell_state
    initialize right_state = cell_state

    compute one-dimensional differences along dir

    for each interior face with enough stencil:
      d_m2 = q[i-1] - q[i-2]
      d_m1 = q[i]   - q[i-1]
      d_p1 = q[i+1] - q[i]
      d_p2 = q[i+2] - q[i+1]

      corr_left, corr_right = WENO5(d_m2, d_m1, d_p1, d_p2)

      state_left_of_face  = q[i]   + corr_right
      state_right_of_face = q[i+1] - corr_left

    handle boundary-adjacent faces using ghost/boundary values

  compute face fluxes from reconstructed states
  update residual by conservative flux difference
```

Pay close attention to the indexing convention in the target code. The current code stores the two interface states in temporary columns `20+l` and `25+l`; another version may use named arrays such as `ql` and `qr`, which is less error-prone.

## 15. Recommended Refactoring When Porting

If the target version allows cleanup, consider isolating the WENO logic more clearly than this legacy Fortran layout:

```fortran
      call weno5_reconstruct(n,dm2,dm1,dp1,dp2,corr_l,corr_r)
```

instead of overwriting `x1` and `x2` in place. This makes the data flow clearer and reduces mistakes when adding boundary handling.

However, if the goal is to reproduce this project's behavior exactly, keep the in-place `xlim2` calling convention and port the directional flux routines as closely as possible.

## 16. Summary

To implement the fifth-order spatial discretization in the same manner as this project:

- Use `RKAP0 = 2.0` or an equivalent switch to select WENO-5.
- Add a WENO-5 routine equivalent to `xlim2`.
- In each inviscid flux direction, build a five-cell stencil through four neighboring first differences.
- Use WENO-Z-like nonlinear weights to reconstruct left/right interface corrections.
- Apply those corrections to cell-centered states.
- Feed reconstructed states into the existing flux routine.
- Preserve conservative residual differencing.
- Handle boundaries and block interfaces explicitly.
- Verify smooth-case fifth-order convergence and shock-case robustness.

