# Implementing IDDES and the `isub` Keyword in CFL3D-Style Code

This document explains how this project implements IDDES-family turbulence-model behavior and the `isub` keyword, so the same approach can be ported into another version of the project.

The implementation is concentrated in the turbulence-model source terms and DES length-scale routines. It does not replace the flow solver, residual assembly, or inviscid flux discretization. The core change is that DES/DDES/IDDES replace the usual RANS turbulence length scale with a hybrid RANS/LES length scale based on wall distance, shielding functions, and a grid/subgrid scale `delta`.

## 1. Relevant Files

The main files are:

- `source/cfl3d/libs/readkey.F`
  - Reads keyword values such as `ides`, `izdes`, `isub`, `icdes`, `cdes`, `cdesw`, `cdese`, `zface`, and `iddes`.

- `source/cfl3d/libs/spalart.F`
  - Implements Spalart-Allmaras turbulence model variants.
  - Contains SA-DES, SA-DDES, SA-MDDES, SA-IDDES, EDDES, and ZDES/DES mode branches.

- `source/cfl3d/libs/twoeqn.F`
  - Implements SST and other two-equation turbulence model variants.
  - Contains SST-DES, SST-DDES, SST-MDDES, SST-IDDES, simplified IDDES, original DDES-SST, and ZDES mode branches.

- `source/cfl3d/libs/deltaDES.F`
  - Computes the DES subgrid scale `delta` for `isub = 1,2,3,4`, and provides the fallback/default max-cell-length scale.

- `source/cfl3d/libs/deltaDES_isub5.F`
  - Computes the extended `isub = 5` subgrid scale using a VTM/FKH correction based on a precomputed `VTM_data` field.

The same files also exist under `source/cfl3d/dist/` in this repository. The active build may use either `libs` or `dist` depending on the build configuration, but the local modified project appears to have the relevant newer code in `source/cfl3d/libs/`.

## 2. Keyword Overview

The DES-related keywords are read from the keyword block near the top of the CFL3D input file.

Example from this project:

```text
ides 6
isub 5
icdes 1
edvislim 100000
```

The important keywords for this document are:

- `ides`: global DES-family model selector.
- `izdes`: zonal DES-family model selector, stored per block.
- `isub`: subgrid length-scale selector.
- `icdes`: selects single or blended/adaptive DES constants.
- `cdes`: single DES constant.
- `cdesw`: DES constant for the wall/k-omega side of SST blending.
- `cdese`: DES constant for the outer/k-epsilon side of SST blending.
- `iddes`: secondary SST IDDES variant switch.
- `zface`: interface location for some ZDES modes.

## 3. Reading and Defaults

In `readkey.F`, the defaults are:

```fortran
      ides = 0
      cdes = 0.65
      isub = 0
      icdes = 0
      cdesw = 0.4
      cdese = 0.3
      iddes = 0
```

The keyword parser then reads user input:

```fortran
      else if (inpstr(lc1:lc2).eq.'ides') then
         read(inpstr(lc2:lcl),*) ides

      else if (inpstr(lc1:lc2).eq.'isub') then
         read(inpstr(lc2:lcl),*) isub

      else if (inpstr(lc1:lc2).eq.'icdes') then
         read(inpstr(lc2:lcl),*) icdes

      else if (inpstr(lc1:lc2).eq.'cdesw') then
         read(inpstr(lc2:lcl),*) cdesw

      else if (inpstr(lc1:lc2).eq.'cdese') then
         read(inpstr(lc2:lcl),*) cdese
```

When porting, preserve the defaults first. Then add parser entries for each keyword if the target version does not already have them.

## 4. `ides` and `izdes`

`ides` activates DES-family behavior globally. `izdes` activates related zonal DES behavior on a block basis.

The code prevents using both at once:

```fortran
      if (izdes.ge.1) then
        if (ides.ge.1) then
          write(iunit11,'('' izdes or ides ? '')')
          call termn8(...)
        end if
      end if
```

DES should be run in 3D. The code warns when `ides > 0` and `i2d = 1`.

For SA, the practical `ides`/`izdes` options are:

| Value | SA behavior |
|---:|---|
| `0` | DES off |
| `1` | SA-DES |
| `2` | SA-DDES |
| `3` | SA-MDDES |
| `4` | SA-IDDES |
| `5` | SA EDDES / ZDES mode-II style branch |
| `6` | SA DES mode-III branch |

For SST/two-equation models, the practical options are:

| Value | SST/two-equation behavior |
|---:|---|
| `0` | DES off |
| `1` | DES-SST |
| `2` | DDES-SST |
| `3` | MDDES-SST |
| `4` | IDDES-SST |
| `5` | Simplified IDDES-SST |
| `6` | Original DDES-SST / Xiao-DDES style branch |
| `7` | ZDES mode-III with interface `zface` |
| `8` | ZDES mode-III with Xiao-DDES style behavior |

This document focuses on IDDES, especially `ides = 4` and the closely related SST `ides = 5`.

## 5. Core IDDES Idea

The project implements IDDES by changing the turbulence length scale used in the model source terms.

For SA, the wall distance used by the SA model is replaced by a hybrid distance-like scale stored in `damp1`.

For SST/two-equation models, the turbulence length scale is replaced by a hybrid scale stored in `xlscale`.

The hybrid form blends between:

- A RANS length scale near walls and attached boundary layers.
- A DES/LES length scale based on `cdes * delta` away from walls or in separated regions.

The blending is controlled by shielding/elevating functions such as `fdn`, `fdt`, `fb`, `fe`, `fe1`, and `fe2`.

## 6. SA-IDDES Implementation

The SA-IDDES branch is in `spalart.F`:

```fortran
c  IDDES-SA
      elseif ( ides .eq. 4 .or. izdes .eq. 4 ) then
```

For each cell, the code computes local grid lengths:

```fortran
      deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
      deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
      deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
      deltamax = max(deltaj,deltak,deltai)
```

It computes wall distance:

```fortran
      dist = abs(smin(j,k,i))
```

It computes a velocity-gradient magnitude:

```fortran
      velterm = sum_of_all_9_velocity_gradient_components_squared
```

It computes an IDDES alpha variable:

```fortran
      alphan = 0.25 - dist/deltamax
      ttn = alphan**2
```

It computes the elevating function component:

```fortran
      if (alphan .ge. 0.) then
         fe1 = 2.*exp(-11.09*ttn)
      else
         fe1 = 2.*exp(-9.*ttn)
      end if
```

It computes SA eddy viscosity variables:

```fortran
      chi = turre(j,k,i)/fnu(j,k,i)
      fv1 = chi**3/(chi**3+cv1**3)
      nut = fv1*turre(j,k,i)
```

It computes DDES/IDDES shielding components:

```fortran
      rdt = nut/(sqrt(velterm)*akarman*akarman*dist*dist*re)
      fdt = 1.0 - tanh((8.0*rdt)**3)

      rdl = fnu(j,k,i)/(sqrt(velterm)*akarman*akarman*dist*dist*re)
      ctt = 1.63
      cll = 3.55
      ftn = tanh((ctt*ctt*rdt)**3)
      fln = tanh((cll*cll*rdl)**10)
      fe2 = 1. - max(ftn,fln)
```

Then:

```fortran
      fe = max(0.,fe1-1.) * fai1 * fe2
      fb = min(2.*exp(-9.*ttn),1.)
      fdn = max(1.-fdt,fb)
```

The subgrid length scale `delta` is selected by `isub`:

```fortran
      if (isub .eq. 0) then
         delta = deltamax
      elseif (isub .le. 4) then
         call deltaDES(...)
      elseif (isub .eq. 5) then
         call deltaDES_isub5(...)
      endif
```

Then the code computes a near-wall DES length scale `deltani`:

```fortran
      cw = 0.15
      deltamin = min(deltaj,deltak,deltai)

      if (isub .ne. 0) then
         deltani = cw*max(dist,deltamax)
         deltani = max(deltani,deltamin)
         deltani = min(deltani,delta)
      else
         deltani = cw*max(dist,delta)
         deltani = max(deltani,deltamin)
         deltani = min(deltani,delta)
      endif
```

Finally, SA-IDDES stores the hybrid wall-distance scale in `damp1`:

```fortran
      damp1(j,k,i) = fdn*(1.+fe)*dist
     +             + (1.-fdn)*cdes*fai1*deltani
```

Later in `spalart.F`, this modified `damp1` is used where the SA model would normally use wall distance.

## 7. SST-IDDES Implementation

The SST-IDDES branch is in `twoeqn.F`:

```fortran
c  IDDES-SST
      elseif ((ides .eq. 4 .or. ides .eq. 5).or.
     .        (izdes .eq. 4 .or. izdes .eq. 5)) then
```

The code computes local grid lengths and maximum grid length:

```fortran
      deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
      deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
      deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
      deltamax = max(deltaj,deltak,deltai)
```

It computes the SST RANS length scale:

```fortran
      ell = sqrt(turre(j,k,i,2))/(cmuc1*turre(j,k,i,1)*re)
```

Here `turre(:,:,:,1)` and `turre(:,:,:,2)` are the two turbulence variables for the active two-equation model. In SST-like use, the expression behaves as a k-omega turbulence length scale.

It computes wall distance:

```fortran
      dist = abs(smin(j,k,i))
```

It computes a shielding variable:

```fortran
      rdt = vist3d(j,k,i)/(q(j,k,i,1)*
     +      sqrt(velterm)*vk*vk*dist*dist*re)
```

The `iddes` keyword changes some SST constants:

```fortran
      if (iddes .eq. 1) then
         cdl = 20.
      else
         cdl = 8.
      end if
```

Then:

```fortran
      fdt = 1.0 - tanh((cdl*rdt)**3)
      alphan = 0.25 - dist/deltamax
      ttn = alphan**2
      fb = min(2.*exp(-9.*ttn),1.)
      fdn = max(1.-fdt,fb)
```

The subgrid length scale `delta` is selected by `isub`:

```fortran
      if (isub .eq. 0) then
         delta = deltamax
      elseif (isub .le. 4) then
         call deltaDES(...)
      elseif (isub .eq. 5) then
         call deltaDES_isub5(...)
      endif
```

The elevating function is:

```fortran
      if (alphan .ge. 0.) then
         fe1 = 2.*exp(-11.09*ttn)
      else
         fe1 = 2.*exp(-9.*ttn)
      end if

      rdl = fnu(j,k,i)/(q(j,k,i,1)*
     +      sqrt(velterm)*vk*vk*dist*dist*re)

      if (iddes .eq. 1) then
         ctt = 1.87
         cll = 5.0
      else
         ctt = 1.63
         cll = 3.55
      end if

      ftn = tanh((ctt*ctt*rdt)**3)
      fln = tanh((cll*cll*rdl)**10)
      fe2 = 1. - max(ftn,fln)
```

For full SST-IDDES:

```fortran
      fe = fe2*max(0.,fe1-1.)
```

For simplified IDDES, selected by `ides = 5` or `izdes = 5`:

```fortran
      fe = 0.
```

The near-wall DES scale `deltan` is then computed:

```fortran
      cw = 0.15

      if (iddes .eq. 1) then
         deltan = cw*max(dist,delta)
         deltan = min(deltan,delta)
      else
         deltamin = min(deltaj,deltak,deltai)
         if (isub .ne. 0) then
            deltan = cw*max(dist,deltamax)
            deltan = max(deltan,deltamin)
            deltan = min(deltan,delta)
         else
            deltan = cw*max(dist,delta)
            deltan = max(deltan,deltamin)
            deltan = min(deltan,delta)
         endif
      endif
```

The local DES constant is selected:

```fortran
      if (icdes .eq. 1) then
         cdesn = cdesw*blend(j,k,i) + cdese*(1.-blend(j,k,i))
      else
         cdesn = cdes
      end if
```

Finally, SST-IDDES stores the hybrid turbulence length scale in `xlscale`:

```fortran
      xlscale(j,k,i) = fdn*(1.+fe)*ell
     +               + (1.-fdn)*cdesn*deltan
```

Later, the turbulence kinetic energy destruction term uses `xlscale`:

```fortran
      if (ides.gt.0 .or. izdes.ge.1) then
         dk = (turre(j,k,i,2)**1.5)/xlscale(j,k,i)
      end if
```

This is the key SST-IDDES coupling point.

## 8. The `isub` Keyword

`isub` selects how the DES subgrid scale `delta` is computed. It is meaningful when `ides > 0` or `izdes > 0`.

The available options in this project are:

| `isub` | Meaning |
|---:|---|
| `0` | Use maximum local cell length scale, approximately `max(deltai,deltaj,deltak)`. This is the default. |
| `1` | Use cube root of cell volume, `vol**(1/3)`. |
| `2` | Use a vorticity-direction-weighted length scale based on face-normal metric lengths. |
| `3` | Use Spalart 2015 geometry-based length scale. |
| `4` | Use Spalart 2015 geometry-based length scale plus VTM/FKH correction. |
| `5` | Use `deltaDES_isub5`, an extended VTM/FKH corrected length scale with precomputed and spatially averaged `VTM_data`. |

## 9. `isub = 0`: Default Max Cell Length

When `isub = 0`, the IDDES branches usually avoid `deltaDES` and use:

```fortran
      delta = deltamax
```

where:

```fortran
      deltaj = 2.*vol/(sj_minus_plus_face_area_sum)
      deltak = 2.*vol/(sk_minus_plus_face_area_sum)
      deltai = 2.*vol/(si_minus_plus_face_area_sum)
      deltamax = max(deltaj,deltak,deltai)
```

This is the simplest DES length scale.

## 10. `isub = 1`: Cube Root Volume Scale

In `deltaDES.F`:

```fortran
      if (isub .eq. 1) then
         delta = (vol(j,k,i))**(1./3.)
      end if
```

This is a simple isotropic cell-size estimate.

The turbulence routines also print:

```fortran
      subgrid scale is cube root of the cell
```

when DES/ZDES is active and `isub = 1`.

## 11. `isub = 2`: Vorticity-Direction Weighted Scale

For `isub = 2`, `deltaDES.F` computes face-normal length scales:

```fortran
      deltaj = 2.*vol/(sj face area sum)
      deltak = 2.*vol/(sk face area sum)
      deltai = 2.*vol/(si face area sum)
```

It computes vorticity components from velocity gradients:

```fortran
      wvi = ux_zy_difference = ux(j,k,i,8) - ux(j,k,i,6)
      wvj = ux_xz_difference = ux(j,k,i,3) - ux(j,k,i,7)
      wvk = ux_yx_difference = ux(j,k,i,4) - ux(j,k,i,2)
```

Then it normalizes the vorticity direction and forms:

```fortran
      delta = wni*deltaj*deltak
     +      + wnj*deltak*deltai
     +      + wnk*deltai*deltaj
```

This makes `delta` depend on the orientation of the resolved rotation/vorticity relative to the grid.

## 12. `isub = 3`: Spalart 2015 Geometry Scale

For `isub = 3`, `deltaDES.F` computes a geometry-based length scale inspired by:

```text
Spalart, Flow Turbulence Combustion (2015) 95:709-737
```

The code:

1. Computes the cell center from the eight vertices.
2. Builds vectors from the cell center to each vertex.
3. Computes vorticity direction.
4. Projects vertex-center vectors using cross products with the vorticity direction.
5. Finds the maximum distance among the projected vectors.
6. Sets:

```fortran
      delta = centermax/SQRT(3.)
```

This better accounts for anisotropic and skewed cells than a simple cube-root volume scale.

## 13. `isub = 4`: Spalart 2015 Scale With VTM/FKH Correction

For `isub = 4`, the code first computes the same Spalart 2015 geometry scale as `isub = 3`.

Then it computes a VTM quantity from strain and vorticity tensors, applies a limiter, and forms:

```fortran
      FKHVTM = max(FKHmin,
     +             min(FKHmax,
     +                 FKHmin + ((FKHmax-FKHmin)/(FKHa2-FKHa1))
     +                *(VTM-FKHa1)))
```

with:

```fortran
      FKHmax = 1.0
      FKHmin = 0.1
      FKHa1 = 0.15
      FKHa2 = 0.3
```

Then:

```fortran
      delta = delta*FKHVTM
```

There is additional mode-dependent protection:

```fortran
      if (ides .eq. 2) then
         if (fd .lt. 1.0-epsilon) FKHVTM = 1.0
      elseif (ides .eq. 4) then
         if (fd .gt. epsilon) FKHVTM = 1.0
      endif
```

This means the VTM/FKH correction is disabled in certain RANS/DES shielding regions, depending on the active DES mode.

## 14. `isub = 5`: Extended VTM/FKH Scale

For `isub = 5`, the code does not use `deltaDES.F`. It calls:

```fortran
      call deltaDES_isub5(...)
```

The SA and SST routines allocate a `VTM_data` array before the DES length-scale loop:

```fortran
      if (isub .eq. 5) then
         allocate(VTM_data(jdim-1,kdim-1,idim-1))
         ...
         VTM_data(j,k,i) = VTM
      end if
```

Then `deltaDES_isub5.F`:

1. Averages `VTM_data` around the local cell. The averaging stencil changes near boundaries.
2. Computes the same Spalart 2015-style geometry length scale:

```fortran
      delta = centermax/SQRT(3.)
```

3. Computes `FKHVTM` using averaged `VTM_AVE`.
4. Applies the same mode-dependent protections for `ides = 2` and `ides = 4`.
5. Returns:

```fortran
      delta = delta*FKHVTM
```

At the end of the turbulence routine, the code deallocates:

```fortran
      if (isub .eq. 5) then
         deallocate(VTM_data)
      end if
```

Your current input uses:

```text
isub 5
```

so this project uses the extended VTM/FKH-corrected DES grid length scale.

## 15. `icdes` Interaction

The `icdes` keyword controls how `CDES` is applied, especially in SST.

For `icdes = 0`:

```fortran
      cdesn = cdes
```

For `icdes = 1`, SST DES/IDDES uses:

```fortran
      cdesn = cdesw*blend(j,k,i) + cdese*(1.-blend(j,k,i))
```

where `blend` is the SST blending function.

Defaults are:

```fortran
      cdesw = 0.4
      cdese = 0.3
```

In SA code there are `icdes = 2` adaptive-constant blocks, but in the active SA-IDDES formula the final length scale still uses `cdes`:

```fortran
      damp1 = ... + (1.-fdn)*cdes*fai1*deltani
```

So treat `icdes = 2` as experimental/incomplete unless the target version deliberately completes and verifies it.

## 16. Porting Checklist

To implement IDDES and `isub` in another project version:

1. Add or verify keyword parsing:
   - `ides`
   - `izdes`
   - `isub`
   - `icdes`
   - `cdes`
   - `cdesw`
   - `cdese`
   - `iddes`
   - `zface`

2. Preserve default values:
   - `ides = 0`
   - `isub = 0`
   - `icdes = 0`
   - `cdes = 0.65`
   - `cdesw = 0.4`
   - `cdese = 0.3`

3. Add consistency checks:
   - Do not allow both `ides >= 1` and `izdes >= 1`.
   - Warn or stop for 2D DES usage, depending on project policy.
   - Ensure ZDES is active consistently across required directions/blocks.

4. Implement `deltaDES`:
   - `isub = 1`: cube-root volume.
   - `isub = 2`: vorticity-direction weighted scale.
   - `isub = 3`: Spalart 2015 geometry scale.
   - `isub = 4`: Spalart 2015 geometry scale plus VTM/FKH correction.
   - fallback/default: max face-normal cell length.

5. Implement `deltaDES_isub5`:
   - Allocate and fill `VTM_data`.
   - Average `VTM_data` with boundary-aware stencils.
   - Compute geometry scale.
   - Apply FKH/VTM correction.
   - Deallocate `VTM_data`.

6. Add SA-IDDES branch:
   - Trigger on `ides = 4` or `izdes = 4`.
   - Compute `fdn`, `fe`, `delta`, `deltani`.
   - Store hybrid wall distance in `damp1`.
   - Ensure downstream SA source terms use `damp1`.

7. Add SST-IDDES branch:
   - Trigger on `ides = 4` or `izdes = 4`.
   - Trigger simplified IDDES on `ides = 5` or `izdes = 5`.
   - Compute `fdn`, `fe`, `delta`, `deltan`.
   - Compute `cdesn` using `icdes`.
   - Store hybrid length scale in `xlscale`.
   - Ensure the k-equation destruction term uses `xlscale`.

8. Verify with cases:
   - RANS baseline: `ides = 0`.
   - DES/DDES sanity: `ides = 1` and `ides = 2`.
   - IDDES: `ides = 4`.
   - Simplified SST IDDES: `ides = 5`.
   - Each `isub` option, especially `isub = 4` and `isub = 5`.
   - 3D separated-flow case where DES/IDDES should activate.

## 17. Common Pitfalls

- Do not compute `delta` before velocity gradients and cell metrics are available.
- Do not forget the `isub = 5` allocation and deallocation of `VTM_data`.
- Do not use `isub = 5` without computing `VTM_data` first.
- Do not apply IDDES in only one coordinate direction; it is a turbulence-model length-scale change, not a directional flux scheme.
- Do not confuse `ides = 4` IDDES with `ides = 6`, which in this project is an original DDES-SST / mode-III style branch.
- Do not enable both `ides` and `izdes`.
- Be careful with `icdes = 1`: it relies on the SST `blend` field being valid.

## 18. Summary

This project achieves IDDES by:

1. Reading DES-family keywords from `readkey.F`.
2. Selecting DES/IDDES branches in `spalart.F` or `twoeqn.F`.
3. Computing a grid/subgrid scale `delta` according to `isub`.
4. Computing shielding and elevating functions such as `fdn` and `fe`.
5. Building a hybrid RANS/LES length scale:
   - `damp1` for SA.
   - `xlscale` for SST/two-equation models.
6. Feeding that hybrid length scale into the turbulence-model source terms.

The `isub` keyword is the main control over how the LES grid scale `delta` is computed. For a faithful port of this project’s behavior, implement all `isub = 0..5` paths and preserve the `ides = 4` IDDES branches in both SA and SST turbulence routines.

