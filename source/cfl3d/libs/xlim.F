c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c
c  The CFL3D platform is licensed under the Apache License, Version 2.0
c  (the "License"); you may not use this file except in compliance with the
c  License. You may obtain a copy of the License at
c  http://www.apache.org/licenses/LICENSE-2.0.
c
c  Unless required by applicable law or agreed to in writing, software
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
c  License for the specific language governing permissions and limitations
c  under the License.
c  ---------------------------------------------------------------------------
c
      subroutine xlim(xkap,n,x1,x2,xc,iswi,npts,leq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Perform monotone interpolations to the interfaces
c     of the cells.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(n),x2(n),xc(n)
c
      common /chk/ ichk
c
c     on input:
c        x1 is left interface gradient del(-)
c        x2 is right interface gradient del(+) of cell
c        xc is cell center value
c     on output:
c        x2 is gradient for right interface construction
c        x1 is gradient for left interface construction
c     for unlimited fully-upwind (kappa=-1;iswi=0)
c        x1(output)=x2(input)
c        x2(output)=x1(input)
c
      ibias = 1
      eps   = 1.0e-06
      phi   = (3.-xkap)/(1.-xkap+eps)
      eps2  = 0.5e-06
      zeroc = 0.
c
      if (real(xkap).lt.-1.0) then
c
c      first order
c
cdir$ ivdep
      do 8200 izz=1,n
      x2(izz) = 0.
      x1(izz) = 0.
 8200 continue
      return
      end if
c
c      bias gradients by average values (leq = 1 or 5)
c
      if (ibias.eq.1) then
      if (leq.eq.1 .or. leq.eq.5) then
cdir$ ivdep
      do 8000 izz=1,n
      x2avg   = xc(izz) + x2(izz)*0.5
      x1avg   = xc(izz) - x1(izz)*0.5
      x2(izz) = x2(izz)/x2avg
      x1(izz) = x1(izz)/x1avg
 8000 continue
      end if
      end if
c
      if (iswi.eq.-1) then
c
c     od - second order
c
cdir$ ivdep
      do 1000 izz=1,n
      x3      = ccabs(x1(izz))
      x4      = ccabs(x2(izz))
      x2(izz) = ccvmgt(x1(izz)*0.5e0,x2(izz),(real(x3).lt.real(x4)))
      x1(izz) = ccvmgt(x2(izz),x1(izz),(real(x3).lt.real(x4)))
      x1(izz) = x2(izz)
 1000 continue
      else if (iswi.eq.1) then
c
c      smooth (van albeda) limiter
c
      term    = xkap*4.e0
cdir$ ivdep
      do 1001 izz=1,n
      x3      = x1(izz)*x1(izz)+x2(izz)*x2(izz)+eps
      x4      = (eps2+x1(izz)*x2(izz))/x3
      x4      = x4*0.5e0
      x3      = (x1(izz)+x2(izz))*x4
      x4      = (x2(izz)-x1(izz))*x4*term*x4
      x2(izz) = x3+x4
      x1(izz) = x3-x4
 1001 continue
      else if (iswi.eq.0) then
c
c      unlimited  (kappa scheme)
c
      term    = xkap*0.25e0
cdir$ ivdep
      do 1002 izz=1,n
      x3      = .25e0*(x1(izz)+x2(izz))
      x4      = term*(x2(izz)-x1(izz))
      x2(izz) = x3+x4
      x1(izz) = x3-x4
 1002 continue
      else if (iswi.eq.2) then
c
c      limited  (min-mod scheme)
c
      term    = xkap*0.25e0
cdir$ ivdep
      do 1003 izz=1,n
      x4      = x1(izz)*x2(izz)
      x1(izz) = ccvmgt(zeroc,x1(izz),(real(x4).lt.0.e0))
      x2(izz) = ccvmgt(zeroc,x2(izz),(real(x4).lt.0.e0))
      x4      = phi*x1(izz)
      x4      = ccabs(x4)
      x3      = phi*x2(izz)
      x3      = ccabs(x3)
      x5      = ccabs(x2(izz))
      x2(izz) = ccvmgt(x1(izz)*phi,x2(izz),(real(x4).lt.real(x5)))
      x5      = ccabs(x1(izz))
      x1(izz) = ccvmgt(x2(izz)*phi,x1(izz),(real(x3).lt.real(x5)))
      x3      = .25e0*(x1(izz)+x2(izz))
      x4      = term*(x2(izz)-x1(izz))
      x2(izz) = x3+x4
      x1(izz) = x3-x4
 1003 continue
      else if (iswi.eq.3 .or. iswi.eq.4) then
c
c     tuned k=1/3 limiter - Spekreijse - Venkat                   AIAA-90-0429
c
      delx    = 10./float(npts)
      eps2    = delx**3
cdir$ ivdep
      do 7500 izz=1,n
      t3      = x1(izz)*x1(izz)
      t4      = x2(izz)*x2(izz)
      t5      = x1(izz)*x2(izz)
      t6      = x1(izz)+x2(izz)
      term    = 0.5*(t5+eps2)/(2.*(t3 + t4) - t5 + 3.*eps2)
      x2(izz) = (x2(izz)+t6)*term
      x1(izz) = (x1(izz)+t6)*term
 7500 continue
      end if
c
c      bias gradients by average values (leq = 1 or 5)
c
      if (ibias.eq.1) then
      if (leq.eq.1 .or. leq.eq.5) then
cdir$ ivdep
      do 9000 izz=1,n
      x2(izz) = x2(izz)*xc(izz)
      x1(izz) = x1(izz)*xc(izz)
 9000 continue
      end if
      end if
c
c     cap density and pressure
c     - ensures that they stay positive
c     - limits their maxima to their stagnation values
c
      if (ichk.eq.2) then
         call prolim(n,x1,x2,xc,leq)
      end if
      if (ichk.eq.3) then
         call prolim2(n,x1,x2,xc,leq)
      end if
c
      return
      end
c
      subroutine xlim2(n,x1,x2,x3,x4)
c
c***********************************************************************
c     Purpose:  Perform fifth-order WENO-Z interface reconstruction.
c
c     On input:
c        x1 is del(-2), x2 is del(-1), x3 is del(+1), x4 is del(+2).
c     On output:
c        x2 is the correction for right-interface construction,
c        x1 is the correction for left-interface construction.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(n),x2(n),x3(n),x4(n)
c
      eps = 1.0e-02
c
cdir$ ivdep
      do 1000 izz=1,n
      t1 = x1(izz)
      t2 = x2(izz)
      t3 = x3(izz)
      t4 = x4(izz)
c
      ql0 =  5.0e0/6.0e0*t2 - 1.0e0/3.0e0*t1
      ql1 =  1.0e0/3.0e0*t3 + 1.0e0/6.0e0*t2
      ql2 = -1.0e0/6.0e0*t4 + 2.0e0/3.0e0*t3
c
      qr0 =  2.0e0/3.0e0*t2 - 1.0e0/6.0e0*t1
      qr1 =  1.0e0/6.0e0*t3 + 1.0e0/3.0e0*t2
      qr2 = -1.0e0/3.0e0*t4 + 5.0e0/6.0e0*t3
c
      beta0 = 13.0e0/12.0e0*(t2-t1)**2
     .      + 0.25e0*(3.0e0*t2-t1)**2
      beta1 = 13.0e0/12.0e0*(t3-t2)**2
     .      + 0.25e0*(t3+t2)**2
      beta2 = 13.0e0/12.0e0*(t4-t3)**2
     .      + 0.25e0*(t4-3.0e0*t3)**2
c
      tau = ccabs(beta0-beta2)
      b0  = 1.0e0 + tau/(eps+beta0)
      b1  = 1.0e0 + tau/(eps+beta1)
      b2  = 1.0e0 + tau/(eps+beta2)
c
      a0  = 0.1e0*b0
      a1  = 0.6e0*b1
      a2  = 0.3e0*b2
      asum = a0+a1+a2
      wl0 = a0/asum
      wl1 = a1/asum
      wl2 = a2/asum
c
      a0  = 0.3e0*b0
      a1  = 0.6e0*b1
      a2  = 0.1e0*b2
      asum = a0+a1+a2
      wr0 = a0/asum
      wr1 = a1/asum
      wr2 = a2/asum
c
      x2(izz) = wl0*ql0 + wl1*ql1 + wl2*ql2
      x1(izz) = wr0*qr0 + wr1*qr1 + wr2*qr2
 1000 continue
c
      return
      end
