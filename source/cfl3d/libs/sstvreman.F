c  ---------------------------------------------------------------------------
c  SST-Vreman hybrid eddy-viscosity model, version 1.
c
c  Standard SST remains responsible for the k-omega transport equations.
c  This file computes a separate viscosity for the mean-flow equations:
c
c       mu_t,H = (1-gamma)*mu_t,SST + gamma*mu_t,Vreman
c
c  The code variables are nondimensional dynamic viscosities.  Consequently
c  the Vreman kinematic viscosity is multiplied by rho*re before blending.
c  ---------------------------------------------------------------------------

      subroutine sstvreman_field(jdim,kdim,idim,q,x,y,z,sj,sk,si,vol,
     + ux,smin,zksav,fnu,vistsst,vishyb,blank)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension q(jdim,kdim,idim,5),x(jdim,kdim,idim),
     + y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     + si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension smin(jdim-1,kdim-1,idim-1)
      dimension zksav(jdim,kdim,idim,2)
      dimension fnu(0:jdim,0:kdim,0:idim)
      dimension vistsst(jdim,kdim,idim),vishyb(jdim,kdim,idim)
      dimension blank(jdim,kdim,idim)
      common /hybrid/ hyb_cv,ihyb_vreman,ihyb_shield,
     .                ihyb_resolution,ihyb_diagnostics,ihyb_test
      common /des/ cdes,ides,cddes,cdesw,cdese,zface,
     .             izdes,isub,icdes,iddes
      common /reyue/ reue,tinf,ivisc(3)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev

      re=reue/xmach
      do i=1,idim-1
        do k=1,kdim-1
          do j=1,jdim-1
            call sstvreman_cell(j,k,i,jdim,kdim,idim,q,x,y,z,sj,sk,
     +        si,vol,ux,smin,zksav,fnu(j,k,i),vistsst,re,xmutv,fd,fres,
     +        ghyb,xlrans,xlgrid,aalpha,bbeta,di,dj,dk)
            if (real(blank(j,k,i)) .le. 0.) then
              vishyb(j,k,i)=vistsst(j,k,i)
            else
              vishyb(j,k,i)=(1.-ghyb)*vistsst(j,k,i)+ghyb*xmutv
              vishyb(j,k,i)=ccmaxrc(0.,vishyb(j,k,i))
              vishyb(j,k,i)=ccmin(vishyb(j,k,i),edvislim)
            end if
          enddo
        enddo
      enddo

c     Fill the unused outer cell planes by constant extrapolation.  Exact
c     interface values are supplied through the normal CFL3D boundary arrays.
      do i=1,idim-1
        do k=1,kdim-1
          vishyb(jdim,k,i)=vishyb(jdim-1,k,i)
        enddo
        do j=1,jdim
          vishyb(j,kdim,i)=vishyb(min(j,jdim-1),kdim-1,i)
        enddo
      enddo
      do k=1,kdim
        do j=1,jdim
          vishyb(j,k,idim)=vishyb(min(j,jdim-1),min(k,kdim-1),idim-1)
        enddo
      enddo
      return
      end


      subroutine sstvreman_cell(j,k,i,jdim,kdim,idim,q,x,y,z,sj,sk,
     + si,vol,ux,smin,zksav,fnu,vistsst,re,xmutv,fd,fres,ghyb,
     + xlrans,xlgrid,aalpha,bbeta,deltai,deltaj,deltak)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension q(jdim,kdim,idim,5),x(jdim,kdim,idim),
     + y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     + si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension smin(jdim-1,kdim-1,idim-1)
      dimension zksav(jdim,kdim,idim,2)
      dimension vistsst(jdim,kdim,idim)
      dimension ex(3),ey(3),ez(3),gu(3),gv(3),gw(3)
      common /hybrid/ hyb_cv,ihyb_vreman,ihyb_shield,
     .                ihyb_resolution,ihyb_diagnostics,ihyb_test
      common /des/ cdes,ides,cddes,cdesw,cdese,zface,
     .             izdes,isub,icdes,iddes

      epsa=1.e-30
      deltaj=2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4)+epsa)
      deltak=2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4)+epsa)
      deltai=2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4)+epsa)

c     Grid-line directions connect the centers of opposite cell faces.
      ex(1)=.25*((x(j,k,i+1)+x(j+1,k,i+1)+x(j,k+1,i+1)+
     + x(j+1,k+1,i+1))-(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+
     + x(j+1,k+1,i)))
      ex(2)=.25*((y(j,k,i+1)+y(j+1,k,i+1)+y(j,k+1,i+1)+
     + y(j+1,k+1,i+1))-(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+
     + y(j+1,k+1,i)))
      ex(3)=.25*((z(j,k,i+1)+z(j+1,k,i+1)+z(j,k+1,i+1)+
     + z(j+1,k+1,i+1))-(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+
     + z(j+1,k+1,i)))
      ey(1)=.25*((x(j+1,k,i)+x(j+1,k+1,i)+x(j+1,k,i+1)+
     + x(j+1,k+1,i+1))-(x(j,k,i)+x(j,k+1,i)+x(j,k,i+1)+
     + x(j,k+1,i+1)))
      ey(2)=.25*((y(j+1,k,i)+y(j+1,k+1,i)+y(j+1,k,i+1)+
     + y(j+1,k+1,i+1))-(y(j,k,i)+y(j,k+1,i)+y(j,k,i+1)+
     + y(j,k+1,i+1)))
      ey(3)=.25*((z(j+1,k,i)+z(j+1,k+1,i)+z(j+1,k,i+1)+
     + z(j+1,k+1,i+1))-(z(j,k,i)+z(j,k+1,i)+z(j,k,i+1)+
     + z(j,k+1,i+1)))
      ez(1)=.25*((x(j,k+1,i)+x(j+1,k+1,i)+x(j,k+1,i+1)+
     + x(j+1,k+1,i+1))-(x(j,k,i)+x(j+1,k,i)+x(j,k,i+1)+
     + x(j+1,k,i+1)))
      ez(2)=.25*((y(j,k+1,i)+y(j+1,k+1,i)+y(j,k+1,i+1)+
     + y(j+1,k+1,i+1))-(y(j,k,i)+y(j+1,k,i)+y(j,k,i+1)+
     + y(j+1,k,i+1)))
      ez(3)=.25*((z(j,k+1,i)+z(j+1,k+1,i)+z(j,k+1,i+1)+
     + z(j+1,k+1,i+1))-(z(j,k,i)+z(j+1,k,i)+z(j,k,i+1)+
     + z(j+1,k,i+1)))
      call sstvreman_unit(ex)
      call sstvreman_unit(ey)
      call sstvreman_unit(ez)

c     Directional derivatives in the local filter basis.  ux ordering is
c     ux,uy,uz,vx,vy,vz,wx,wy,wz.
      gu(1)=ex(1)*ux(j,k,i,1)+ex(2)*ux(j,k,i,2)+
     +      ex(3)*ux(j,k,i,3)
      gu(2)=ey(1)*ux(j,k,i,1)+ey(2)*ux(j,k,i,2)+
     +      ey(3)*ux(j,k,i,3)
      gu(3)=ez(1)*ux(j,k,i,1)+ez(2)*ux(j,k,i,2)+
     +      ez(3)*ux(j,k,i,3)
      gv(1)=ex(1)*ux(j,k,i,4)+ex(2)*ux(j,k,i,5)+
     +      ex(3)*ux(j,k,i,6)
      gv(2)=ey(1)*ux(j,k,i,4)+ey(2)*ux(j,k,i,5)+
     +      ey(3)*ux(j,k,i,6)
      gv(3)=ez(1)*ux(j,k,i,4)+ez(2)*ux(j,k,i,5)+
     +      ez(3)*ux(j,k,i,6)
      gw(1)=ex(1)*ux(j,k,i,7)+ex(2)*ux(j,k,i,8)+
     +      ex(3)*ux(j,k,i,9)
      gw(2)=ey(1)*ux(j,k,i,7)+ey(2)*ux(j,k,i,8)+
     +      ey(3)*ux(j,k,i,9)
      gw(3)=ez(1)*ux(j,k,i,7)+ez(2)*ux(j,k,i,8)+
     +      ez(3)*ux(j,k,i,9)

      b11=deltai*deltai*gu(1)*gu(1)+
     +    deltaj*deltaj*gu(2)*gu(2)+deltak*deltak*gu(3)*gu(3)
      b22=deltai*deltai*gv(1)*gv(1)+
     +    deltaj*deltaj*gv(2)*gv(2)+deltak*deltak*gv(3)*gv(3)
      b33=deltai*deltai*gw(1)*gw(1)+
     +    deltaj*deltaj*gw(2)*gw(2)+deltak*deltak*gw(3)*gw(3)
      b12=deltai*deltai*gu(1)*gv(1)+
     +    deltaj*deltaj*gu(2)*gv(2)+deltak*deltak*gu(3)*gv(3)
      b13=deltai*deltai*gu(1)*gw(1)+
     +    deltaj*deltaj*gu(2)*gw(2)+deltak*deltak*gu(3)*gw(3)
      b23=deltai*deltai*gv(1)*gw(1)+
     +    deltaj*deltaj*gv(2)*gw(2)+deltak*deltak*gv(3)*gw(3)
      aalpha=0.
      do n=1,9
        aalpha=aalpha+ux(j,k,i,n)*ux(j,k,i,n)
      enddo
      bbeta=b11*b22-b12*b12+b11*b33-b13*b13+b22*b33-b23*b23
      bbeta=ccmaxrc(0.,bbeta)
      xmutv=0.
      if (real(aalpha).gt.1.e-30 .and. real(bbeta).gt.1.e-30) then
        xmutv=q(j,k,i,1)*re*hyb_cv*sqrt(bbeta/aalpha)
      endif

      dist=ccabs(smin(j,k,i))
      velterm=aalpha+epsa
      rd=(vistsst(j,k,i)+fnu)/
     +   (q(j,k,i,1)*sqrt(velterm)*.41*.41*dist*dist*re+epsa)
      fd=1.-cctanh((20.*rd)*(20.*rd)*(20.*rd))
      fd=ccmaxrc(0.,ccminrc(1.,fd))
      if (ihyb_shield .eq. 0) fd=1.

      xksafe=ccmax(zksav(j,k,i,2),1.e-30)
      osafe=ccmax(zksav(j,k,i,1),1.e-30)
      xlrans=sqrt(xksafe)/(.09*osafe*re+epsa)
      dh=ccmax(deltai,ccmax(deltaj,deltak))
      xlgrid=cdes*dh
      fres=ccmaxrc(0.,1.-xlgrid/(xlrans+epsa))
      fres=ccmaxrc(0.,ccminrc(1.,fres))
      if (ihyb_resolution .eq. 0) fres=1.
      ghyb=ccmaxrc(0.,ccminrc(1.,fd*fres))

c     User-declared laminar regions never activate the LES branch.
      if (real(smin(j,k,i)) .lt. 0.) ghyb=0.
      if (ihyb_test .eq. 1) ghyb=0.
      if (ihyb_test .eq. 2) ghyb=1.
      if (ihyb_test .eq. 3) xmutv=vistsst(j,k,i)
      return
      end


      subroutine sstvreman_bvis(jdim,kdim,idim,vistsst,vishyb,
     + vj0,vk0,vi0)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension vistsst(jdim,kdim,idim),vishyb(jdim,kdim,idim)
      dimension vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     + vi0(jdim,kdim,1,4)
      jlast=jdim-1
      klast=kdim-1
      ilast=idim-1

c     Preserve CFL3D's boundary treatment and extrapolate only the hybrid
c     correction from the adjacent interior cell into ghost planes 1 and 3.
      do i=1,ilast
        do k=1,klast
          vj0(k,i,1,1)=ccmaxrc(0.,vj0(k,i,1,1)+vishyb(1,k,i)-
     +      vistsst(1,k,i))
          vj0(k,i,1,3)=ccmaxrc(0.,vj0(k,i,1,3)+vishyb(jlast,k,i)-
     +      vistsst(jlast,k,i))
        enddo
        do j=1,jlast
          vk0(j,i,1,1)=ccmaxrc(0.,vk0(j,i,1,1)+vishyb(j,1,i)-
     +      vistsst(j,1,i))
          vk0(j,i,1,3)=ccmaxrc(0.,vk0(j,i,1,3)+vishyb(j,klast,i)-
     +      vistsst(j,klast,i))
        enddo
      enddo
      do k=1,klast
        do j=1,jlast
          vi0(j,k,1,1)=ccmaxrc(0.,vi0(j,k,1,1)+vishyb(j,k,1)-
     +      vistsst(j,k,1))
          vi0(j,k,1,3)=ccmaxrc(0.,vi0(j,k,1,3)+vishyb(j,k,ilast)-
     +      vistsst(j,k,ilast))
        enddo
      enddo
      return
      end


      subroutine sstvreman_coll(vistf,vistc,jdim,kdim,idim,
     + jj2,kk2,ii2)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension vistf(jdim,kdim,idim),vistc(jj2,kk2,ii2)
      kk=0
      do k=1,kdim-1,2
        kk=kk+1
        jj=0
        do j=1,jdim-1,2
          jj=jj+1
          ii=0
          do i=1,idim-1,2
            ii=ii+1
            vistc(jj,kk,ii)=.125*(vistf(j,k,i)+vistf(j,k,i+1)+
     +       vistf(j+1,k,i)+vistf(j+1,k,i+1)+vistf(j,k+1,i)+
     +       vistf(j,k+1,i+1)+vistf(j+1,k+1,i)+
     +       vistf(j+1,k+1,i+1))
          enddo
        enddo
      enddo
      do i=1,ii2-1
        do k=1,kk2-1
          vistc(jj2,k,i)=vistc(jj2-1,k,i)
        enddo
        do j=1,jj2
          vistc(j,kk2,i)=vistc(min(j,jj2-1),kk2-1,i)
        enddo
      enddo
      do k=1,kk2
        do j=1,jj2
          vistc(j,k,ii2)=vistc(min(j,jj2-1),min(k,kk2-1),ii2-1)
        enddo
      enddo
      return
      end


      subroutine sstvreman_unit(e)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension e(3)
      en=sqrt(e(1)*e(1)+e(2)*e(2)+e(3)*e(3)+1.e-30)
      e(1)=e(1)/en
      e(2)=e(2)/en
      e(3)=e(3)/en
      return
      end
