c  DES grid-length definitions used by hybrid turbulence models.
      subroutine deltaDES(isub,ides,j,k,i,jdim,kdim,idim,x,y,z,
     + sj,sk,si,vol,ux,fd,delta)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     + sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     + si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),
     + ux(jdim-1,kdim-1,idim-1,9),xv(8),yv(8),zv(8)
      deltaj=2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
      deltak=2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
      deltai=2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
      delta=ccmax(deltaj,ccmax(deltak,deltai))
      if (isub .eq. 1) then
        delta=vol(j,k,i)**(1./3.)
        return
      endif
      wx=ux(j,k,i,8)-ux(j,k,i,6)
      wy=ux(j,k,i,3)-ux(j,k,i,7)
      wz=ux(j,k,i,4)-ux(j,k,i,2)
      wraw=sqrt(wx*wx+wy*wy+wz*wz)
      if (isub .eq. 2 .and. real(wraw) .lt. 1.e-15) return
      wmag=sqrt(wx*wx+wy*wy+wz*wz+1.e-30)
      wx=wx/wmag
      wy=wy/wmag
      wz=wz/wmag
      if (isub .eq. 2) then
        delta=sqrt(ccabs(wx)*deltaj*deltak+
     +             ccabs(wy)*deltak*deltai+
     +             ccabs(wz)*deltai*deltaj)
        return
      endif
      xv(1)=x(j,k,i)
      yv(1)=y(j,k,i)
      zv(1)=z(j,k,i)
      xv(2)=x(j+1,k,i)
      yv(2)=y(j+1,k,i)
      zv(2)=z(j+1,k,i)
      xv(3)=x(j,k+1,i)
      yv(3)=y(j,k+1,i)
      zv(3)=z(j,k+1,i)
      xv(4)=x(j+1,k+1,i)
      yv(4)=y(j+1,k+1,i)
      zv(4)=z(j+1,k+1,i)
      xv(5)=x(j,k,i+1)
      yv(5)=y(j,k,i+1)
      zv(5)=z(j,k,i+1)
      xv(6)=x(j+1,k,i+1)
      yv(6)=y(j+1,k,i+1)
      zv(6)=z(j+1,k,i+1)
      xv(7)=x(j,k+1,i+1)
      yv(7)=y(j,k+1,i+1)
      zv(7)=z(j,k+1,i+1)
      xv(8)=x(j+1,k+1,i+1)
      yv(8)=y(j+1,k+1,i+1)
      zv(8)=z(j+1,k+1,i+1)
      xc=0.
      yc=0.
      zc=0.
      do n=1,8
        xc=xc+xv(n)/8.
        yc=yc+yv(n)/8.
        zc=zc+zv(n)/8.
      enddo
      centermax=0.
      do n=1,8
        rx=xv(n)-xc
        ry=yv(n)-yc
        rz=zv(n)-zc
        rdotw=rx*wx+ry*wy+rz*wz
        rperp=sqrt(ccabs(rx*rx+ry*ry+rz*rz-rdotw*rdotw))
        centermax=ccmax(centermax,rperp)
      enddo
      delta=centermax/sqrt(3.)
      if (isub .eq. 4) then
        call des_vtm_cell(ux,j,k,i,jdim,kdim,idim,vtm)
        call des_fkh(vtm,ides,fd,fkhvtm)
        delta=delta*fkhvtm
      endif
      return
      end

      subroutine des_vtm_cell(ux,j,k,i,jdim,kdim,idim,vtm)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension ux(jdim-1,kdim-1,idim-1,9),ug(9)
      do n=1,9
        ug(n)=ux(j,k,i,n)
      enddo
      call des_vtm(ug,vtm)
      return
      end

      subroutine des_vtm(ug,vtm)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension ug(9)
      s11=ug(1)
      s22=ug(5)
      s33=ug(9)
      s12=.5*(ug(2)+ug(4))
      s13=.5*(ug(3)+ug(7))
      s23=.5*(ug(6)+ug(8))
      wx=ug(8)-ug(6)
      wy=ug(3)-ug(7)
      wz=ug(4)-ug(2)
      swx=s11*wx+s12*wy+s13*wz
      swy=s12*wx+s22*wy+s23*wz
      swz=s13*wx+s23*wy+s33*wz
      sn=sqrt(s11*s11+s22*s22+s33*s33+
     + 2.*(s12*s12+s13*s13+s23*s23)+1.e-30)
      wn=sqrt(wx*wx+wy*wy+wz*wz+1.e-30)
      vtm=sqrt(swx*swx+swy*swy+swz*swz)/(sn*wn+1.e-30)
      return
      end

      subroutine des_fkh(vtm,ides,fd,fkhvtm)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      fkhvtm=.1+.9*(vtm-.15)/.15
      fkhvtm=ccmaxrc(.1,ccminrc(1.,fkhvtm))
      if (ides .eq. 2 .and. real(fd) .lt. 1.-1.e-6) fkhvtm=1.
      if (ides .eq. 4 .and. real(fd) .gt. 1.e-6) fkhvtm=1.
      return
      end
