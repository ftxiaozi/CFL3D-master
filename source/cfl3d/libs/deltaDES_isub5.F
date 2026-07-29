c  Extended VTM/FKH scale using a boundary-aware 3x3x3 VTM average.
      subroutine deltaDES_isub5(ides,j,k,i,jdim,kdim,idim,x,y,z,
     + sj,sk,si,vol,ux,fd,vtm_data,delta)
#ifdef CMPLX
      implicit complex(a-h,o-z)
#endif
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),
     + sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),
     + si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),
     + ux(jdim-1,kdim-1,idim-1,9),
     + vtm_data(jdim-1,kdim-1,idim-1)
      vtmave=0.
      nvtm=0
      do ii=max(1,i-1),min(idim-1,i+1)
        do kk=max(1,k-1),min(kdim-1,k+1)
          do jj=max(1,j-1),min(jdim-1,j+1)
            vtmave=vtmave+vtm_data(jj,kk,ii)
            nvtm=nvtm+1
          enddo
        enddo
      enddo
      vtmave=vtmave/float(nvtm)
      call deltaDES(3,ides,j,k,i,jdim,kdim,idim,x,y,z,sj,sk,si,
     + vol,ux,fd,delta)
      call des_fkh(vtmave,ides,fd,fkhvtm)
      delta=delta*fkhvtm
      return
      end
