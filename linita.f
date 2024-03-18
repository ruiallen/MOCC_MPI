c********************************************************************ZLB
c 11-5-94 modified by pcs
c linita initializes channel angular momenta
c for "type-a" parity channels.
c shared arrays
c--------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c--------------------------------------------
      subroutine linita(jmin,jstep,jv,erel,cent,nch,ndat,npar,npot)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      integer nch,npar,jstep
czlb      parameter(nch=3,npar=64)
      double precision twomu,e,erel,cent,scale,step,
     . l(npar),jv(npar),jmin,jval
czlb
      dimension erel(nch),cent(npar,nch)
      common/pot/twomu,e,scale
c      common/pot/twomu,e,erel(nch),cent(npar,nch),scale
czlb
c
       jval=jmin

       step=dfloat(jstep)
       jv(1)=jval
       jval=jval+jstep
       do i20=1,nch
          l(1)=jv(1)
          cent(1,i20)=l(1)*(l(1)+1.0d0)
       enddo

       return
       end
