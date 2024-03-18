c********************************************************************ZLB
c 11-5-94 modified by pcs for n-channel problem
c potval returns the potentials u(npar,nch,nch)
c at radius r, for an nch channel problem.
c h is the step size.
c all other parameters are in common block POT.
c------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c------------------------------------------------
      subroutine potval(r,h,u,erel,cent,unit,lam,coupindx,
     1                  ra,pot,ypot,alpha,z,z2,
     2                  ea,delta,cinf,impot,nch,ndat,npar,npot,
     3                  xi,r2av,jstart)
      implicit real*8 (a-h,o-z)
      integer nch, npar, coupindx
czlb      parameter(nch=3,npar=64,ndat=31)
czlb        parameter(npot=nch*(nch+1)/2)
      double precision v(npar,nch,nch)
c shared arrays
      double precision u(npar,nch,nch),
     . twomu,e,erel,cent,scale,unit,r,lam,lamm,jstart
czlb
      dimension erel(nch),cent(npar,nch)
      dimension unit(nch,nch)
      dimension lam(nch),coupindx(nch,nch)
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension ra(ndat)
      dimension alpha(nch),z(nch),z2(nch)
      dimension ea(nch),cinf(npot),impot(nch)
c
      common/pot/twomu,e,scale
      common/fit/rs,rl,nprint
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb
c
c      common/pot/twomu,e,erel(nch),cent(npar,nch),scale
c      common/une/unit(nch,nch)
c      common/lamm/lam(nch),coupindx(nch,nch)
czlb
c potential is rescaled by alpha parameter.
       call dpot(r,v,ra,pot,ypot,alpha,z,z2,
     1           ea,delta,cinf,impot,nch,ndat,npar,npot,
     2           xi,r2av)
C  CALCULATE POTENTIAL AT R = T AND SCALE BY -H**2/6
c note: v already contains value of threshold energies
      do 10 k=1,npar
c rotational couplings
c       v(k,1,2)=(2.0d0/twomu/r**2*(cent(k,1))**0.5)*v(k,1,2)
c       v(k,2,1)=v(k,1,2)
c       v(k,2,3)=(2.0d0/twomu/r**2*(cent(k,1))**0.5)*v(k,2,3)
c       v(k,3,2)=v(k,2,3)
       do 10 i=1,nch
       do 10 j=1,nch
          if(coupindx(i,j).eq.2) then
             lamm=min(lam(i),lam(j))
             if((dfloat(k-1)+jstart).ge.lamm) then
                v(k,i,j)=2.0d0/twomu/r**2*v(k,i,j)
     .           *((dfloat(k-1)+jstart-lamm)
     .           *(dfloat(k-1)+jstart+lamm+1.0d0))**0.5
c           if(r.lt.0.055) print*,dfloat(k-1)+jstart,i,j,lamm,
c     .            ((dfloat(k-1)+jstart-lamm)
c     .          *(dfloat(k-1)+jstart+lamm+1.0d0))**0.5,
c     .           v(k,i,j)
             else
               v(k,i,j)=0.0d0
             end if
          end if
          if(coupindx(i,j).gt.2)  then
             v(k,i,j)=0.0d0
c          print*,k,i,j
          end if
cnew 12/22/2004
c pcs, to ensure J(J+1)-lambda^2 >0, i.e. no contributions from these
c      partial waves due to radial coupling
          if((coupindx(i,j).eq.1).and.(cent(k,i).lt.lam(i)**2))
     .         v(k,i,j)=0.0d0
c
c all potentials go to zero NOT THRESHOLD Energy. Put that into erel(i)
      u(k,i,j)=-scale*(twomu*(v(k,i,j)-erel(i)*unit(i,j))+
     1unit(i,j)*(cent(k,i)-lam(i)**2)/(r*r))
      u(k,j,i)=u(k,i,j)
cnew 12/22/2004
10    continue
c      write(23,*)r,v(1,5,5)+420.0d0/(r*r)/twomu,
c     .   v(1,5,5)+650.0d0/(r*r)/twomu
c      write(20,*)r,v(1,1,1)-erel(1)+e,v(1,2,2)-erel(2)+e
c      write(21,*)r,v(1,3,3)-erel(3)+e,v(1,4,4)-erel(4)+e
c      write(22,*)r,v(1,3,1),v(1,3,2)
c      write(26,*)r,v(1,3,4),v(1,3,5)
c      write(27,*)r,v(1,4,5)
      return
      end
c********************************************************************ZLB
