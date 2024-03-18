c********************************************************************ZLB
c  2-9-95  modified for SiH2+, two channels, pcs
c 11-5-94  modified for arbitrary number of channels, pcs
c          for high energies increase xmax
c----------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c---------------------------------------------------------
      subroutine asymp(rn,rki,eta,ei,jv,rjl,rjlp,rnl,rnlp,
     1                 unit,nch,ndat,npar,npot,xmax,nofg)
      implicit real*8 (a-h,o-z)
      integer i,nch,npar,l,ln,kfn,mode1,ifail
      double precision rho,deta,xmin,xmax
czlb      parameter(nch=3,npar=64,xmax=2000.0d0)
      double precision jv(1),unit,rn,a,b,t,rtki
c shared arrays
czlb
         dimension unit(nch,nch)
c         common/une/unit(nch,nch)
czlb
      double precision rki(nch),eta(nch),ei(nch),rjl(npar,nch,nch)
     1,rjlp(npar,nch,nch),rnl(npar,nch,nch),rnlp(npar,nch,nch)
c local arrays
       dimension l(nch)
czlb       double precision fc(3001,nch),gc(3001,nch),fcp(3001,nch),
czlb     .   gcp(3001,nch),f1(3001),f2(3001),g1(3001),g2(3001)
       double precision fc(nofg,nch),gc(nofg,nch),fcp(nofg,nch),
     .   gcp(nofg,nch),f1(nofg),f2(nofg),g1(nofg),g2(nofg)

c asymp return regular and irregular asymptotic function matrices.

         xmin=0.d0
         do 112 i=1,nch
         kfn=0
         rho=rn*rki(i)
         mode1=1

         if(eta(i).eq.0.0d0)  then
         kfn=1
         t=rho
         a=1.0d0
         b=-1.0d0
         else
         kfn=0
         t=1.0d0
         b=1.0d0
         a=0.d0
         end if

         deta=eta(i)

         call coulfg(rho,deta,xmin,xmax,f1,g1,f2,g2,mode1,kfn,ifail)

         maxl=dint(xmax)
      


         do 111 k=1,maxl-1
         fc(k,i)=f1(k)*t
         fcp(k,i)=a*f1(k)+t*f2(k)
         gc(k,i)=g1(k)*t*b
         gcp(k,i)=b*(a*g1(k)+t*g2(k))
111      continue
112      continue

30       continue
c construct zero-order asymptotic functions

         do 313 i=1,nch
313      l(i)=idint(jv(1))

         do 412 i=1,nch
         do 412 j=1,nch
         rtki=dsqrt(rki(j))
         ln=l(j)+1
         rjl(1,i,j)=unit(i,j)*fc(ln,j)/rtki
         rjlp(1,i,j)=unit(i,j)*fcp(ln,j)*rtki
         rnl(1,i,j)=unit(i,j)*gc(ln,j)/rtki
         rnlp(1,i,j)=unit(i,j)*gcp(ln,j)*rtki

412      continue

c         do k=1,npar
c         write(15,*)rjl(k,1,1),rjlp(k,1,1),rnl(k,1,1),rnlp(k,1,1)
c         end do

         return
         end
