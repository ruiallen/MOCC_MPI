c********************************************************************ZLB
c apr 12 2022 Subrouine to output dcs
c Redesigned for MPI over partial wave
c--------------------------------------------
c num_p: number of processors, must be identical to Jstep
c--------------------------------------------
      subroutine diffcs(efreal,efreal2,efim,efim2,ofreal,ofreal2,
     & ofim,ofim2,rmat,imat,rki,p,p2,lstart,nch,npass,numa,num_p,ii)
      implicit real*8 (a-h,o-z)
      
      include 'mpif.h'
      integer a,i,n,nn,lstart
      integer npass,numa,num_p,nch,my_id,ii,index,l,ll
crw      parameter(nch=3,npar=64)
      double precision efreal,efreal2,efim,efim2
crw
      double precision ofreal,ofreal2,ofim,ofim2
      double precision rmat,imat,rki,p,p2
      dimension efreal(numa,nch,nch),efreal2(numa,nch,nch)
      dimension efim(numa,nch,nch),efim2(numa,nch,nch)
      dimension ofreal(numa,nch,nch),ofreal2(numa,nch,nch)
      dimension ofim(numa,nch,nch),ofim2(numa,nch,nch)
      dimension rmat(npass/num_p,nch,nch),imat(npass/num_p,nch,nch)
      dimension rki(nch),p(numa,0:npass),p2(numa,0:npass)
czlb
c
      
      do a = 1,numa
 
      do i = 0,npass/num_p-1,2
c elastic for n = 1
      l = i+lstart
c      write(*,*) l,rmat(i,1,1),p(a,l)
      efreal(a,1,1)=efreal(a,1,1)+
     &1./(2.*rki(1))*(2*l+1)*(rmat(i+1,1,1)-1)*p(a,l)
      efim(a,1,1)=efim(a,1,1)+
     &1./(2.*rki(1))*(2*l+1)*(imat(i+1,1,1))*p(a,l)
      efreal2(a,1,1)=efreal2(a,1,1)+
     &1./2./rki(1)*(2*l+1)*(rmat(i+1,1,1)-1)*p2(a,l)
      efim2(a,1,1)=efim2(a,1,1)+
     &1./2./rki(1)*(2*l+1)*imat(i+1,1,1)*p2(a,l)
c elastic for n >1
      do n = 2,nch
      efreal(a,n,n)=efreal(a,n,n)+
     &1./(2.*rki(n))*(2*l+1)*(rmat(i+1,n,n)-1)*p(a,l)
      efim(a,n,n)=efim(a,n,n)+
     & 1./(2.*rki(n))*(2*l+1)*imat(i+1,n,n)*p(a,l)
      efreal2(a,n,n)=efreal2(a,n,n)+
     &1./(2.*rki(n))*(2*l+1)*(rmat(i+1,n,n)-1)*p2(a,l)
      efim2(a,n,n)=efim2(a,n,n)+
     &1./(2.*rki(n))*(2*l+1)*imat(i+1,n,n)*p2(a,l)
      do nn = 1,n-1
c all inelastic case
      efreal(a,nn,n)=efreal(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(rmat(i+1,nn,n))*p(a,l)
      efim(a,nn,n)=efim(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(imat(i+1,nn,n))*p(a,l)
      efreal2(a,nn,n)=efreal2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(rmat(i+1,nn,n))*p2(a,l)
      efim2(a,nn,n)=efim2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(imat(i+1,nn,n))*p2(a,l)
      enddo
      enddo
      enddo
      do i = 1,npass/num_p-1,2
c elastic for n = 1
      l = i+lstart
c      write(*,*) l,rmat(i,1,1),p(a,l)
      ofreal(a,1,1)=ofreal(a,1,1)+
     &1./(2.*rki(1))*(2*l+1)*(rmat(i+1,1,1)-1)*p(a,l)
      ofim(a,1,1)=ofim(a,1,1)+
     &1./(2.*rki(1))*(2*l+1)*(imat(i+1,1,1))*p(a,l)
      ofreal2(a,1,1)=ofreal2(a,1,1)+
     &1./2./rki(1)*(2*l+1)*(rmat(i+1,1,1)-1)*p2(a,l)
      ofim2(a,1,1)=ofim2(a,1,1)+
     &1./2./rki(1)*(2*l+1)*imat(i+1,1,1)*p2(a,l)
c elastic for n >1
      do n = 2,nch
      ofreal(a,n,n)=ofreal(a,n,n)+
     &1./(2.*rki(n))*(2*l+1)*(rmat(i+1,n,n)-1)*p(a,l)
      ofim(a,n,n)=ofim(a,n,n)+
     & 1./(2.*rki(n))*(2*l+1)*(imat(i+1,n,n))*p(a,l)
      ofreal2(a,n,n)=ofreal2(a,n,n)+
     &1./2./rki(n)*(2*l+1)*(rmat(i+1,n,n)-1)*p2(a,l)
      ofim2(a,n,n)=ofim2(a,n,n)+
     &1./2./rki(n)*(2*l+1)*imat(i+1,n,n)*p2(a,l)
      do nn = 1,n-1
c all inelastic case
      ofreal(a,nn,n)=ofreal(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(rmat(i+1,nn,n))*p(a,l)
      ofim(a,nn,n)=ofim(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(imat(i+1,nn,n))*p(a,l)
      ofreal2(a,nn,n)=ofreal2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*rmat(i+1,nn,n)*p2(a,l)
      ofim2(a,nn,n)=ofim2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*l+1)*(imat(i+1,nn,n))*p2(a,l)
      enddo
      enddo
      enddo
      enddo

       return
       end
