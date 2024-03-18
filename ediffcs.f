c********************************************************************ZLB
c apr 12 2022 Subrouine to output dcs
c Redesigned for MPI over partial wave
c--------------------------------------------
c num_p: number of processors, must be identical to Jstep
c--------------------------------------------
      subroutine ediffcs(efreal,efreal2,efim,efim2,rmat,imat,
     & rki,p,p2,nch,npass,numa,num_p,my_id)
      implicit real*8 (a-h,o-z)
      
      include 'mpif.h'
      integer a,i,n,nn
      integer npass,numa,num_p,nch,my_id
crw      parameter(nch=3,npar=64)
      double precision efreal,efreal2,efim,efim2
crw
      double precision rmat,imat,rki,p,p2
      dimension efreal(numa,nch,nch),efreal2(numa,nch,nch)
      dimension efim(numa,nch,nch),efim2(numa,nch,nch)
      dimension rmat(0:npass/num_p,nch,nch),imat(0:npass/num_p,nch,nch)
      dimension rki(nch),p(numa,0:npass),p2(numa,0:npass)
czlb
c
      
      do a = 1,2
 
      do i = 0,npass/num_p-1,1
c elastic for n = 1
      write(*,*) 'even',my_id+i*num_p,rmat(i,1,1),p(a,i)
      efreal(a,1,1)=efreal(a,1,1)+
     &1./(2.*rki(1))*(2*i+1)*(rmat(i,1,1)-1)*p(a,i)
      efim(a,1,1)=efim(a,1,1)+
     &1./(2.*rki(1))*(2*i+1)*(imat(i,1,1))*p(a,i)
      efreal2(a,1,1)=efreal2(a,1,1)+
     &1./2./rki(1)*(2*i+1)*(rmat(i,1,1)-1)*p2(a,i)
      efim2(a,1,1)=efim2(a,1,1)+
     &1./2./rki(1)*(2*i+1)*imat(i,1,1)*p2(a,i)
c elastic for n >1
      do n = 2,nch
      efreal(a,n,n)=efreal(a,n,n)+
     &1./(2.*rki(n))*(2*i+1)*(rmat(i,n,n)-1)*p(a,i)
      efim(a,n,n)=efim(a,n,n)+
     & 1./(2.*rki(n))*(2*i+1)*imat(i,n,n)*p(a,i)
      efreal2(a,n,n)=efreal2(a,n,n)+
     &1./(2.*rki(n))*(2*i+1)*(rmat(i,n,n)-1)*p2(a,i)
      efim2(a,n,n)=efim2(a,n,n)+
     &1./(2.*rki(n))*(2*i+1)*imat(i,n,n)*p2(a,i)
      do nn = 1,n-1
c all inelastic case
      efreal(a,nn,n)=efreal(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(rmat(i,nn,n))*p(a,i)
      efim(a,nn,n)=efim(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(imat(i,nn,n))*p(a,i)
      efreal2(a,nn,n)=efreal2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(rmat(i,nn,n))*p2(a,i)
      efim2(a,nn,n)=efim2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(imat(i,nn,n))*p2(a,i)
      enddo
      enddo
      enddo
      enddo

       return
       end
