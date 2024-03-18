c********************************************************************ZLB
c apr 12 2022 Subrouine to output dcs
c Redesigned for MPI over partial wave
c--------------------------------------------
c num_p: number of processors, must be identical to Jstep
c--------------------------------------------
      subroutine odiffcs(ofreal,ofreal2,ofim,ofim2,rmat,imat,
     & rki,p,p2,nch,npass,numa,num_p,my_id)
      implicit none
      
      include 'mpif.h'
      integer a,i,n,nn
      integer npass,numa,num_p,nch,my_id
crw      parameter(nch=3,npar=64)
      double precision ofreal,ofreal2,ofim,ofim2
crw
      double precision rmat,imat,rki,p,p2
      dimension ofreal(numa,nch,nch),ofreal2(numa,nch,nch)
      dimension ofim(numa,nch,nch),ofim2(numa,nch,nch)
      dimension rmat(0:npass/num_p,nch,nch),imat(0:npass/num_p,nch,nch)
      dimension rki(nch),p(numa,0:npass),p2(numa,0:npass)
czlb
c
      
      do a = 1,2
      do i = 0,npass/num_p-1,1
c elastic for n = 1
      write(*,*) "odd",my_id+i*num_p,rmat(i,1,1),p(a,my_id+i*num_p)
      ofreal(a,1,1)=ofreal(a,1,1)+
     &1./(2.*rki(1))*(2*i+1)*(rmat(i,1,1)-1)*p(a,my_id+i*num_p)
      ofim(a,1,1)=ofim(a,1,1)+
     &1./(2.*rki(1))*(2*i+1)*(imat(i,1,1))*p(a,my_id+i*num_p)
      ofreal2(a,1,1)=ofreal2(a,1,1)+
     &1./2./rki(1)*(2*i+1)*(rmat(i,1,1)-1)*p2(a,my_id+i*num_p)
      ofim2(a,1,1)=ofim2(a,1,1)+
     &1./2./rki(1)*(2*i+1)*imat(i,1,1)*p2(a,my_id+i*num_p)
c elastic for n >1
      do n = 2,nch
      ofreal(a,n,n)=ofreal(a,n,n)+
     &1./(2.*rki(n))*(2*i+1)*(rmat(i,n,n)-1)*p(a,my_id+i*num_p)
      ofim(a,n,n)=ofim(a,n,n)+
     & 1./(2.*rki(n))*(2*i+1)*(imat(i,n,n))*p(a,my_id+i*num_p)
      ofreal2(a,n,n)=ofreal2(a,n,n)+
     &1./2./rki(n)*(2*i+1)*(rmat(i,n,n)-1)*p2(a,my_id+i*num_p)
      ofim2(a,n,n)=ofim2(a,n,n)+
     &1./2./rki(n)*(2*i+1)*imat(i,n,n)*p2(a,my_id+i*num_p)
      do nn = 1,n-1
c all inelastic case
      ofreal(a,nn,n)=ofreal(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(rmat(i,nn,n))*p(a,my_id+i*num_p)
      ofim(a,nn,n)=ofim(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(imat(i,nn,n))*p(a,my_id+i*num_p)
      ofreal2(a,nn,n)=ofreal2(a,nn,n)+1./(2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*rmat(i,nn,n)*p2(a,my_id+i*num_p)
      ofim2(a,nn,n)=ofim2(a,nn,n)+1./(-2.*SQRT(rki(n)*rki(nn)))
     & *(2*i+1)*(imat(i,nn,n))*p2(a,my_id+i*num_p)
      enddo
      enddo
      enddo
      enddo

       return
       end
