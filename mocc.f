      program mocc
      use mpi
      integer nch, ndat,npar
     . ,npass,jstep, nprint,ne
     . ,i,j,k,ii,kcount,npot,ici,
     . i99, j99,jj,kk
      real*8 rs,rl
      double precision astep
      double precision jmin,xmax
      double precision x1,x2,diff_param
      common nch, npar, ndat, npot
      common/pot/twomu,e,scale
      common/count/lcount
      common/fit/rs,rl,nprint
czlb        06/27/04
      integer my_id,num_p,num_e_point,ierr
      integer an_id,lstart,lend
      integer STATUS(MPI_STATUS_SIZE)
      real*8 delta,twomu,scale,e,
     . emu,rmu,h,rstart,rn,estart,jin1,jin

c     file handling
      integer jnumber
      CHARACTER(len=255) :: cwd,outpath



      
      real*8,allocatable,dimension(:,:):: pot,ypot
      real*8,allocatable,dimension(:,:,:):: yn, rjl,rnl,
     & rjlp, rnlp
      real*8,allocatable,dimension(:):: rki,eta,gamma
      real*8,allocatable,dimension(:):: ei,ke,ev,erel,erelk
      real*8,allocatable,dimension(:):: sigt
      real*8,allocatable,dimension(:,:):: sig
      real*8,allocatable,dimension(:):: g,jv
      real*8,allocatable,dimension(:,:):: cent,unit
      real*8,allocatable,dimension(:,:,:):: fkmat,smat,parsig
      real*8,allocatable,dimension(:):: ea,cinf
      real*8,allocatable,dimension(:):: alpha,z,z2
      real*8,allocatable,dimension(:):: lam
      integer,allocatable,dimension(:,:):: coupindx
      integer,allocatable,dimension(:):: impot
      real*8,allocatable,dimension(:):: r,xi,r2av
      real*8,allocatable,dimension(:):: l
      real*8,allocatable,dimension(:):: parwave


# these params are for differential
      real*8,allocatable,dimension(:,:):: p,p2
      real*8,allocatable,dimension(:):: cosines,cosines2
      real*8,allocatable,dimension(:):: sines
      real*8,allocatable,dimension(:,:,:)::rmat,imat

      real*8,allocatable,dimension(:,:,:)::freal,fim
      real*8,allocatable,dimension(:,:,:)::freal2,fim2

      real*8,allocatable,dimension(:,:,:)::ofreal,ofim
      real*8,allocatable,dimension(:,:,:)::ofreal2,ofim2
      real*8,allocatable,dimension(:,:,:)::efreal,efim
      real*8,allocatable,dimension(:,:,:)::efreal2,efim2

      real*8,allocatable,dimension(:,:,:)::fdreal,fdim
      real*8,allocatable,dimension(:,:,:)::fexreal,fexim
      real*8,allocatable,dimension(:,:,:)::fdreal2,fdim2
      real*8,allocatable,dimension(:,:,:)::fexreal2,fexim2


      real*8,allocatable,dimension(:):: x,w,angles

      real*8,allocatable,dimension(:,:,:):: indist_diff
      real*8,allocatable,dimension(:,:,:):: dist_diff
      real*8,allocatable,dimension(:,:,:):: cxdiff
c these are for test only
      real*8,allocatable,dimension(:,:,:):: cxdiff2,dist_diff2

c single diff is inelastic only
      real*8,allocatable,dimension(:,:,:):: single_diff

      real*8,allocatable,dimension(:,:,:):: smr, smi




      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,num_p,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
      call get_environment_variable('PWD',cwd)
      open(unit = 777,file = 'ang920.dat')
      open(unit = 33,file = 'energy')
      open(unit=3, file='input.capture')
      open(unit=16, file='pot.diabatic',status='old')
      open(unit=7,file='diffgs.txt')
      open(unit=900,file='erealt.txt')
      open(unit=901,file='erealt2.txt')
      open(unit=902,file='eimt.txt')
      open(unit=903,file='eimt2.txt')
      open(unit=904,file='orealt.txt')
      open(unit=905,file='orealt2.txt')
      open(unit=906,file='oimt.txt')
      open(unit=907,file='oimt2.txt')

      open(unit=908,file='rmat')
      open(unit=909,file='imat')

      read(3,*) nch,ndat,npar,numa
      numa = 919
      npot=nch*(nch+1)/2

      allocate (cinf(npot))




      do i=1,npot
            cinf(i)=0.0d0
      end do

      allocate (yn(npar,nch,nch))
      allocate (rjl(npar,nch,nch))
      allocate (rnl(npar,nch,nch))
      allocate (rjlp(npar,nch,nch))
      allocate (rnlp(npar,nch,nch))
      allocate (rki(nch))
      allocate (eta(nch))
      allocate (gamma(nch))

      allocate (ke(nch))
      allocate (ei(nch))
      allocate (sig(nch,nch))
      allocate (sigt(nch))

      allocate (erelk(nch))

      allocate (g(nch))

      allocate (jv(npar))

      allocate (erel(nch))
      allocate (cent(npar,nch))

      allocate (unit(nch,nch))

      allocate (fkmat(npar,nch,nch))
      allocate (smat(npar,nch,nch))
      allocate (parsig(npar,nch,nch))

      allocate (ea(nch))

      allocate (impot(nch))

      allocate (alpha(nch))
      allocate (z(nch))
      allocate (z2(nch))

      allocate (lam(nch))
      allocate (coupindx(nch,nch))

      allocate (pot(ndat,npot))
      allocate (ypot(ndat,npot))

      allocate (r(ndat))

      allocate (xi(nch))
      allocate (r2av(nch))

      read(3,*) rstart,jmin,npass,jstep,nprint
      read(3,*) ne,h,rn,emu
      read(3,*) estart, he
      read(3,*) delta, rs, rl
      
      read(3,*) (ea(i),alpha(i),xi(i),r2av(i),z(i),z2(i),
     &            lam(i),g(i),impot(i),i=1,nch)

      read(3,*) (ei(i),i=1,nch)
      read(3,*) ici

      read(3,*) (j,cinf(j),i=1,ici)

      
      allocate(ev(ne))

      jstep = num_p
      nofg=npar*npass+101
      xmax=nofg-1

      allocate (cosines(numa))
      allocate (cosines2(numa))
      allocate (sines(numa))

      allocate (p(numa,0:npass))
      allocate (p2(numa,0:npass))

      allocate (angles(numa))

      
c      x1=-1.d0
c      x2= 1.d0
c      call gauleg(x1,x2,x,w,numa)
      read(777,*) (angles(i),i=1,numa)
      pi=4.0d0*datan(1.0d0)
      do j=1,numa
         cosines(j)=cos(angles(j))
         cosines2(j)=-1*cosines(j)
         sines(j)=sin(angles(j))
      end do


      do i=1,nch
        gamma(i)=z(i)*z2(i)
        do j=1,nch
          coupindx(i,j)=0
        enddo
      enddo

      do i=1,nch-1
        do j = i+1,nch
            coupindx(i,j)=dabs(dint(lam(i)-lam(j)))+1
            coupindx(j,i)=coupindx(i,j)
        enddo
      enddo



883      format('# state-to-state electron capture cross sections')
879      format('#    E              sigma_ij (10^-16 cm2)')
884      format('#  (eV/u)',8x,30(i2,i2,10x))
878      format('# Partial cross sections (a_0^2/g)')
877      format('#    J    sigma_ij(J)')
c *******

876      format('#',8x,30(i2,i2,10x))
      

      call input(r,pot,ypot,alpha,z,z2,nch,ndat,npar,npot,xi,r2av)

      call legendre(numa,npass,cosines,p)
      call legendre(numa,npass,cosines2,p2)

      do 556 kount=1,ne

      allocate (smr(npar,nch,nch))
      allocate (smi(npar,nch,nch))

      allocate (ofreal(numa,nch,nch))
      allocate (ofim(numa,nch,nch))
      allocate (ofreal2(numa,nch,nch))
      allocate (ofim2(numa,nch,nch))
      allocate (efreal(numa,nch,nch))
      allocate (efreal2(numa,nch,nch))
      allocate (efim(numa,nch,nch))
      allocate (efim2(numa,nch,nch))
      allocate (rmat(npass/num_p,nch,nch))
      allocate (imat(npass/num_p,nch,nch))

        write(*,*) 'E point =', kount

      ev(kount)=estart+dfloat(kount-1)*he

        e=ev(kount)/27.2113962d0

        rmu=1822.88732D0*emu

        pi=4.0d0*datan(1.0d0)
        DO 130 I=1,NCH
          erel(i) = e-ei(i)
          ke(i)=erel(i)*27.2113962d0/emu

          erelk(i)=erel(i)

          if(erel(i).lt.0.d0) then
            erelk(i)=dabs(erel(i))
          end if

          RKI(I)=DSQRT(2.D0*RMU*(ERELk(I)))
130       ETA(I)=RMU*GAMMA(I)/RKI(I)
          TWOMU=RMU+RMU
          scale=h*h/6.0d0
          do 201 i=1,nch
            sigt(i)=0.0d0
          do 201 j=1,nch
            sig(i,j)=0.0d0
          unit(i,j)=0.0d0
201       unit(j,j)=1.0d0

          jin = jmin+my_id*npass/num_p
          lstart = jin
            
           do 555 i = 1,npass/num_p
c      do 555 i = 1,npass,2
           
           call linita(jin,jstep,jv,erel,cent,nch,ndat,npar,npot)
           jin1=jin
           jin=jv(1)+1

           call asymp(rn,rki,eta,ei,jv,rjl,rjlp,rnl,rnlp,
     &           unit,nch,ndat,npar,npot,xmax,nofg)
            
           call ylog(my_id,h,rstart,rn,yn,erel,cent,unit,lam,coupindx,
     &          r,pot,ypot,alpha,z,z2,
     &          ea,delta,cinf,impot,nch,ndat,npar,npot,
     &          xi,r2av,jin1)
            


            call kmatrix(yn,rjl,rjlp,rnl,rnlp,jv,rki,nprint,
     &             unit,fkmat,smat,parsig,nch,ndat,npar,npot,
     & smr,smi)


      do i99=1,nch
      do j99=1,nch
       rmat(i,i99,j99) = smr(1,i99,j99)
       imat(i,i99,j99) = smi(1,i99,j99)
      
      enddo
      enddo
      
c      write(*,*) jv(1),rmat(i-1,1,1)
      

        do i99=1,npar
        do j99=1,nch
        do k99=1,nch
          if(parsig(i99,j99,k99).lt.1.0d-99) then
            parsig(i99,j99,k99)=0.0d0
          end if
        enddo
        enddo
        enddo

        

555   continue
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      
      if (my_id .ne. 0) then
      call MPI_SEND(lstart,1,MPI_INTEGER,0,11,MPI_COMM_WORLD,ierr)
      call MPI_SEND(rmat,npass/num_p*nch*nch,MPI_DOUBLE_PRECISION,0,1,
     & MPI_COMM_WORLD,ierr)
      call MPI_SEND(imat,npass/num_p*nch*nch,MPI_DOUBLE_PRECISION,0,2,
     & MPI_COMM_WORLD,ierr)

      else
c finish calculation with data on process 0
c      write(*,*) lstart
      do i = 1,npass/num_p
      write(908,*) i+lstart
      write(909,*) i+lstart
      do jj = 1,nch
      write(908,455) (rmat(i,jj,j),j=1,nch)
      write(909,455) (imat(i,jj,j),j=1,nch)
      enddo
      enddo

      call diffcs(efreal,efreal2,efim,efim2,ofreal,ofreal2,
     & ofim,ofim2,rmat,imat,rki,p,p2,lstart,nch,npass,numa,num_p,0)
      do ii=1,num_p-1
            
c finish calculation with data on other process, one by one
      call MPI_RECV(lstart,1,MPI_INTEGER,ii,11,
     & MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(rmat,npass/num_p*nch*nch,MPI_DOUBLE_PRECISION,ii,1,
     & MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(imat,npass/num_p*nch*nch,MPI_DOUBLE_PRECISION,ii,
     & 2,MPI_COMM_WORLD,status,ierr)

      do i = 1,npass/num_p
      write(908,*) i+lstart
      write(909,*) i+lstart
      do jj = 1,nch
        write(908,455) (rmat(i,jj,j),j=1,nch)
        write(909,455) (imat(i,jj,j),j=1,nch)
      enddo
      enddo
c      write(*,*) lstart
c now index should start from npass/num_p, how to modify index?????????/
      call diffcs(efreal,efreal2,efim,efim2,ofreal,ofreal2,
     & ofim,ofim2,rmat,imat,rki,p,p2,lstart,nch,npass,numa,num_p,ii)
      enddo
      endif

91    format('# Energy at ',1pe14.7,' eV/u')
92    format(10(1pe14.7))
93    format('angles at ',1pe14.7)





      if (my_id .eq. 0) then
      

      write(900,*) 'Energy at ', ke(1)
      write(901,*) 'Energy at ', ke(1)
      write(902,*) 'Energy at ', ke(1)
      write(903,*) 'Energy at ', ke(1)
      write(904,*) 'Energy at ', ke(1)
      write(905,*) 'Energy at ', ke(1)
      write(906,*) 'Energy at ', ke(1)
      write(907,*) 'Energy at ', ke(1)

      do a = 1,numa
      write(900,*) angles(a)
      write(901,*) angles(a)
      write(902,*) angles(a)
      write(903,*) angles(a)
      write(904,*) angles(a)
      write(905,*) angles(a)
      write(906,*) angles(a)
      write(907,*) angles(a)
      do i = 1,nch
            write(900,455) (efreal(a,i,j),j=1,nch)
            write(901,455) (efreal2(a,i,j),j=1,nch)
            write(902,455) (efim(a,i,j),j=1,nch)
            write(903,455) (efim2(a,i,j),j=1,nch)

            write(904,455) (ofreal(a,i,j),j=1,nch)
            write(905,455) (ofreal2(a,i,j),j=1,nch)
            write(906,455) (ofim(a,i,j),j=1,nch)
            write(907,455) (ofim2(a,i,j),j=1,nch)
      enddo
      enddo
      endif
      




455     format(20(1pe18.7))
999     format(i6,30(1pe14.7))


      deallocate (smr)
      deallocate (smi)
      deallocate (efreal)
      deallocate (efim)
      deallocate (efreal2)
      deallocate (efim2)
      deallocate (ofreal)
      deallocate (ofim)
      deallocate (ofreal2)
      deallocate (ofim2)
      deallocate (rmat)
      deallocate (imat)

556     continue



      deallocate (yn)
      deallocate (rjl)
      deallocate (rnl)
      deallocate (rjlp)
      deallocate (rnlp)
      deallocate (rki)
      deallocate (eta)
      deallocate (gamma)
      deallocate (erelk)

      deallocate (g)

      deallocate (jv)

      deallocate (erel)
      deallocate (cent)

      deallocate (unit)

      deallocate (fkmat)
      deallocate (smat)

      deallocate (ea)

      deallocate (impot)

      deallocate (alpha)
      deallocate (z)
      deallocate (z2)

      deallocate (lam)
      deallocate (coupindx)

      deallocate (pot)
      deallocate (ypot)
      deallocate (sigt)
      deallocate (sig)

      deallocate (r)

      deallocate (xi)
      deallocate (r2av)

      call MPI_FINALIZE(ierr)

      stop 'MOCC NORMAL EXITS'
      end
