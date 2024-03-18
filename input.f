c********************************************************************ZLB
c 12-16-03 generalized, pcs
c 11-5-94 modified by pcs
c         subroutine to input diabatic potential matrix
c         from fort.17. Note required serial structure
c         of fort.17. Structure is diagonal elements first
c         (1,1 2,2 3,3 ...) followed by off-diagonal elements
c         (1,2 1,3 ... 2,3 ...). Subroutine also determines
c         spline fit.
c-----------------------------------------------------------
c ndat   number of potential data points per channel
c npot   size of potential array, includes diagonal
c        elements and upper triangle
c pot    array of diabatic potential data
c r      vector of internuclear distances
c ypot   array of spline data
c---------------------------------------------------
        subroutine input(r,pot,ypot,alpha,z,z2,nch,ndat,npar,npot,
     .                   xi,r2av)
        implicit real*8 (a-h,o-z)
        integer ndat, npot,nch,zc,zn
czlb        parameter(ndat=31, nch=3)
czlb        parameter(npot=nch*(nch+1)/2)
c input reads in x, y values spline nodes
czlb
        dimension pot(ndat,npot),ypot(ndat,npot)
        dimension r(ndat)
        dimension alpha(nch),z(nch),z2(nch)
        dimension potm(nch,nch)
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb
c
c        common/dat/pot(ndat,npot),ypot(ndat,npot)
c        common/grid/r(ndat)
c        common/fit2/alpha(nch),z(nch),z2(nch)
czlb
        double precision pott(ndat),ypott(ndat),
     .     yp1(npot),yp2(npot),z,alpha,z2
     .    ,r,pot,ypot,ypp2,ypp1,dumb,tmp
c
c set all first and last point derivatives to zero
c
        do 223 j=1,npot
           yp1(j)=0.0d0
           yp2(j)=0.0d0
223     continue
c
c **** read potential (fort.17) and coupling (fort.16)
c **** data, may require modification
czlb
czlb read potentials from file <pot.diabatic>
czlb
        do 23 i=1,ndat
          read(16,*)
          read(16,100) r(i)
          do j=1,nch
            read(16,*) (potm(j,k),k=1,j)
          enddo
czlb
czlb make matrix symmatric
czlb
          do j=1,nch
            do k=1,j
              potm(k,j)=potm(j,k)
            enddo
          enddo
c
          do j=1,nch
            pot(i,j)=potm(j,j)
          enddo
c
          jmin=nch+1
          jmax=nch+(nch-1)
          do j=1,nch
            do k=j+1,nch
              pot(i,jmin+k-j-1)=potm(j,k)
c              write(*,*) jmin+k-j-1
            enddo
            jmin0=jmin
            jmax0=jmax
            jmin=jmax0+1
            jmax=jmax0+(jmax0-jmin0)
          enddo
c
czlb           read(17,*)r(i),(pot(i,j),j=1,nch)
czlb           read(16,*)r(i),(pot(i,j),j=1+nch,npot)
c        write(51,*) r(i), pot(i,1),pot(i,2),pot(i,3)
23      continue
100     format(40x,E16.8)
c
c determine first point derivatives for all potentials
c and couplings and last point for couplings
c
        do 24 j=1,npot
           yp1(j)=(pot(2,j)-pot(1,j))/(r(2)-r(1))
           yp2(j)=(pot(ndat,j)-pot(ndat-1,j))
     .            /(r(ndat)-r(ndat-1))
24      continue
c
c determine last point derivatives for potentials only
c using asymptotic forms
c

        zc=0
        do 53 j=1,nch
c    check if Coulomb channel
           if (z2(j).ne.0.0d0) then
              yp2(j)=-(z(j)*z2(j))/r(ndat)**2
              zc=zc+1
c              print*,'j=',j,'zc=',zc
              goto 53
           else
              zn=j

           end if
c    for ion-neutral polarization potential
czlb           yp2(j)=2.0d0*alpha(j)*z(j)**2/r(ndat)**5
           yp2(j)=2.0d0*alpha(j)*z(j)**2/r(ndat)**5
     .           -3.0*xi(j)*r2av(j)/r(ndat)**4
c           print*,'dev',alpha(j),z(j),yp2(j)
 53     continue
c    number of  ion-ion Coulomb channels


c223     continue
c
c get spline coefficients
c for pot1
        do 52 j=1,npot
        ypp1=yp1(j)
        ypp2=yp2(j)
        do 50 i=1,ndat
50      pott(i)=pot(i,j)
        call spline(r,pott,ndat,ypp1,ypp2,ypott)
        do 51 i=1,ndat
51      ypot(i,j)=ypott(i)
52      continue

        return
c pot(i,j), ypot(i,j) contains the spline parameters at node i
c potential j.
        end
