c********************************************************************ZLB
c 12-16-03 generalized, pcs
c  3-4-96 modified for all potentials to go to zero asymptotically
c  2-9-95 modified for SiH2+, two channels, pcs
c 11-5-94 modified by pcs
c         subroutine to return diabatic potential at
c         a given internuclear distance x.
c         Some editing of short range fits is required.
c         See below. If outgoing channels are not coulomb,
c         i.e., for a singly charged system, the long range
c         fits will also need modification.
c----------------------------------------------------------
c nch    number of channels
c ndat   number of potential data points per channel
c npar   number of partial waves per pass
c npot   size of potential array, includes diagonal
c        elements and upper triangle
c-----------------------------------------------------------
        subroutine dpot(x,v,r,pot,ypot,alpha,z,z2,
     1                  ea,delta,cinf,impot,nch,ndat,npar,npot,
     2                  xi,r2av)
      implicit real*8 (a-h,o-z)
        integer nch, ndat,npot,npar,i,j,k,l,impot,
     .        zc,zn,nprint
czlb        parameter(nch=3,ndat=31,npar=64)
czlb        parameter(npot=nch*(nch+1)/2)
czlb
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension r(ndat)
      dimension ea(nch),cinf(npot),impot(nch)
      dimension alpha(nch),z(nch),z2(nch)
      common/fit/rs,rl,nprint
czlb      06/27/04
      dimension xi(nch),r2av(nch)
czlb
c
c
c      common/dat/pot(ndat,npot),ypot(ndat,npot)
c      common/grid/r(ndat)
c      common/fit/rs,rl,nprint
c      common/fit1/ea(nch),delta,cinf(npot),impot(nch)
c      common/fit2/alpha(nch),z(nch),z2(nch)
czlb
c
        double precision v(npar,nch,nch),u(npot)
        double precision pott(ndat),ypott(ndat)
     .   ,y,r,pot,ypot,a,b,dem,delta,rs,rl,alpha,
     .    z,x,ea,z2,cinf
c
c rename spline fit data
c
       do 11 j=1,npot
       do 10 i=1,ndat
        ypott(i)=ypot(i,j)
10      pott(i)=pot(i,j)
c
c if internuclear distance is outside of spline range
c skip spline call
c
        if((x.gt.rl).or.(j.gt.nch .and. x.gt.rl)) then
          y=0.0d0
          goto 12
        end if
c
c spline fit all couplings and potentials
c
        call splint(r,pott,ypott,ndat,x,y)
12      u(j)=y
c
c get short range or long range fits
c
c short range potentials of form a*exp(-b*r) + Vmin
c pott(ndat)=Vmin is last data point
c  or asymptotic energy for coulomb channels  and
c  potential at equilibrium distance for incoming
c  channel. Adjust as required.
c
       if(x.lt.rs) then
c
c   potentials
c
         if(j.le.nch) then
         b=-(pott(2)-pott(1))/(r(2)-r(1))/(pott(1)-pott(impot(j)))
         a=(pott(1)-pott(impot(j)))*dexp(b*rs)
         u(j)=a*dexp(-b*x)+pott(impot(j))
c         write(20, *) 'coulomb',a,b,pott(ndat),u(j)
         end if
c
c
c short range form for diabatic off diagonal potential
c  which goes to zero at R=0. use form a*r^2+b*R.
c  coefficients determined by Cramers' rule.
         if(j.gt.nch) then
          dem=r(1)**2*r(2)-r(2)**2*r(1)
          a=(r(2)*pott(1)-r(1)*pott(2))/dem
          b=(r(1)**2*pott(2)-r(2)**2*pott(1))/dem
          u(j)=a*x**2+b*x
czlb
c           r0=r(1)*0.95
c           p0=pott(1)*0.95d0
c
c           a=(pott(1)/sin(0.5*r(1))-p0/sin(0.5*r0))/(r(1)-r0)
c           b=p0/sin(0.5*r0)-a*r0
c           u(j)=(a*x+b)*sin(0.5*x)
czlb

c        write(20, *) 'off0',a,b,dem, u(j)
         end if
       end if
c        write(21,*) 'x,u',x,u(j)
11      continue
c
c add constants so that all channels go to zero
c at R=infinity. delta is total electronic energy
c of channel 1 at infinity. ea(i) are excitation
c energies of each channel above channel 1 at
c R=infinity. Must use ab initio data.
c
c Also create v(k,i,j) potential matrix for indices i,j
c and partial wave k
c
        do 30 k=1,npar
c diagonal diabats
           do 32 i=1,nch
              v(k,i,i)=u(i)+delta-ea(i)
32         continue
c
c    now make long-range fits for the diabatic potentials
c
           zc=0
           if(x.gt.rl) then
              do 31 i=1,nch
c    check if Coulomb channel
                 if (z2(i).ne.0.0d0) then
                    v(k,i,i)=(z(i)*z2(i))/x
                    zc=zc+1
                    goto 31
                 else
                    zn=i
                 end if
c    for ion-neutral polarization potential
czlb              v(k,i,i)=-alpha(i)*z(i)**2/x**4/2.0d0
               v(k,i,i)=-alpha(i)*z(i)**2/x**4/2.0d0
     .                  +xi(i)*r2av(i)/x**3
31            continue
c    for ion-ion Coulomb channels
              if (zc.gt.0) then
czlb                 print*, 'Coulomb channels,zn=', zn,zc
              end if
c    for couplings
              do 36 j=1+nch,npot
cwy   cccccccccccccccccccccccccccccccccccccccccccccccccc
                if (abs(pot(ndat,j)).lt.1.d-8) then
                  u(j)=0.d0
cwy   cccccccccccccccccccccccccccccccccccccccccccccccccc
	          else
                  b=-(pot(ndat,j)-pot(ndat-1,j))
     .            /(r(ndat)-r(ndat-1))/pot(ndat,j)
c               b=0.0
czlb             a=(pot(ndat,j)-cinf(j))*dexp(b*rl)
czlb             u(j)=a*dexp(-b*x)+cinf(j)
c
                 a=(pot(ndat,j)-cinf(j))
                 u(j)=a*dexp(-b*(x-rl))+cinf(j)
cwy   ccccccccccccccccccccccccccccccccccccccccccccccccccc
	          end if
cwy   ccccccccccccccccccccccccccccccccccccccccccccccccccc

c              write(6, 40) j,a,b,pot(ndat,j),u(j)
36            continue
           end if
40         format(i3,1x,4e15.5)
c off-diagonal diabatic couplings
           l=nch+1
           do 33 i=1,nch-1
              do 34 j=i+1,nch
                v(k,i,j)=u(l)
                v(k,j,i)=u(l)
              l=l+1
34            continue
33         continue
30      continue
        if(nprint.eq.3) then
           write(60,999) x,(v(1,i,i),i=1,nch)
             do 35 i=1,nch-1
                write(60+i,999) x,(v(1,i,j),j=i+1,nch)
35           continue
c        write(41,999) x,v(1,1,2)
c        write(23,999) x,v(1,3,4)
         do ii=1,nch
         do jj=1,nch
           if(abs(v(1,ii,jj)).lt.1.0d-99) v(1,ii,jj) =0.0d0
         enddo
         enddo
         write(400,999) x,v(1,1,1),v(1,2,2),v(1,3,3),v(1,4,4),v(1,5,5),
     4                                               v(1,6,6),v(1,7,7)
         write(100,999) x,v(1,2,1),
     1                    v(1,3,1),v(1,3,2),
     2                    v(1,4,1),v(1,4,2),v(1,4,3)
         write(200,999) x,v(1,5,1),v(1,5,2),v(1,5,3),v(1,5,4),
     3                    v(1,6,1),v(1,6,2),v(1,6,3),v(1,6,4),v(1,6,5)
         write(300,999) x,v(1,7,1),v(1,7,2),v(1,7,3),v(1,7,4),v(1,7,5),
     4                                               v(1,7,6)
        end if
999     format(14g14.6)
        return
        end
