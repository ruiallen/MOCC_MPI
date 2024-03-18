c********************************************************************ZLB
c11-5-94  modified by pcs
c ylog is the vectorized log-deritative integrator,
c RSTART is starting value
c YN is the NxN log-deritative at RN. H is the step size.
c local arrays
c--------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c-------------------------------------------------------
      subroutine ylog(my_id,h,rstart,rn,yn,erel,cent,unit,lam,coupindx,
     1                r,pot,ypot,alpha,z,z2,
     2                ea,delta,cinf,impot,nch,ndat,npar,npot,
     3                xi,r2av,jstart)
      implicit real*8 (a-h,o-z)
      integer nch,npar,nol,n4,n1,m,max,my_id,i
czlb      parameter(nch=3,npar=64,ndat=31)
czlb      parameter(npot=nch*(nch+1)/2)
      double precision work1(npar,nch,nch),work2(npar,nch,nch)
     1,work3(npar,nch,nch),zne(npar,nch,nch),zno(npar,nch,nch),sum(npar)
c shared arrays
      double precision unit,re,ro,rn,rstart,h
      double precision v(npar,nch,nch),yn(npar,nch,nch)
      double precision lam,jstart
      integer coupindx
      dimension lam(nch),coupindx(nch,nch)
      dimension pot(ndat,npot),ypot(ndat,npot)
      dimension r(ndat)
      dimension alpha(nch),z(nch),z2(nch)
      dimension ea(nch),cinf(npot),impot(nch)
czlb        06/27/04
      dimension xi(nch),r2av(nch)
czlb
      dimension unit(nch,nch)
      dimension erel(nch),cent(npar,nch)
c      common/une/unit(nch,nch)
czlb      data nol/npar/,n4/npar/
      nol = npar
      n4  = npar
czlb

       n1=nch
       do 150 i=1,nch
       do 150 j=1,nch
       do 149 k=1,npar
149    zne(k,i,j)=unit(i,j)*1.0d20
150    continue
       re=rstart
       max=idnint((rn-rstart)/(2.0d0*h))
       do 200 m=1,max
       ro=re+h
c	if(ro.gt.50.d0) go to 200



       call potval(ro,h,v,erel,cent,unit,lam,coupindx,
     1             r,pot,ypot,alpha,z,z2,
     2             ea,delta,cinf,impot,nch,ndat,npar,npot,
     3             xi,r2av,jstart)
       do 160 i=1,nch
       do 160 j=1,nch
       do 159 k=1,npar
       work1(k,i,j)=unit(i,j)
       work2(k,i,j)=unit(i,j)
       work3(k,i,j)=unit(i,j)
       zno(k,i,j)=unit(i,j)+zne(k,i,j)
159    work3(k,i,j)=unit(i,j)+v(k,i,j)
160    continue
       call gaussv(zno,work1,sum,nch,nol,n1,n4)
       call gaussv(work3,work2,sum,nch,nol,n1,n4)

	 do 171 i=1,nch
       do 171 j=1,nch
       do 168 k=1,npar
168    work3(k,i,j)=0.0d0
       do 170 n=1,nch
       do 169 k=1,npar
169    work3(k,i,j)=work3(k,i,j)+zne(k,i,n)*work1(k,n,j)-v(k,i,n)*
     1work2(k,n,j)*8.0d0

170    continue
       do 172 k=1,npar
172    zno(k,i,j)=work3(k,i,j)
171    continue
       re=ro+h


       call potval(re,h,v,erel,cent,unit,lam,coupindx,
     1             r,pot,ypot,alpha,z,z2,
     2             ea,delta,cinf,impot,nch,ndat,npar,npot,
     3             xi,r2av,jstart)

       do 180 i=1,nch
       do 180 j=1,nch
       do 179 k=1,npar
       work1(k,i,j)=unit(i,j)
       work2(k,i,j)=unit(i,j)
179    zne(k,i,j)=unit(i,j)+zno(k,i,j)
180    continue

       call gaussv(zne,work1,sum,nch,nol,n1,n4)

       do 191 i=1,nch
       do 191 j=1,nch
       do 188 k=1,npar
188    work2(k,i,j)=0.0d0
       do 190 n=1,nch
       do 189 k=1,npar
189    work2(k,i,j)=work2(k,i,j)+zno(k,i,n)*work1(k,n,j)
     1-4.0d0*v(k,i,n)*unit(n,j)
190    continue

       do 192 k=1,npar
192    zne(k,i,j)=work2(k,i,j)
191    continue

200    continue

c last integration point

        do 201 i=1,nch
        do 201 j=1,nch
        do 202 k=1,npar
        zne(k,i,j)=(zne(k,i,j)+2.0d0*v(k,i,j))/h
202     yn(k,i,j)=zne(k,i,j)
201     continue
c       write(21,*) yn(1,1,1),yn(1,2,2)
c       write(21,*) yn(1,1,2),yn(1,2,1)
       return
       end
