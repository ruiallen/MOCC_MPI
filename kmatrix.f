      subroutine kmatrix(ylog,rjl,rjlp,rnl,rnlp,jv,rki,nprint,
     1                   unit,fkmat,smat,parsig,nch,ndat,npar,
     1                   npot,smr,smi)
      implicit real*8 (a-h,o-z)
      integer nch,npar,nol,n4,nprint
czlb      parameter(nch=3,npar=64)
c shared arrays
       double precision ylog(npar,nch,nch),rjl(npar,nch,nch)
     1,rnl(npar,nch,nch),rjlp(npar,nch,nch),rnlp(npar,nch,nch)
     1,rki(nch)
      double precision jv(npar),unit,fkmat,smat,parsig
czlb
      dimension unit(nch,nch)
      dimension fkmat(npar,nch,nch),smat(npar,nch,nch),
     1          parsig(npar,nch,nch)
      dimension smr(npar,nch,nch),smi(npar,nch,nch)
c
c      common/une/unit(nch,nch)
c      common/fmtrix/ fkmat(npar,nch,nch),smat(npar,nch,nch)
c     . ,parsig(npar,nch,nch)
c       common/smatrix/ smr(npar,nch,nch),smi(npar,nch,nch)
czlb
c local arrays
       double precision fnum(npar,nch,nch),fden(npar,nch,nch),sum(npar)
     . ,yr(npar,nch,nch),yl(npar,nch,nch)
     . ,u(npar,nch,nch),smr,smi
     . ,term1ij,term2ij,x,sqs,pi
czlb       data nol/npar/,n4/npar/
       nol = npar
       n4  = npar
czlb
       pi=4.0d0*datan(1.0d0)
       n1=nch
       do 350 k=1,npar
       do 340 i=1,nch
       do 340 j=1,nch
       term1ij=0.0d0
       term2ij=0.0d0
       do 329 m=1,nch
      term1ij=term1ij+ylog(k,i,m)*rjl(k,m,j)-rjlp(k,i,m)*unit(m,j)
329   term2ij=term2ij+rnlp(k,i,m)*unit(j,m)-ylog(k,i,m)*rnl(k,m,j)
       fnum(k,i,j)=term1ij
340    fden(k,i,j)=term2ij
350    continue
       call gaussv(fden,fnum,sum,nch,nol,n1,n4)
c*****k-matrix=w.  now calculate the s-matrix.
       do 4 j=1,nch
       do 4 i=1,nch
       x=unit(i,j)
       do 6 k=1,npar
       fkmat(k,i,j)=fnum(k,i,j)
       yr(k,i,j)=x
6      sum(k)=x
       do 7 m=1,nch
       do 7 k=1,npar
7      sum(k)= sum(k) + fnum(k,i,m)*fnum(k,m,j)
       do 5 k=1,npar
5      u(k,i,j)=sum(k)
4      continue
       call  gaussv(u,yr,sum,nch,nol,n1,n4)

       do 12 j=1,nch
       do 12 i=1,nch
       do 13 k=1,npar
13     sum(k)=0.0d0
       do 14 m=1,nch
       do 14 k=1,npar
14     sum(k)= sum(k) + 2.0d0*yr(k,i,m)*fnum(k,m,j)
       do 15 k=1,npar
15     yl(k,i,j)=sum(k)
12     continue

       do 16 j=1,nch
       do 17 i=1,nch
       do 17 k=1,npar
17     yr(k,i,j) = -unit(i,j) +2.0d0*yr(k,i,j)
16     continue
c*****yr=re(s) and yl=im(s)

c...   calculate cross sections
       do 52 i=1,nch
       do 52 j=1,nch
       do 53 k=1,npar
       sqs=(unit(i,j)-yr(k,i,j))**2+yl(k,i,j)**2
       smat(k,i,j)=sqs
       smr(k,i,j)=yr(k,i,j)
       smi(k,i,j)=yl(k,i,j)
53     parsig(k,i,j)=(sqs)*pi*(2.0d0*jv(k)+1.0d0)/(rki(i)**2)
52     continue
c get output routine
       if(nprint.eq.2) call outdat(fkmat,smat,parsig,smr,smi,unit,
     1                             nch,ndat,npar,npot,jv)
       return
       end
