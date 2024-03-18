c********************************************************************ZLB
c 11-5-94  modified for arbitrary number of channels, pcs
c      output kmatrix and S*S matrix
c      modify format statements 39 and 45 as required.
c---------------------------------------------------------
c nch    number of channels
c npar   number of partial waves per pass
c--------------------------------------------------------
       subroutine outdat(fkmat,smat,parsig,smr,smi,unit,
     1                   nch,ndat,npar,npot,jv)
       implicit real*8 (a-h,o-z)
       integer nch,npar
czlb       parameter(nch=3,npar=64)
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
       double precision  jv(1),unit,fkmat,smat,parsig,smr,smi

       do 30 k=1,npar
       write(2,33)jv(k)
       write(2,36)
       do 29 i=1,nch
       write(2,34)(fkmat(k,i,j),j=1,nch)
29       continue
       write(2,35)
       do 28 i=1,nch
      write(2,34)(smat(k,i,j),j=1,nch)
28        continue
c
       write(2,44)
       do 18 i=1,nch
       write(2,34)(smr(k,i,j),j=1,nch)
18        continue
c
       write(2,46)
       do 19 i=1,nch
       write(2,34)(smi(k,i,j),j=1,nch)
19        continue
30     continue
33     format( /,' J= ',f10.1)
34     format(10g12.3)
35     format( ' - - - - - T*T MATRIX - - - - - -')
36      format( ' - - - - - K MATRIX - - - - - - -')
44     format( ' - - - - - Re S MATRIX - - - - - -')
46      format( ' - - - - - Im S MATRIX - - - - - - -')
c output charge transfer cross section
c sigma(3-1),sigma(3-2)
       write(4,39)
39     format( '                        ')

       do 40 k=1,npar

       write(4,33)jv(k)
c********writes may require modification*************
       write(4,41)(smat(k,nch,i),i=1,nch)
c     %,(parsig(k,nch,i),i=nch-1,1,-1)
c output graph data
c       write(7,42)jv(k),(parsig(k,nch,i),i=nch-1,1,-1)
40     continue
c *********formats may require modification **********
41     format( ' S*S ',10g12.4)
c42     format(6g15.5)

c45     format( ' TOTAL ',/, ' SIGMA(2-1)'
c     %,g15.4,
c     %' TOTAL ',g15.4)
       return
       end
