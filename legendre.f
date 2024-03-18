c********************************************************************ZLB
c********************************************************************RWANG
c Ruihan Wang Sept 9th, 2019
c function to calculate Legendre polynomial
c
      subroutine legendre(m,n,x,p)
      implicit none

      integer m
      integer n

      integer i
      integer j
      double precision p(m,0:n)
      double precision x(m)

      if ( n .lt. 0 ) then
        return
      end if

      do i = 1,m
        do j = 0,n
            p(i,j) = 0.0d0
      enddo
      enddo


      do i = 1, m
        p(i,0) = 1.0D+00
      end do

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, m
        p(i,1) = x(i)
      end do

      do j = 2, n
        do i = 1, m

          p(i,j) = ( dble ( 2 * j - 1 ) * x(i) * p(i,j-1)
     &             - dble (     j - 1 ) *        p(i,j-2) )
     &             / dble (     j     )

         end do
      end do

      return
      end
