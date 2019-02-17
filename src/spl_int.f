C FILE: SPL_INT.F
C f2py -c -m spl_int  spl_int.f
C import spl_int
C print spl_int.spl_evf.__doc__
      SUBROUTINE SPL_INT(xi,yi,n,k)
C
C     CALCULATE THE CURVATURE
C     Cubic spline representation for arbitary knot spacing.
C
      INTEGER N
      DOUBLE PRECISION x(N),y(N),k(N),u(N-2),v(N-2),xi(N),yi(N)
Cf2py intent(in) xi,yi,n,
Cf2py intent(out) xi,yi,k
Cf2py depend(n) xi,yi,k
      do 100 i=1,N-1
         x(i) = xi(i+1) - xi(i)
         y(i) = (yi(i+1) - yi(i))/x(i)
 100  enddo
c     # Gaussian Elimination
      u(1) = 2 * (x(2) + x(1))
      v(1) = 6 * (y(2) - y(1))
      do 200 i=2,N-2
         u(i) = 2 * (x(i) + x(i+1)) - x(i)**2/u(i-1)
         v(i) = 6 * (y(i+1) - y(i)) - x(i)*v(i-1)/u(i-1)
 200  enddo
c     # Back-substitution
c     Natural spline
      k(N) = 0
      k(1) = 0
      do 300 i=N-2,1,-1
         k(i+1) = (v(i)-x(i+1)*k(i+2)) / u(i)
 300  enddo
      END

      SUBROUTINE SPL_INTf(xi,yi,n,b,k,d)
C
C     CALCULATE THE CURVATURE fast
C     Returns all coefficients for each interval
C
      INTEGER N
      DOUBLE PRECISION xi(N),yi(N)
      DOUBLE PRECISION h(N),dy(N),b(N-1),k(N),d(N-1),u(N-2),v(N-2)
Cf2py intent(in) xi,yi,n
Cf2py intent(out) xi,yi,b,k,d
Cf2py depend(n) xi,yi,k
Cf2py depend(n-1) b,d
      do 100 i=1,N-1
         h(i) = xi(i+1) - xi(i)
         dy(i) = (yi(i+1) - yi(i))/h(i)
 100  enddo
c     # Gaussian Elimination
      u(1) = 2 * (h(2) + h(1))
      v(1) = 6 * (dy(2) - dy(1))
      do 200 i=2,N-2
         u(i) = 2 * (h(i) + h(i+1)) - h(i)**2/u(i-1)
         v(i) = 6 * (dy(i+1) - dy(i)) - h(i)*v(i-1)/u(i-1)
 200  enddo
c     # Back-substitution
      k(N) = 0
      do 300 i=N-2,1,-1
         k(i+1) = (v(i)-h(i+1)*k(i+2))/u(i)
         b(i+1) = -k(i+2)*h(i+1)/6. - k(i+1)*h(i+1)/3 + dy(i+1)
C         b(i+1) = dy(i+1)
         d(i+1) = (k(i+2)-k(i+1))/6/h(i+1)
 300  enddo
      k(1) = 0
c      b(1) = -k(2)*h(1)/6. - k(1)*h(1)/3 + dy(1)/h(1)
      b(1) = dy(1)
      d(1) = (k(2)-k(1))/6/h(1)
C     y'(x) = b + 2cx + 3dx**2
C     b(N) = b(N-1) + k(N-1)*h(N-1) + 3*d(N-1)*h(N-1)**2
      END

      SUBROUTINE SPL_EQ_INT(xi,yi,n,k)
C
C     CALCULATE THE CURVATURE
C     Uniform knot spacing
C
      INTEGER N
      DOUBLE PRECISION x(N),y(N),k(N),u(N-2),v(N-2),xi(N),yi(N)
Cf2py intent(in) xi,yi,n,
Cf2py intent(out) xi,yi,k
Cf2py depend(n) xi,yi,k
      
      do 100 i=1,N-1
C         x(i) = 1 ! unit grid xi(i+1) - xi(i)
         y(i) = yi(i+1) - yi(i)
 100  enddo
c     # Gaussian Elimination
      u(1) = 4     !2 * (x(2) + x(1))
      v(1) = 6 * (y(2) - y(1))
      do 200 i=2,N-2
         u(i) = 4 - 1. /u(i-1) ! 2 * (x(i) + x(i+1)) - x(i)**2/u(i-1)
         v(i) = 6 * (y(i+1) - y(i)) - v(i-1)/u(i-1)
 200  enddo
c     # Back-substitution
      k(N) = 0
      k(1) = 0
      do 300 i=N-2,1,-1
c         k(i+1) = (v(i)-x(i+1)*k(i+2))/u(i)
         k(i+1) = (v(i)-k(i+2))/u(i)
 300  enddo
      END

      SUBROUTINE SPL_EV(x,y,k,n,xx,nn,yy)
C
C     CALCULATE THE SPLINE for equidistant grid
C
C     must be sorted
      INTEGER N,nn
      DOUBLE PRECISION x(N),y(N),k(N),xx(nn),yy(nn),h,dx
Cf2py intent(in) x,y,n,k
Cf2py intent(out) yy
Cf2py depend(n) x,y,k
Cf2py depend(nn) xx,yy
      m=1
c      write(6,*) xx(m), x(1)
      do 100 i=1,nn
         ! find the left knot x
         do while (x(m) .le. xx(i) .and. m .lt. N)
            m=m+1
         enddo
         m=m-1
         h = x(m+1) - x(m) !x(i+1)-x(i)
         dx = (xx(i)-x(m))
         yy(i) = y(m) + dx* (-k(m+1)*h/6. - k(m)*h/3 + (y(m+1)-y(m))/h
     1                       + dx * (k(m)/2 + dx * (k(m+1)-k(m))/6/h))
c         if (i .lt. 10) write(6,*) m,yy(i),xx(i),k(m),dx
 100  enddo
      END

      SUBROUTINE SPL_EVf(x,y,b,k,d,n,xx,nn,yy)
C
C     CALCULATE THE SPLINE
C
C     must be sorted
      INTEGER N,nn
      DOUBLE PRECISION x(N),y(N),b(N-1),k(N),d(N-1),xx(nn),yy(nn)
      DOUBLE PRECISION dx
Cf2py intent(in) x,y,b,k,d,n,xx,nn
Cf2py intent(out) yy
Cf2py depend(n) x,y,k, b,d
Cf2py depend(nn) xx,yy
      m = 1
c      write(6,*) xx(m), x(1)
      do 100 i=1,nn
         ! find the left knot x
         do while (x(m+1) .lt. xx(i) .and. m .lt. N-1)
            m = m + 1
         enddo
c         m = m - 1
c         h = x(m+1) - x(m) !x(i+1)-x(i)
         dx = xx(i) - x(m)
         yy(i) = y(m) + dx * (b(m) + dx * (k(m)/2 + dx * d(m)))
c         if (i .lt. 10) write(6,*) 'ok',m,yy(i),b(m),
c     1   (y(m+1)-y(m))/h

c         if (i .lt. 10) write(6,*) m,yy(i),xx(i),k(m),dx
 100  enddo
      END

      SUBROUTINE SPL_EQ_EV(x,y,k,n,xx,nn,yy)
C
C     CALCULATE THE SPLINE for equidistant grid
C
C     must be sorted
      INTEGER N,nn
      DOUBLE PRECISION x(N),y(N),k(N),xx(nn),yy(nn),h,dx
Cf2py intent(in) x,y,n,k
Cf2py intent(out) yy
Cf2py depend(n) x,y,k
Cf2py depend(nn) xx,yy
      h = x(2) - x(1) !x(i+1)-x(i)
      do 100 i=1,nn
         dx = (xx(i)/h-i+1)
         yy(i) = y(i) + dx* (-k(i+1)/6. - k(i)/3 + (y(i+1)-y(i)) +
     1                        dx * (k(i)/2 + dx * (k(i+1)-k(i))/6))
 100  enddo
      END
C END FILE SPL_INT.F
