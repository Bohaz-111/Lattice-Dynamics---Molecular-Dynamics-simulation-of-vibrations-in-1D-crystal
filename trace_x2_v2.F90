program trace_x2_v2
      implicit real(8) (a-h,o-z)
!
!     ---------------------------------------------------------------
!     Compute Tr[x^2 exp(-beta H)] for
!       H = p^2/2 + 1/2 * Omega^2 * x^2 + g * x^4
!     ---------------------------------------------------------------
!
      allocatable :: h(:,:),e(:),w(:),q(:)
      real(8) :: qmin,qmax,omega0,g, temp,beta
      integer :: n, ierr, j, k
      real(8) :: Z, trace_x2, x2k, fac
!
!     Parameters as requested
      omega0 = sqrt(3.69d0)
      g      = 4.3d0
      temp   = 1.3d0
      beta   = 1.0d0/temp
!
      read (5,*) qmax, n
      qmin = -qmax
!
      allocate (h(n,n),e(n),w(n),q(n))
      call sindvr (0.5d0, qmin, qmax, n, h, w, q)
!
      do j = 1,n
         h(j,j) = h(j,j) + 0.5d0*(omega0*q(j))**2 + g*q(j)**4
      enddo
!
      call symevp (1, h, n, n, e, ierr)
      if (ierr .ne. 0) stop 'symevp failed'
!
!     Partition function and trace
      Z = 0.0d0
      trace_x2 = 0.0d0
      do k = 1, n
         fac = exp(-beta*e(k))
         Z = Z + fac
!        Eigenvector k is column k of h (real); x operator is diagonal: x(j)
         x2k = 0.0d0
         do j = 1, n
            x2k = x2k + (h(j,k)*h(j,k)) * (q(j)*q(j))
         enddo
         trace_x2 = trace_x2 + fac * x2k
      enddo
!
      write (6,'(1x,A,ES22.12)') 'Tr[x^2 e^{-beta H}] = ', trace_x2
      write (6,'(1x,A,ES22.12)') '<x^2> (normalized)   = ', trace_x2 / Z
!
      deallocate (h,e,w,q)
      stop

end program trace_x2_v2

      subroutine sindvr (on2m,a,b,n,t,w,x)
      implicit real(8) (a-h,o-z)
!
!     ----------------------------------------------------------------- 
!     n point sine DVR in a < x < b
!     D.T.Colbert and W.H.Miller, JCP 96 (1992) 1982, eq.(A6)
!     ----------------------------------------------------------------- 
!
      dimension t(n,n),w(n),x(n)
      dimension s(2*n)
!
      m = n+1
      pi = acos(-1.0d0)
      alfa = 0.5d0*on2m*(pi/(b-a))**2
      beta = pi/(2*m)
      cosa = 1.0d0
      sina = 0.0d0
      cosb = cos(beta)
      sinb = sin(beta)
      dx = (b-a)/m
      wt = sqrt(dx)
      t0 = alfa*(2*m**2+1)/3.0d0
      do k = 1,2*n
         alfa = -alfa
         temp = cosa*sinb+sina*cosb
         cosa = cosa*cosb-sina*sinb
         sina = temp
         s(k) = alfa/sina**2
      enddo
      do j = 2,n
         do i = 1,j-1
            t(i,j) = s(j-i)-s(j+i)
            t(j,i) = t(i,j)
         enddo
      enddo
      do j = 1,n
         t(j,j) = t0-s(j+j)
         w(j) = wt
         x(j) = a+j*dx
      enddo
      return
      end 

      subroutine symevp (iv,a,lda,n,d,ierr)
      implicit real(8) (a-h,o-z)
!
!     ----------------------------------------------------------------- 
!     Real symmetric matrix eigenvalue problem.
!
!     On entry: 
!     a contains a real symmetric matrix. 
!     only the LOWER triangle is used.
!     
!     On return: 
!     d contains the eigenvalues of a in ascending order.
!     a is overwritten with the corresponding eigenvectors if iv>0.
!     ----------------------------------------------------------------- 
!
      dimension a(lda,n),d(n)
      allocatable :: e(:)
!
      allocate (e(n))
      call rstred (a,lda,n,d,e)
      if (iv .gt. 0) then
         call rstqlv (d,e,n,a,lda,ierr)
      else
         call rstqli (d,e,n,ierr)
      endif
      deallocate (e)
      return
      end

      subroutine rstred (a,lda,n,d,e)
      implicit real(8) (a-h,o-z)
!
!     ----------------------------------------------------------------- 
!     Reduction of a real symmetric matrix to tridiagonal form.
!
!     Taken from Numerical Recipes.
!     ----------------------------------------------------------------- 
!
      dimension a(lda,n),d(n),e(n)
!
      do i = n,2,-1
         l = i-1
         h = 0.d0
         scale = 0.d0
         if (l .gt. 1) then
            do k = 1,l
               scale = scale+abs(a(i,k))
            enddo
            if (scale .le. 0.d0) then
               e(i) = a(i,l)
            else
               do k = 1,l
                  a(i,k) = a(i,k)/scale
                  h = h+a(i,k)**2
               enddo
               f = a(i,l)
               g = -sign(sqrt(h),f)
               e(i) = scale*g
               h = h-f*g
               a(i,l) = f-g
               f = 0.d0
               do j = 1,l
                  a(j,i) = a(i,j)/h
                  g = 0.d0
                  do k = 1,j
                     g = g+a(j,k)*a(i,k)
                  enddo
                  do k = j+1,l
                     g = g+a(k,j)*a(i,k)
                  enddo
                  e(j) = g/h
                  f = f+e(j)*a(i,j)
               enddo
               hh = f/(h+h)
               do j = 1,l
                  f = a(i,j)
                  g = e(j)-hh*f
                  e(j) = g
                  do k = 1,j
                     a(j,k) = a(j,k)-f*e(k)-g*a(i,k)
                  enddo
               enddo
            endif
         else
            e(i) = a(i,l)
         endif
         d(i) = h
      enddo
      d(1) = 0.d0
      do i = 1,n
         l = i-1
         if (abs(d(i)) .gt. 0.d0) then
            do j = 1,l
               g = 0.d0 
               do k = 1,l
                  g = g+a(i,k)*a(k,j)
               enddo
               do k = 1,l
                  a(k,j) = a(k,j)-g*a(k,i)
               enddo
            enddo
         endif
         d(i) = a(i,i)
         a(i,i) = 1.d0
         do j = 1,l
            a(i,j) = 0.d0
            a(j,i) = 0.d0
         enddo
      enddo
      e(1) = 0.d0
      return
      end

      subroutine rstqli (d,e,n,ierr)
      implicit real(8) (a-h,o-z)
!
!     ----------------------------------------------------------------- 
!     Eigenvalues of a real symmetric tridiagonal matrix. 
!
!     Taken from Numerical Recipes.
!     ----------------------------------------------------------------- 
!
      dimension d(n),e(n)
!
      do i = 2,n
         e(i-1) = e(i)
      enddo
      e(n) = 0.d0
      do l = 1,n
         iter = 0
   1     do m = l,n-1
            dd = abs(d(m))+abs(d(m+1))
            ee = abs(e(m))
            if (ee+dd .le. dd) goto 2
         enddo
         m = n
   2     if (m .ne. l) then
            if (iter .eq. 30) then
               ierr = 1
               return
            endif 
            iter = iter+1
            g = (d(l+1)-d(l))/(2.d0*e(l))
            r = sqrt(1.d0+g**2)
            g = d(m)-d(l)+e(l)/(g+sign(r,g))
            s = 1.d0
            c = 1.d0
            p = 0.d0
            do i = m-1,l,-1
               f = s*e(i)
               b = c*e(i)
               r = sqrt(f**2+g**2)
               e(i+1) = r
               if (abs(r) .le. 0.d0) then
                  d(i+1) = d(i+1)-p
                  e(m) = 0.d0
                  goto 1
               endif
               s = f/r
               c = g/r
               g = d(i+1)-p
               r = (d(i)-g)*s+2.d0*c*b
               p = s*r
               d(i+1) = g+p
               g = c*r-b
            enddo
            d(l) = d(l)-p
            e(l) = g
            e(m) = 0.d0
            goto 1
         endif
      enddo
      do j = 1,n-1
         k = j
         do i = j+1,n
            if (d(i) .lt. d(k)) k = i
         enddo
         if (k .ne. j) then
            swap = d(k)
            d(k) = d(j)
            d(j) = swap
         endif
      enddo
      ierr = 0
      return
      end

      subroutine rstqlv (d,e,n,v,ldv,ierr)
      implicit real(8) (a-h,o-z)
!
!     ----------------------------------------------------------------- 
!     Eigenvalues and eigenvectors of a real symmetric tridiagonal
!     matrix. (v must be initialized to a unit matrix on entry if
!     this routine is not called by symevp.)
!
!     Taken from Numerical Recipes.
!     ----------------------------------------------------------------- 
!
      dimension d(n),e(n),v(ldv,n)
!
      nv = min(ldv,n)
      do i = 2,n
         e(i-1) = e(i)
      enddo
      e(n) = 0.d0
      do l = 1,n
         iter = 0
   1     do m = l,n-1
            dd = abs(d(m))+abs(d(m+1))
            ee = abs(e(m))
            if (ee+dd .le. dd) goto 2
         enddo
         m = n
   2     if (m .ne. l) then
            if (iter .eq. 30) then
               ierr = 1
               return
            endif
            iter = iter+1
            g = (d(l+1)-d(l))/(2.d0*e(l))
            r = sqrt(1.d0+g**2)
            if (abs(g-r) .gt. abs(g+r)) r = -r
            g = d(m)-d(l)+e(l)/(g+r)
            s = 1.d0
            c = 1.d0
            p = 0.d0
            do i = m-1,l,-1
               f = s*e(i)
               b = c*e(i)
               r = sqrt(f**2+g**2)
               e(i+1) = r
               if (abs(r) .le. 0.d0) then
                  d(i+1) = d(i+1)-p
                  e(m) = 0.d0
                  goto 1
               endif
               s = f/r
               c = g/r
               g = d(i+1)-p
               r = (d(i)-g)*s+2.d0*c*b
               p = s*r
               d(i+1) = g+p
               g = c*r-b
               do k = 1,nv
                  f = v(k,i+1)
                  v(k,i+1) = s*v(k,i)+c*f
                  v(k,i) = c*v(k,i)-s*f
               enddo
            enddo
            d(l) = d(l)-p
            e(l) = g
            e(m) = 0.d0
            goto 1
         endif
      enddo
      do j = 1,n-1
         k = j
         do i = j+1,n
            if (d(i) .lt. d(k)) k = i
         enddo
         if (k .ne. j) then
            swap = d(k)
            d(k) = d(j)
            d(j) = swap
            do i = 1,nv
               swap = v(i,k)
               v(i,k) = v(i,j)
               v(i,j) = swap
            enddo
         endif
      enddo
      ierr = 0
      return
      end


