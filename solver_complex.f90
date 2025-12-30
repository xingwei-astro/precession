! gfortran -o solver_complex solver_complex.f90 nag.f

program main
implicit none
integer n
parameter (n=2)
double complex a(n,n), a_inv(n,n), b(n), c(n)

a(1,1)=(1.d0, 0.d0)
a(1,2)=(0.d0, 1.d0)
a(2,1)=(0.e0, -1.d0)
a(2,2)=(-1.d0, 1.d0)
call mat_inv(n,a,a_inv)
write(6,*) a_inv
b(1)=(1.d0, 0.d0)
b(2)=(1.d0, -1.d0)
call mat_mul(n,a_inv,b,c)
write(6,*) c
end program main

subroutine mat_inv(n,a,b)
implicit none
integer i, j, n
double complex a(n,n), b(n,n)
double precision RINT(n), DETR, DETI
integer IDETE, IFAIL
do i=1, n
 do j=1, n
  b(i,j)=(0.d0, 0.d0)
  if(i.eq.j) b(i,j)=(1.d0, 0.d0)
 enddo
enddo
call F03AHF(n,a,n,DETR,DETI,IDETE,RINT,IFAIL)
call F04AKF(n,n,a,n,RINT,b,n)
end subroutine mat_inv

subroutine mat_mul(n,a,b,c)
implicit none
integer i, k, n
double complex a(n,n), b(n), c(n)
do i=1,n
 c(i)=0.d0
 do k=1,n
  c(i)=c(i)+a(i,k)*b(k)
 enddo
enddo
end subroutine mat_mul