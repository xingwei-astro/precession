program main
implicit none
double precision x0, x1, y, k
integer i
x0=1.d0
do i=1, 100000
 call diff(x0,k)
 x1=x0-y(x0)/k
 if(abs(y(x1)).lt.1.d-6) then
  write(6,*) x1, y(x1)
  write(6,*) 2.d0/sqrt(1.d0+x1**2), 4.d0/sqrt(4.d0+x1**2)
  stop
 endif
 x0=x1
enddo
end program main

function y(x)
double precision x, y
y=2.d0/sqrt(1.d0+x**2) + 4.d0/sqrt(4.d0+x**2)-1.d0
end function y

subroutine diff(x,k)
double precision x, k, y
k=(y(x+1.d-6)-y(x-1.d-6))/(2.d0*1.d-6)
end subroutine diff