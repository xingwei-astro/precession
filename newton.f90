program main
implicit none
double precision pi, Pr
double precision x0, x1, y, k
integer i
pi=acos(-1.d0)
Pr=0.1
x0=1.d0
do i=1, 100000
 call diff(x0,k)
 x1=x0-y(x0)/k
 if(abs(y(x1)).lt.1.d-6) then
  write(6,*) x1, y(x1), x1/sqrt(2.d0)
!  write(6,*) 2*pi/sqrt(pi**2+x1**2), -4*pi/sqrt(4*pi**2+x1**2)
  write(6,*) 2*pi/sqrt(pi**2+x1**2)*sqrt((1+Pr)/(1-Pr)), &
             -4*pi/sqrt(4*pi**2+x1**2)*sqrt((1+Pr)/(1-Pr))
  stop
 endif
 x0=x1
enddo
end program main

function y(x)
double precision x, y
double precision pi, Pr
pi=acos(-1.d0)
Pr=0.1
!y = 2*pi/sqrt(pi**2+x**2) + 4*pi/sqrt(4*pi**2+x**2) - 1
y = 2*pi/sqrt(pi**2+x**2)*sqrt((1+Pr)/(1-Pr)) &
  + 4*pi/sqrt(4*pi**2+x**2)*sqrt((1+Pr)/(1-Pr)) - 1
end function y

subroutine diff(x,k)
double precision x, k, y
k=(y(x+1.d-6)-y(x-1.d-6))/(2.d0*1.d-6)
end subroutine diff