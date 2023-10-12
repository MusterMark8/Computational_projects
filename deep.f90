module fun
implicit none
contains

function pesi(h,d,n,i)
integer, intent(in) :: n,d,i
integer, dimension(d) :: h
real, dimension(n,h(i)) :: pesi
end function

function forward(a,b,n)
integer, intent(in) :: n
real,dimension(n) :: a, b
real :: forward

forward =tanh(dot_product(a,b))

end function

function backward(a,b,d,n)
integer, intent(in) :: n
real, dimension(n) :: a,b
real :: d, backward

backward= d*(1-(tanh(dot_product(a,b)))**2)

end function

function Err(f,r,n)
integer, intent(in) :: n
real, dimension(n) :: f,r
real :: Err

Err=.5*dot_product((f-r),(f-r))

end function

end module fun



Program lie
use fun
implicit none
!                    dim   esempi   strati
integer, parameter :: n=5, m=16, d=3
integer :: i, j, k
integer, dimension(d) :: h, delta
real :: ran, dom, u, v, delta1, delta2, ft, e
real, dimension(n,m) :: x
real, dimension(n,h) :: p1
real, dimension(m) :: r, f
real, dimension(h) :: w, y, yt
real, dimension(n) :: t

h = [3, 3, 1]

x(:,1)=[1,1,1,1,1]
x(:,2)=[1,1,1,-1,1]
x(:,3)=[1,1,-1,1,1]
x(:,4)=[1,1,-1,-1,1]
x(:,5)=[1,-1,1,1,1]
x(:,6)=[1,-1,1,-1,1]
x(:,7)=[1,-1,-1,1,1]
x(:,8)=[1,-1,-1,-1,1]
x(:,9)=[-1,1,1,1,1]
x(:,10)=[-1,1,1,-1,1]
x(:,11)=[-1,1,-1,1,1]
x(:,12)=[-1,1,-1,-1,1]
x(:,13)=[-1,-1,1,1,1]
x(:,14)=[-1,-1,1,-1,1]
x(:,15)=[-1,-1,-1,1,1]
x(:,16)=[-1,-1,-1,-1,1]

r=1

do i=1, m
	do j=1, n-1
		r(i)=r(i)*x(j,i)
	end do
end do


e=0.01

!INIZIALIZZAZIONE PESI
	
do k=1,d
	do j=1,n
		do i=1,h(k)-1
			call random_number(ran)
			pesi(h,d,n,k)(j,i)=ran-.5
		end do
	end do
end do


do j=1,n
	do i=1,h-1
		call random_number(ran)
			p(j,i)=ran-.5
	end do
end do

do k=1,h
	call random_number(ran)
		w(k)=ran-.5
end do 

!APPRENDIMENTO

	y(h) = 1  

do j=1, 10**6
	
	call random_number(ran) 
	i=int(m*ran) + 1
	do l=1,d
		do k=1,h(l)-1 
			y(k)=forward(x(:,i),p(:,k),n)
		end do
	f(i)=forward(y,w,h)
	delta2 = r(i)-f(i)
	!w=w+e*backward(w,y,delta2,h)*y
	delta1 = backward(w,y,delta2,h)
	do k=1,h-1
		p(:,k)=p(:,k)+e*backward(p(:,k),x(:,i),w(k)*delta1,n)*x(:,i)
	end do	
	w=w+e*backward(w,y,delta2,h)*y
	

	if( mod(j,100)==0 ) then
		write(7,*) j, Err(f,r,m)
		write(8,*) j, w
	end if	
	
	if( Err(f,r,m)< 0.1 ) then
		exit
	end if

	
end do

print*, Err(f,r,m)

!TEST
	
	yt(h)=1
	
do i=1,m
	t=x(:,i)
	print*, t(1), t(2), t(3), t(4)
	do k=1,h-1
		yt(k)=forward(t,p(:,k),n)
	end do
	print*, nint(forward(yt,w,h))
end do



!print*, NINT(ft)

!end do

end program lie