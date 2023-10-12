module fun
implicit none
contains

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

Err=.5*dot_product((f-r),(f-r))/n*100

end function

end module fun



program conn
use fun
implicit none

!q nueroni hidden, n matrice 10x10 +1 (soglia), m esempi

integer, parameter :: q=7, n=50, m=16000, t=4000, h=10, esp=6   
integer :: i, j, k, a, b, c, z, ai, aj
real :: ran, delta1, delta2, e=0.001, oldE
real, dimension(q,q) :: es
real, dimension(n,m) :: x
real, dimension(n,t) :: tx
real, dimension(n,h-1) :: p
real, dimension(m) :: r, f
real, dimension(t) :: tr, tf
real, dimension(h) :: w, y, ty, yt
real, dimension(5,15) :: ma

	x(n,:)=1
	tx(n,:)=1
	
do j=1,m
	if(j>m/2) then
		r(j)=-1
	else
		r(j)=1
	end if
end do

do j=1,t
	if(j>t/2) then
		tr(j)=-1
	else
		tr(j)=1
	end if
end do

open(13,file="fort.13",status='old')
do i=1,q
	do j=1,q
		read(13,*) x(q*(i-1)+j,:)
	end do
end do

close(13)
!print*, x(:,2501)
	
open(14,file="fort.14",status='old')
do i=1,q
	do j=1,q
		read(14,*) tx(q*(i-1)+j,:)
	end do
end do

close(14)

!CROSS VALIDATION

!open(1,file="AAAvalid.txt",status='old')
!do i=1,15

!		read(1,*) ma(:,i)

!end do

!close(1)

!do i=1,15

!		write(3,*) ma(:,i)

!end do
	
!inizializzare i pesi
	
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
	
	y(h)=1  
	ty(h)=1

	print*, 'inizio'
	
do j=1, 10**esp


	
	call random_number(ran) 
	i=int(m*ran) + 1
	!print*, i
	do k=1,h-1  
		y(k)=forward(x(:,i),p(:,k),n)
	end do
	f(i)=forward(y,w,h)
	delta2 = r(i)-f(i)
	delta1 = backward(w,y,delta2,h)
	w=w+e*backward(w,y,delta2,h)*y
	do k=1,h-1
		p(:,k)=p(:,k)+e*backward(p(:,k),x(:,i),w(k)*delta1,n)*x(:,i)
	end do	
	!w=w+e*backward(w,y,delta2,h)*y
	
	

	if( mod(j,int(10**(esp-2)))==0 ) then
		
		oldE=Err(f,r,m)
		
		do i=1,m
		
			do k=1,h-1  
				y(k)=forward(x(:,i),p(:,k),n)
			end do
			f(i)=forward(y,w,h)
			
		end do
		
		do i=1,t
		
			do k=1,h-1  
				ty(k)=forward(tx(:,i),p(:,k),n)
			end do
			tf(i)=forward(ty,w,h)
			
		end do
		
		write(1,*) j, Err(f,r,m)
		write(2,*) j, Err(tf,tr,t)                       !attenzione all'ordine
	end if	
	
	if( Err(f,r,m)<0.1 ) then
		print*, 'err', j
		exit
	end if
	
	!if(e**2*dot_product(backward(w,y,delta2,h)*y,backward(w,y,delta2,h)*y)<10**(-6)) then
	!	print*, 'delta', j
	!	exit
	!end if
	
	!if( abs(oldE - Err(f,r,m))< 10.0**(-8)) then
	!	print*, 'derr', j
	!	exit
	!end if
		

	
end do

print*, h, Err(f,r,m), Err(tf,tr,t) 

do
	print*, 'inserire matrice'

	x(n,:)=1
	do i=1,q
		read*, es(i,:)
	end do
	
	do i=1,q
		do j=1,q
			x(q*(i-1)+j,1)=es(i,j)
		end do
	end do
	
		yt(h)=1
	
	do k=1,h-1
		yt(k)=forward(x(:,1),p(:,k),n)
	end do
	ft(1)=forward(yt,w,h)
	
	print*, ft(1)
	

end do





end program