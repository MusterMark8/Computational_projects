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

Err=.5*dot_product((f-r),(f-r))

end function

end module fun



program conn
use fun
implicit none

!q nueroni hidden, n matrice 10x10 +1 (soglia), m esempi

integer, parameter :: q=6, n=37, m=5000, h=20, nstep=10, esp=6, con=4  !    
integer :: i, j, k, a, b, c, z, ai, aj
real :: ran, rand, dom, ness, u, v, delta1, delta2, e=0.01
real, dimension(q,q,m) :: es, tes
real, dimension(q/2,q/2,m/2,con) :: esd, tesd
real, dimension(n,m) :: x, tx
real, dimension(n,h-1) :: p
real, dimension(m) :: r, f, ft, tf, tr
real, dimension(h) :: w, y, yt, ty
real, dimension(n) :: t

!generare m esempi                                                         
	
	es=-1
	esd = -1
	tes = -1
	tesd = -1
	
do k=1,m/8							!esempi giusti
                                    !posizione casuale nella matrice
	call random_number(ran)
	i=int(q*ran)+1
	call random_number(rand)
	j=int(q*rand)+1
	es(i,j,k)=1
	
	call random_number(dom)
	do a=1, int(nstep*dom)+1
		call random_number(ness)
		b= int(4*ness)+1
			if (b==1 .and. j/=q) j=j+1
			if (b==2 .and. j/=1) j=j-1
			if (b==3 .and. i/=q) i=i+1
			if (b==4 .and. i/=1) i=i-1
			
		es(i,j,k)=1
	end do	
	
	r(k)=1

end do

do k=1,m/8							!4 componenti
                                    !posizione casuale nella matrice
                                    	
    do z=1,4
    
		call random_number(ran)
		i=int((q/2-1)*ran) +1
		call random_number(rand)
		j=int((q/2-1)*rand)+1
		esd(i,j,k,z)=1
	
		call random_number(dom)
		do a=1, int(nstep/4*dom)+1
			call random_number(ness)
			b= int(4*ness)+1
			
				if (b==1 .and. j/=q/2) j=j+1
				if (b==2 .and. j/=1) j=j-1
				if (b==3 .and. i/=q/2) i=i+1
				if (b==4 .and. i/=1) i=i-1
			
			esd(i,j,k,z)=1
		end do	
	
	end do
	
	
	do z=1,4
			if (z==1) then
				ai=0
				aj=0
				b=1
			end if
			if (z==2) then
				ai=0
				aj=1
				b=2
			end if
			if (z==3) then
				ai=1
				aj=0
				b=2
			end if
			if (z==4) then
				ai=1
				aj=1
				b=1
			end if
	do i=1,q/2
		do j=1,q/2
			es(i,j,m/8*z+k)=esd(i,j,k,z)
			es(i+ai*q/2,j+aj*q/2,5*m/8+k)=esd(i,j,k,z)
			es(i+ai*q/2,j+aj*q/2,5*m/8+m/8*b+k)=esd(i,j,k,z)
			
		end do
	end do
	
				r(m/8*z+k)=1
				r(5*m/8+k)=-1
				r(5*m/8+m/8*b+k)=-1
	
	end do
	

end do

print*, sum(r)

!generare campione test

	
do k=1,m/2							!esempi giusti
                                    !posizione casuale nella matrice
	call random_number(ran)
	i=int(q*ran)+1
	call random_number(rand)
	j=int(q*rand)+1
	tes(i,j,k)=1
	
	call random_number(dom)
	do a=1, int(nstep*dom)+1
		call random_number(ness)
		b= int(4*ness)+1
			if (b==1 .and. j/=q) j=j+1
			if (b==2 .and. j/=1) j=j-1
			if (b==3 .and. i/=q) i=i+1
			if (b==4 .and. i/=1) i=i-1
			
		tes(i,j,k)=1
	end do
	
	tr(k)=1	

end do

do k=1,m/2							!4 componenti
                                    !posizione casuale nella matrice
                                    	
    do z=1,4
    
		call random_number(ran)
		i=int((q/2-1)*ran) +1
		call random_number(rand)
		j=int((q/2-1)*rand)+1
		tesd(i,j,k,z)=1
	
		call random_number(dom)
		do a=1, int(nstep/2*dom)+4
			call random_number(ness)
			b= int(4*ness)+1
			
				if (b==1 .and. j/=q/2) j=j+1
				if (b==2 .and. j/=1) j=j-1
				if (b==3 .and. i/=q/2) i=i+1
				if (b==4 .and. i/=1) i=i-1
			
			tesd(i,j,k,z)=1
		end do	
	
	end do
	
	tes(:,:,m/2+k)=0
	
	do z=1,4
			if (z==1) then
				ai=0
				aj=0
			end if
			if (z==2) then
				ai=0
				aj=1
			end if
			if (z==3) then
				ai=1
				aj=0
			end if
			if (z==4) then
				ai=1
				aj=1
			end if
	do i=1,q/2
		do j=1,q/2
			tes(i+ai*q/2,j+aj*q/2,m/2+k)=tesd(i,j,k,z)
		end do
	end do
	
	end do
	
	tr(m/2+k)=-1	

end do

print*, 'fatto'
	
!check

do i=1,q
	print*, es(i,:,3)
end do

do i=1,q
	print*, es(i,:,3*m/4 + 3)
end do

!coversione matrice array

x(n,:)=1									!esempi
do i=1,q
	do j=1,q
		x(q*(i-1)+j,:)=es(i,j,:)
	end do
end do

tx(n,:)=1									!test
do i=1,q
	do j=1,q
		tx(q*(i-1)+j,:)=tes(i,j,:)
	end do
end do


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
	
	

	if( mod(j,int(10**(esp-3)))==0 ) then
		
		do i=1,m
		
			do k=1,h-1  
				y(k)=forward(x(:,i),p(:,k),n)
				ty(k)=forward(tx(:,i),p(:,k),n)
			end do
			f(i)=forward(y,w,h)
			tf(i)=forward(y,w,h)
						
		end do
		
		write(1,*) j, Err(f,r,m)
		write(2,*) j, Err(tf,tr,m)                       !attenzione all'ordine
	end if	
	
	if( Err(f,r,m)<0.1 ) then
		print*, 'err', j
		exit
	end if
	
	if(e**2*dot_product(backward(w,y,delta2,h)*y,backward(w,y,delta2,h)*y)<10**(-8)) then
		print*, 'delta', j
		exit
	end if

	
end do

print*, Err(f,r,m)            

!TEST
!		es=-1
	
!do k=1,m
!                                    !posizione casuale nella matrice
!	call random_number(ran)
!	i=int(q*ran)+1
!	call random_number(rand)
!	j=int(q*rand)+1
!	es(i,j,k)=1
!	
!	call random_number(dom)
!	do a=1, int(nstep*dom)+4
!		call random_number(ness)
!		b= int(4*ness)+1
!			if (b==1 .and. j/=q) j=j+1
!			if (b==2 .and. j/=1) j=j-1
!			if (b==3 .and. i/=q) i=i+1
!			if (b==4 .and. i/=1) i=i-1
!			
!		es(i,j,k)=1
!	end do	
!
!end do

print*, 'fatto'


!x(n,:)=1
!do i=1,q
!	do j=1,q
!		x(q*(i-1)+j,:)=es(i,j,:)
!	end do
!end do

!do i=1,m	
!	do k=1,h
!		yt(k)=forward(x(i,:),p(:,k),n)
!	end do
!	ft(i)=forward(yt,w,h)
!end do

!print*, Err(ft,r,m)

!test manuale

do
	print*, 'inserire matrice'

	x(n,:)=1
	do i=1,q
		read*, es(i,:,1)
	end do
	
	do i=1,q
		do j=1,q
			x(q*(i-1)+j,1)=es(i,j,1)
		end do
	end do
	
		yt(h)=1
	
	do k=1,h-1
		yt(k)=forward(x(:,1),p(:,k),n)
	end do
	ft(1)=forward(yt,w,h)
	
	print*, ft(1)
	

end do

end program conn