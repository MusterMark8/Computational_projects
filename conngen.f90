module funz
implicit none
contains

function search(x,n,a,b)
integer, intent(in) :: a, b, n
real, dimension(n,n) :: x
real :: search
integer :: i,j

search = 0

if(a==n .or. b==n) then
	
	if(a==n .and. b==n) then
		search = search - 5
	else
		search = search - 3
	end if

end if

do i=a-1,a+1
	do j=b-1,b+1
		search = search + x(i,j)
	end do
end do

end function

end module funz


program imm
use funz
implicit none

integer, parameter :: q=7, n=50, m=16000, t=4000, nstep=20
integer :: i, j, k, a, b
real :: ran, rand, dom, ness
real, dimension(q,q,m) :: es, copy
real, dimension(q,q,t) :: tes, tcopy
real, dimension(m) :: r
real, dimension(t) :: tr

es=-1
copy=-1

tes=-1
tcopy=-1

!generare m esempi                                                         

do k=1,m/2							!esempi giusti
                                    !posizione casuale nella matrice
	call random_number(ran)
	i=int(q*ran)+1
	call random_number(rand)
	j=int(q*rand)+1
	es(i,j,k)=1
	
	!print*, search(es(i,j,k)
	
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

print*, 'connesse'

do k=1,m/2

	call random_number(ran)
	i=int(q*ran)+1
	call random_number(rand)
	j=int(q*rand)+1
	es(i,j,k+m/2)=1
	
	call random_number(dom)
	do a=1, int(nstep/2*dom)+1
		call random_number(ness)
		b= int(4*ness)+1
			if (b==1 .and. j/=q) j=j+1
			if (b==2 .and. j/=1) j=j-1
			if (b==3 .and. i/=q) i=i+1
			if (b==4 .and. i/=1) i=i-1
			
		es(i,j,k+m/2)=1
	end do
	
	copy = es
	
	do
	
		call random_number(ran)
		i=int(q*ran)+1
		call random_number(rand)
		j=int(q*rand)+1
		
		if(search(es(:,:,k+m/2),q,i,j)==-9) then
	
		es(i,j,k+m/2)=1
		call random_number(dom)
		do a=1, int(nstep/2*dom)+1
		
			call random_number(ness)
			b= int(4*ness)+1
			if (b==1 .and. j/=q .and. search(copy(:,:,k+m/2),q,i,j+1)==-9) j=j+1
			if (b==2 .and. j/=1 .and. search(copy(:,:,k+m/2),q,i,j-1)==-9) j=j-1
			if (b==3 .and. i/=q .and. search(copy(:,:,k+m/2),q,i+1,j)==-9) i=i+1
			if (b==4 .and. i/=1 .and. search(copy(:,:,k+m/2),q,i-1,j)==-9) i=i-1
			
			es(i,j,k+m/2)=1
			
		end do
		
		exit
		end if
		
	end do
	
	r(k+m/2)=-1
	
end do

do i=1,q
	print*, es(i,:,m/2+1)
end do

!generare test
	
do k=1,t/2							!esempi giusti
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

do k=1,t/2

	call random_number(ran)
	i=int(q*ran)+1
	call random_number(rand)
	j=int(q*rand)+1
	tes(i,j,k+t/2)=1
	
	call random_number(dom)
	do a=1, int(nstep/2*dom)+1
		call random_number(ness)
		b= int(4*ness)+1
			if (b==1 .and. j/=q) j=j+1
			if (b==2 .and. j/=1) j=j-1
			if (b==3 .and. i/=q) i=i+1
			if (b==4 .and. i/=1) i=i-1
			
		tes(i,j,k+t/2)=1
	end do
	
	tcopy = tes
	
	do
	
		call random_number(ran)
		i=int(q*ran)+1
		call random_number(rand)
		j=int(q*rand)+1
		
		if(search(tes(:,:,k+t/2),q,i,j)==-9) then
	
		es(i,j,k+m/2)=1
		call random_number(dom)
		do a=1, int(nstep/2*dom)+1
		
			call random_number(ness)
			b= int(4*ness)+1
			if (b==1 .and. j/=q .and. search(tcopy(:,:,k+t/2),q,i,j+1)==-9) j=j+1
			if (b==2 .and. j/=1 .and. search(tcopy(:,:,k+t/2),q,i,j-1)==-9) j=j-1
			if (b==3 .and. i/=q .and. search(tcopy(:,:,k+t/2),q,i+1,j)==-9) i=i+1
			if (b==4 .and. i/=1 .and. search(tcopy(:,:,k+t/2),q,i-1,j)==-9) i=i-1
			
			tes(i,j,k+t/2)=1
			
		end do
		
		exit
		end if
		
	end do
	
	tr(k+t/2)=-1
	
end do

!salva

	do i=1,q
		do j=1,q
			write(13,*) es(i,j,:)
		end do
	end do
	
	do i=1,q
		do j=1,q
			write(14,*) tes(i,j,:)
		end do
	end do	

	do i=1,q
		do j=1,q
			do k=1,m/2
				write(15,*) es(i,j,k)
			end do
		end do
	end do
	
	do i=1,q
		do j=1,q
			do k=m/2 +1, m
				write(16,*) es(i,j,k)
			end do
		end do
	end do

end program imm