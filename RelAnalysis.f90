program main
  implicit none
  integer,parameter::nvar=2
  double precision :: sigmax(nvar),meanx(nvar) ! input mean and variance
  double precision :: beta,x(nvar)

  meanx(1)=-0.8
  meanx(2)=0.8

  sigmax(:)=0.1
  
  call relanalysis(nvar,meanx,sigmax,beta,x)

  print*,beta
  print*,x


end program main



subroutine relanalysis(nvar,meanx,sigmax,beta,x)
  implicit none

  integer,intent(in)::nvar

  double precision,intent(out) :: x(nvar) ! design vector
  double precision,intent(in)  :: sigmax(nvar),meanx(nvar) ! input mean and variance
  integer:: n,nmax,i ! iterations
  double precision tolerance,eps ! tolerance and stopping
  double precision :: u(nvar)
  double precision :: beta,betaold ! MPP
  double precision :: g,grad(nvar) ! limit state function
  double precision :: tempval
  double precision::pfail

  ! User input

!  x=(/-1.8,1.8/) ! Initial point
!  nvar=2

 ! sigmax(1)=0.1
 ! sigmax(2)=0.1

  do i =1,nvar
     x(i)=meanx(i)
  end do

  ! Stopping parameters
  
  tolerance=1.0D-15
  nmax=100 !set max iterations

  !Initialize
  
  n=0 !iteration counter
  eps=1.0d320
  beta=0.0
  betaold=0.0

  do while ((eps.gt.tolerance).and.(n <= nmax))

     if (n.eq.nmax) then
       write(*,'(a,i8,a)')"Finding MPP failed after",n,"   iterations"
     end if


     call get_G(nvar,x,g)    !get G value in x space
     call get_dG(nvar,x,grad) ! get grad in x space

     do i=1,nvar
        grad(i)=grad(i)*sigmaX(i) ! convert grads to u space
     end do

     call transform(nvar,x,meanx,sigmax,x) ! map from x to u
     tempval=(DOT_PRODUCT(grad,x)-g)/DOT_PRODUCT(grad,grad) ! in u space
     
     do i=1,nvar
        u(i)=tempval*grad(i) ! in u space
     end do

     beta=dsqrt(DOT_PRODUCT(u,u)) ! in u space

     eps=abs((beta-betaold)/beta)
     betaold=beta

     call invtransform(nvar,u,meanx,sigmax,x) ! map from u to x
     n=n+1
  end do

  write(*,*),">> Beta found successfully"

!!$
!!$  print*,"g:",g
!!$  print*,"grad:",grad
!!$  print*,"eps:",eps
!!$  print*,"n :",n
!!$  print*,"b*:",beta
!!$  print*,"x*:",x
!!$  print*,"u*:",u
!!$

!!$  ! Gradient in X-Space
!!$  do i=1,nvar
!!$     grad(i)=grad(i)/sigmaX(i) ! convert grads to u space
!!$  end do
!!$  !call get_dG(nvar,x,grad) ! get grad in x space
!!$
!!$  call CDF(-beta,0.d0,1.d0,pfail)
!!$
!!$  print*,"pfail:",pfail

  return
end subroutine relanalysis
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_G(nvar,x,g)
  implicit none

  integer::nvar

  double precision,intent(in)  :: x(nvar)
  double precision, intent(out):: g

  g= -(x(1)-1.0)**2 - x(2)**2 + x(1) + 6

  return
end subroutine get_G
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_dG(nvar,x,dg)
  implicit none

  integer::nvar

  double precision,intent(in)  :: x(nvar)
  double precision, intent(out):: dG(nvar)

  dG(1)= -2.0*x(1) + 3.0
  dG(2)= -2.0*x(2)

  return
end subroutine get_dG
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine transform(nvar,x,meanx,sigmax,u)
  implicit none
  integer,intent(in):: nvar
  double precision,intent(in) ::x(nvar),sigmax(nvar),meanx(nvar)
  double precision,intent(out)::u(nvar)
  integer::i

  do i=1,nvar
     u(i)=(x(i)-meanx(i))/sigmax(i)
  end do

  return
end subroutine transform
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine invtransform(nvar,u,meanx,sigmax,x)
  implicit none
  integer,intent(in):: nvar
  double precision,intent(in) ::u(nvar),sigmax(nvar),meanx(nvar)
  double precision,intent(out)::x(nvar)
  integer::i

  do i=1,nvar
     x(i)=u(i)*sigmax(i)+meanx(i)
  end do

  return
end subroutine invtransform
!++++++++++++++++++++++++++++++++++++++++++++++++
