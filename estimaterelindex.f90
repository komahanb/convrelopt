program main
  
  implicit none

  integer n,nmax,d
  double precision tolerance,x(2),grad(2),hess(2,2),Ihess(2,2)
  double precision :: u(2)
  double precision :: beta,betaold,eps
  double precision :: g

  double precision :: sigmax(2),meanx(2)
  double precision :: val
  integer::nvar,i


  ! User input

  x=(/-0.8,0.8/) ! Initial point
  nvar=2
  sigmax(1)=0.1
  sigmax(2)=0.1

  meanx(1)=-0.8
  meanx(2)=0.8

  ! Stopping parameters
  
  tolerance=1.0D-10
  nmax=1000 !set max iterations

  !Initialize
  
  n=0 !iteration counter
  eps=1.0d320
  beta=0.0
  betaold=0.0

  do while ((eps.gt.tolerance).and.(n <= nmax))

     call get_G(nvar,x,g)    !get G value in x space

     call get_dG(nvar,x,grad) ! get grad in x space
     do i=1,nvar
        grad(i)=grad(i)*sigmaX(i) ! convert grads to u space
     end do

     call transform(nvar,x,meanx,sigmax,x) ! map from x to u
     val=(DOT_PRODUCT(grad,x)-g)/DOT_PRODUCT(grad,grad) ! in u space
     
     do i=1,nvar
        u(i)=val*grad(i) ! in u space
     end do

     beta=dsqrt(DOT_PRODUCT(u,u)) ! in u space

     eps=abs((beta-betaold)/beta)
     betaold=beta

     call invtransform(nvar,u,meanx,sigmax,x) ! map from u to x
     n=n+1
  end do
  print*,"eps:",eps
  print*,"n :",n
  print*,"b*:",beta
  print*,"x*:",x
  print*,"u*:",u
  
end program main
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

