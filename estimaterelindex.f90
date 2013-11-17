program main
  
  implicit none

  integer n,nmax,d
  double precision tolerance,x(2),grad(2),hess(2,2),Ihess(2,2)
  double precision :: u(2)
  double precision :: beta,betaold,eps
  double precision :: g

  double precision :: sigmax(2),meanx(2)
  double precision :: val
  double precision::pfail
  integer::nvar,i


  ! User input

  x=(/-1.8,1.8/) ! Initial point
  nvar=2
  sigmax(1)=0.1
  sigmax(2)=0.1

  meanx(1)=-0.8
  meanx(2)=0.8

  ! Stopping parameters
  
  tolerance=1.0D-15
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

  print*,"g:",g
  print*,"grad:",grad
  
  print*,"eps:",eps
  print*,"n :",n
  print*,"b*:",beta
  print*,"x*:",x
  print*,"u*:",u

  call CDF(-beta,0.d0,1.d0,pfail)

  print*,"pfail:",pfail
  
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
subroutine CDF(xin,xc,st,vout)
  implicit none
  double precision, intent(in)  :: xin,xc,st
  double precision, intent(out) :: vout
  double precision :: vtmp
  !       vout = 0.5d0 * (1.d0 + erf( (xin-xc)/(st*dsqrt(2.d0)) ))
  call ERF_MINE1( (xin-xc)/(st*dsqrt(2.d0)), vtmp )
  vout = 0.5d0 * (1.d0 + vtmp)
end subroutine CDF
!++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ERF_MINE1(xin,yout)
  implicit none
  double precision, intent(in)  :: xin
  double precision, intent(out) :: yout
  integer :: i,k,n
  double precision :: vsum,kai
  ! n is the order of Taylor
  ! Maybe accurate within the range of [-4:4] with n=100
  n = 100
  vsum = 0.d0
  do i=0,n
     kai = 1.d0
     do k=1,i
        kai = kai * (-1.d0) * xin**2 / dble(k)
     end do
     vsum = vsum + kai*xin/dble(2*i+1)
  end do
  yout = vsum*2.d0/(dsqrt(3.141592653589793238d0))

  if(yout.gt.1.d0)write(*,'(a,2e15.5)')'*ERF>1 ',xin,yout-1.d0
end subroutine ERF_MINE1
!++++++++++++++++++++++++++++++++++++++++++++++++
