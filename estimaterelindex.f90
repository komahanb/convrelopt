program main
  
  implicit none

  integer n,nmax,d
  double precision tolerance,x(2),grad(2),hess(2,2),y(2),Ihess(2,2)
  double precision :: u(2)
  double precision :: beta,betaold,eps
  double precision :: g

  double precision :: sigmax(2),meanx(2)
  double precision :: val
  
  n=0 !iteration counter

!  tolerance=1.0D320	

  eps=1.0d320

  nmax=1000	              	!MAXIMUM ITERATIONS

  x=(/-0.8,0.8/)		!INITIAL GUESS

  sigmax(1)=0.1
  sigmax(2)=0.1

  meanx(1)=-0.8
  meanx(2)=0.8
 

  do while ((eps.gt.1.0D-3).and.(n <= 2))

!!$
!!$     if (n.ge.nmax) then
!!$        if (residual>eps) then 
!!$           print *,'Not converged after max iterations'
!!$           stop
!!$        end if
!!$     end if
     
!     call transform(2,x,meanx,sigmax,u)

     call get_G(2,x,g)    
     print*,G,x,u

     call get_dG(2,x,grad)
     print*,grad
!     call transform(2,grad,meanx,sigmax,grad)
!     print*,grad
     call transform(2,x,meanx,sigmax,x)
     print*,x
     grad(1)=grad(1)*sigmaX(1)
     grad(2)=grad(2)*sigmaX(2)

!     print*,g
  !   print*,grad
!stop

 !    x(1)=(x(1)-(-0.8))/sigmax(1)
 !    x(2)=(x(2)-(0.8))/sigmax(2)
!     print*,x
 !    stop

     val=(DOT_PRODUCT(grad,x)-g)/DOT_PRODUCT(grad,grad)

     y=val*grad

     beta=dsqrt(DOT_PRODUCT(y,y))

     print*,beta
!     stop
     eps=(beta-betaold)/beta
     betaold=beta

     print*,y
!stop
     call invtransform(2,y,meanx,sigmax,x)
     print*,x

!     x=y
     n=n+1

  end do

  print*,"beta found:",beta,n+1

!!$  WRITE(*,*)''
!!$  WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!$  WRITE(*,*)''




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


!!$     WRITE(*,'(i4,3f15.4,3f15.4,3f15.4,3f15.4)') n,x(1),x(2),x(3),residual

!     grad=(/18*x(1)+18*x(2)+x(3)*2*x(1)+2*x(3), 18*x(1)+26*x(2)+2*x(3)*x(2), x(1)**2+x(2)**2+2*x(1)-16/)

!!$     hess(1,1)=18+2*x(3)
!!$     hess(1,2)=18
!!$     hess(1,3)=2*x(1)+2
!!$
!!$     hess(2,1)=18
!!$     hess(2,2)=26+2*x(3)
!!$     hess(2,3)=2*x(2)
!!$
!!$     hess(3,1)=2*x(1)+2
!!$     hess(3,2)=2*x(2)
!!$     hess(3,3)=0


!!$     CALL FINDInv(hess,Ihess,3,errorflag)

     
 !    y=x-matmul(ihess,grad)
  !   x=y

   !  residual=norm2(grad) 


!!$
!!$  WRITE(*,*)''
!!$  WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!$  WRITE(*,*)'%%%%              CLASSICAL NEWTONS METHOD                        %%%%'
!!$  WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!$  WRITE(*,*)''
!!$  WRITE(*,*)'Iteration', '      X1','            X2','              X3','         Residual'
!!$  WRITE(*,*)''
!!$


