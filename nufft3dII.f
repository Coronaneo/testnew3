c	subroutine nufft2dII(nj,x,iflag,ns,rt,tol,U,V,xsub,r)
c	implicit none
c
c	integer :: ns,rt,iflag,nj,k(ns*ns,2),j,r,xsub(nj),i,xxsub(nj)
c	real  :: tol
c	real*8 pi,x(nj)
c	parameter (pi=3.141592653589793238462643383279502884197d0)
c	complex*16 fftconst,U(nj,r),V(ns*ns,r)
c	complex*16 M(nj,ns*ns)
c
c	do j = 1,ns
c	   do i = 1,ns
c	      k((j-1)*ns+i,1) = i
c	      k((j-1)*ns+i,1) = j
c	   enddo
c	enddo
c
c	fftconst = iflag*dcmplx(0,1)/ns*2*pi
c
c	do i = 1,nj
c	   do j = 1,ns*ns
c	      M(i,j) = exp(fftconst*(k(j,1)*(x(i,1)-floor(x(i,1)+0.5))+
c	&     k(j,2)*(x(i,2)-floor(x(i,2)+0.5))))
c	   enddo
c	enddo
c
c	call lowrankfac(M,tol,rt,rt,U,V)
c
c	xsub = mod(floor(x+0.5),ns)+1
c	do i = 1,nj
c	   xxsub(i) = xsub(i,2)*ns-ns+xsub(i,1)
c	enddo
c
c	r = size(V,2)
c
c	end subroutine

	subroutine nufft3dIIapp(nj,plan,c,U,V,xxsub,ns,iflag,r,S)
	implicit none
	integer  r,i,j,k,nj,ns,iflag,num
	integer mm
	integer xxsub(nj)
	complex*16 M(r,ns*ns*ns),N(ns*ns*ns,r),S(nj),c(ns*ns*ns),U(nj,r)
        complex*16 cc(ns,ns,ns),c2(ns*ns*ns),V(r,ns*ns*ns)
        complex*16,allocatable :: c1(:,:,:)
	double complex in1, out1
	real*16  time_begin,time_end,countrage,countmax
	dimension in1(ns,ns,ns), out1(ns,ns,ns)
	integer*16 :: plan
	integer FFTW_FORWARD,FFTW_MEASURE
	parameter (FFTW_FORWARD=-1)
	parameter (FFTW_MEASURE=0)

	
        mm=floor(ns/2.0+0.6)
        cc=reshape(c,(/ns,ns,ns/))
        allocate(c1(mm,mm,mm))
        c1=cc(1:mm,1:mm,1:mm)
        cc(1:mm,1:mm,1:mm)=cc(mm+1:ns,mm+1:ns,mm+1:ns)
        cc(mm+1:ns,mm+1:ns,mm+1:ns)=c1
        c1=cc(1:mm,mm+1:ns,1:mm)
        cc(1:mm,mm+1:ns,1:mm)=cc(mm+1:ns,1:mm,mm+1:ns)
        cc(mm+1:ns,1:mm,mm+1:ns)=c1
        c1=cc(mm+1:ns,1:mm,1:mm)
        cc(mm+1:ns,1:mm,1:mm)=cc(1:mm,mm+1:ns,mm+1:ns)
        cc(1:mm,mm+1:ns,mm+1:ns)=c1
        c1=cc(mm+1:ns,mm+1:ns,1:mm)
        cc(mm+1:ns,mm+1:ns,1:mm)=cc(1:mm,1:mm,mm+1:ns)
        cc(1:mm,1:mm,mm+1:ns)=c1
        c2=reshape(cc,(/ns*ns*ns/))
        M=0

	do i = 1,ns*ns*ns
	   do k = 1,r
	      M(k,i) = V(k,i)*c2(i)
              
	   enddo
	enddo
        
	do i = 1,r
	   in1 = reshape(M(i,:),(/ns,ns,ns/))
	   call dfftw_execute_dft(plan, in1, out1)
	   N(:,i) = reshape(out1,(/ns*ns*ns/))
	enddo
	
        
	S = sum(U*N(xxsub,:),2)

	end subroutine
