	program test
	implicit none
        
        integer r
        parameter (r=125)
        integer ms
        parameter (ms=12)
        integer nj
        parameter (nj=ms*ms*ms)
        integer i,iflag,xsub(nj,3),ier,num,j,xxsub(nj)
        integer n1,n2,mt,k1,k2,mm,k3
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(r,ms*ms*ms),V1(r,nj),U2(r,ms*ms*ms),V2(r,nj)
        real*16 re1(nj),re2(nj),time1,time2
        real*16 arr(4)
        real*8 pi,xj(nj),yj(nj),zj(nj)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(ms*ms*ms,r),V(r,nj),cj(nj),S(ms*ms*ms),re(nj)
        complex*16 fk(ms,ms,ms),fk1(ms*ms*ms),uu(ms,ms,ms,r)
        complex*16,allocatable :: NN(:,:,:,:)
        real*8 x(nj,3),eps,error
        double complex in1, out1
        dimension in1(ms,ms,ms), out1(ms,ms,ms)
	integer*16 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001


        iflag=-1
        eps=1E-12
        num=100
        
        open(unit = 10,file = 'U3r1.txt')
        read(10,*) U1
        open(unit = 20,file = 'V3r1.txt')
        read(20,*) V1
        open(unit = 10,file = 'U3i1.txt')
        read(10,*) U2
        open(unit = 10,file = 'V3i1.txt')
        read(10,*) V2
        !open(unit = 10,file = 'Re2r1.txt')
        !read(10,*) re1
        !open(unit = 10,file = 'Re2i1.txt')
        !read(10,*) re2
        !re=dcmplx(re1,re2)


        call dfftw_plan_dft_3d(plan,ms,ms,ms,in1,out1,FFTW_FORWARD,0)

        print *,'start 3D type 1 testing:','nj  =',nj,'ms  =',ms
        print *,'eps             =',eps
        !re=dcmplx(re1,re2)        
        U=transpose(dcmplx(U1,U2))
        
        V=dcmplx(V1,V2)
        V=conjg(V)
        !print *,V(1:5,1)
        uu=reshape(U,(/ms,ms,ms,r/))
        mm=floor(ms/2.0+0.6)
        allocate(NN(mm,mm,mm,r))
        NN=uu(1:mm,1:mm,1:mm,:)
        uu(1:mm,1:mm,1:mm,:)=uu(mm+1:ms,mm+1:ms,mm+1:ms,:)
        uu(mm+1:ms,mm+1:ms,mm+1:ms,:)=NN
        NN=uu(1:mm,mm+1:ms,1:mm,:)
        uu(1:mm,mm+1:ms,1:mm,:)=uu(mm+1:ms,1:mm,mm+1:ms,:)
        uu(mm+1:ms,1:mm,mm+1:ms,:)=NN
        NN=uu(mm+1:ms,1:mm,1:mm,:)
        uu(mm+1:ms,1:mm,1:mm,:)=uu(1:mm,mm+1:ms,mm+1:ms,:)
        uu(1:mm,mm+1:ms,mm+1:ms,:)=NN
        NN=uu(mm+1:ms,mm+1:ms,1:mm,:)
        uu(mm+1:ms,mm+1:ms,1:mm,:)=uu(1:mm,1:mm,mm+1:ms,:)
        uu(1:mm,1:mm,mm+1:ms,:)=NN
        U=reshape(uu,(/ms*ms*ms,r/))
        !U=conjg(U)
        !print *,V(2,:)
        !print *,U(1,:)
        do k3 = -ms/2,(ms-1)/2
	 do k1 = -ms/2, (ms-1)/2
	   do k2 = -ms/2, (ms-1)/2
	      j =  (k1+ms/2+1) + (k2+ms/2)*ms + (k3+ms/2)*ms*ms
	      xj(j) = pi*dcos(-pi*k1/ms)
	      yj(j) = pi*dcos(-pi*k2/ms)
	      zj(j) = pi*dcos(-pi*k3/ms)
	      cj(j) = dcmplx(dsin(pi*j/ms),dcos(pi*j/ms))
	   enddo
	 enddo
        enddo
        x(:,1)=xj
        x(:,2)=yj
        x(:,3)=zj
        xsub=mod(floor(x/2/pi*ms+0.5),ms)+1
        do i = 1,nj
           xxsub(i)=xsub(i,3)*ms*ms-ms*ms+xsub(i,2)*ms-ms+xsub(i,1)
        enddo
        !print *,'xxsub=',size(xxsub),xxsub
        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft3dIapp(nj,plan,cj,U,V,xxsub,ms,-1,r,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(1:5)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num
        
        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft3d1f90(nj,xj,yj,zj,cj,iflag,eps,ms,ms,ms,fk,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
        fk1=reshape(fk,(/ms*ms*ms/))
        !print *,'S/nj(1:5)',S(1:5)/nj
        !print *,'fk1(1:5)',fk1(1:10)
        error=sqrt(real(sum((S/nj-fk1)*conjg(S/nj-fk1))/
     &  sum(fk1*conjg(fk1))))
        print *,' relative error= ',error
        
        call dfftw_destroy_plan(plan)
        

	end program
