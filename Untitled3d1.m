format long
ms = 12;
mt = 12;
mu = 12;
n1 = 12;
n2 = 12;
n3 = 12;
nj = n1*n2*n3;
iflag=-1;
eps=1e-12;
xj=zeros(nj,1);
yj=zeros(nj,1);
zj=zeros(nj,1);
cj=zeros(nj,1);
for k3 = -n3/2:(n3-1)/2
 for k2 = -n2/2:(n2-1)/2
    for k1 = -n1/2:(n1-1)/2
       j =  (k1+n1/2+1) + (k2+n2/2)*n1 + (k3+n3/2)*n1*n2;
       xj(j) = pi*cos(-pi*k1/n1);
       yj(j) = pi*cos(-pi*k2/n2);
       zj(j) = pi*cos(-pi*k3/n3);
       cj(j) = sin(pi*j/n1)+1i*cos(pi*j/n2);
    end
 end
end
%xj=pi*rand(nj,1);
%yj=pi*rand(nj,1);
%cj=pi*rand(nj,1);
x=[xj(:) yj(:) zj(:)];
[k1,k2,k3] = ndgrid(-ms/2:ms/2-1,-ms/2:ms/2-1,-ms/2:ms/2-1);
k = [k1(:) k2(:) k3(:)];


fftconst = iflag*1i/ms*2*pi;
ratiofun = @(k,x)exp(fftconst*k*(x-round(x))');
[U,V] = lowrank(k,x/2/pi*ms,ratiofun,eps,500,500);

xsub = mod(round(x/2/pi*ms),ms)+1;
xxsub = sub2ind([ms ms ms],xsub(:,1),xsub(:,2),xsub(:,3));
spPerm = sparse(xxsub,1:ms^3,ones(1,ms^3),ms^3,ms^3);
r = size(V,2)
[n,ncol] = size(cj);

M = repmat(conj(V),1,ncol).*reshape(repmat(cj,r,1),n,r*ncol);
MM = reshape(spPerm*M,ms,ms,ms,r*ncol);

MMM = fft3(MM);
MMMM=fftshift(fftshift(fftshift(MMM,1),2),3);


fhat = sum(U.*reshape(MMMM,ms*ms*ms,r),2);

fhat = fhat/nj;




fid=fopen('U3r1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',U.');
fclose(fid);
fid=fopen('V3r1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',V.');
fclose(fid);
fid=fopen('U3i1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*U.');
fclose(fid);
fid=fopen('V3i1.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',-1i*V.');
fclose(fid);

