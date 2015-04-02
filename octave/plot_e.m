clear;

ExGSE=load('ExGSEk.txt');
EyGSE=load('EyGSEk.txt');
EzGSE=load('EzGSEk.txt');
E=load("E_LES.txt");

ExLES=load('ExLES.txt');
EyLES=load('EyLES.txt');
EzLES=load('EzLES.txt');

L=20;
nx=20;
ny=20;
nz=20;
ccelec=332.0716;

xm = zeros(nx*ny,1);
ym = zeros(nx*ny,1);
i=1;
for y=0:ny-1
    for x=0:nx-1
        xm(i) = x+0.5;
        ym(i) = y+0.5;
        i = i + 1;
    end
end

% Build 3d interpolation
phi=load('phi.txt');
m=0:nx-1;
ExV=zeros(nx,ny,nz);
EyV=zeros(nx,ny,nz);
EzV=zeros(nx,ny,nz);
ExV2=zeros(nx,ny,nz);
EyV2=zeros(nx,ny,nz);
EzV2=zeros(nx,ny,nz);
phiV=zeros(nx,ny,nz);
i=1;
for z=1:nz
    for x=1:nx
        for y=1:ny
            ExV(x,y,z) = ExLES(i);
            EyV(x,y,z) = EyLES(i);
            EzV(x,y,z) = EzLES(i);
            ExV2(x,y,z) = ExGSE(i);
            EyV2(x,y,z) = EyGSE(i);
            EzV2(x,y,z) = EzGSE(i);
            phiV(x,y,z) = phi(i);
            i = i + 1;
        end
    end
end

%quiver3(E(1:368,1),E(1:368,2),E(1:368,3),E(1:368,4),E(1:368,5),E(1:368,6));

ind = find(E(1:368,3)==10);
figure(1);
%quiver(E(ind,1),E(ind,2),E(ind,4),E(ind,5));
plot(9.5,10.5,'*k');
hold on;
plot(10.5,10.5,'*k');
f=5;
quiver(xm,ym,ExLES(11*nx*ny+1:12*nx*ny)*f,EyLES(11*nx*ny+1:12*nx*ny)*f,0,'b');
quiver(xm,ym,ExGSE(11*nx*ny+1:12*nx*ny)*f,EyGSE(11*nx*ny+1:12*nx*ny)*f,0,'r');
hold off;

xlabel('x');
ylabel('y');
axis([6 14 7 14]);

%find(sum((E(1:368,1:3)-[10 10 10]).^2,2)==0)

EdE = zeros(368,1);
Exdip = 4*pi*1/(L*L*L);
val = 0;
for i=1:368
    x = E(i,1);
    y = E(i,2);
    z = E(i,3);
    %Ex_m1 = interp3(m,m,m,ExV,x-1,y,z);
    %Ex_nn = interp3(m,m,m,ExV,x,y,z);
    %Ex_p1 = interp3(m,m,m,ExV,x+1,y,z);
    %dEx_m1= interp3(m,m,m,ExV,x,y,z) - interp3(m,m,m,ExV,x-1,y,z);
    %dEx_nn= interp3(m,m,m,ExV,x+1,y,z) - interp3(m,m,m,ExV,x,y,z);
    %dEx_p1= interp3(m,m,m,ExV,x+2,y,z) - interp3(m,m,m,ExV,x+1,y,z);
    %EdE(i)=((20/12)*Ex_nn*dEx_nn + (1/12)*(Ex_nn*dEx_p1 + Ex_p1*dEx_nn ...
    %                                       + Ex_nn*dEx_m1 +
    %                                       Ex_m1*dEx_nn));
    
    Ex_m2 = interp3(m,m,m,ExV,x-2,y,z);
    Ex_m1 = interp3(m,m,m,ExV,x-1,y,z);
    Ex_nn = interp3(m,m,m,ExV,x,y,z);
    Ex_p1 = interp3(m,m,m,ExV,x+1,y,z);
    Ex_p2 = interp3(m,m,m,ExV,x+2,y,z);
    Ex_p3 = interp3(m,m,m,ExV,x+3,y,z);
    %dEx_m1 = (10.0*(Ex_nn-Ex_m1) + (Ex_m1-Ex_m2) + (Ex_p1-Ex_nn))/12.0;
    dEx_nn = (10.0*(Ex_p1-Ex_nn) + (Ex_nn-Ex_m1) + (Ex_p2-Ex_p1))/12.0;
    %dEx_p1 = (10.0*(Ex_p2-Ex_p1) + (Ex_p1-Ex_nn) + (Ex_p3-Ex_p2))/12.0;
    %dEx_m1 = (Ex_nn-Ex_m1);
    %dEx_nn = (Ex_p1-Ex_nn);
    %dEx_p1 = (Ex_p2-Ex_p1);
    %EdE(i)=10.0/6.0*Ex_nn*dEx_nn + (Ex_nn*dEx_p1+Ex_p1*dEx_nn+ ...
    %                                Ex_nn*dEx_m1+Ex_m1*dEx_nn)/12.0;
    %EdE(i)=Ex_nn*dEx_nn + dEx_nn*(5/6*Ex_nn + 1/12*Ex_m1 +
    %1/12*Ex_p1);
    EdE(i)=(5/6*Ex_nn + 1/12*Ex_m1 + 1/12*Ex_p1)*dEx_nn;
    val = val + EdE(i);
end
val*ccelec/(8*pi)
EdE=EdE*ccelec/(8*pi);
figure(2);
plot(EdE);
return;

phiLES=zeros(20,1);
phiLES(1) = -interp3(m,m,m,ExV,19,11,11) - Exdip;
for i=2:nx
    phiLES(i) = phiLES(i-1) - interp3(m,m,m,ExV,i-1,11,11) - Exdip;
end
figure(2);
plot((0:19)+0.5,phiLES,'.-b');
hold on;
plot(0:19,interp3(m,m,m,ExV,0:19,11,11),'-.b');
plot(0:19,interp3(m,m,m,phiV,0:19,11,11),'.-r');
plot((0:19)+0.5,interp3(m,m,m,ExV2,0:19,11,11),'-.r');
hold off;
xlabel('x');
ylabel('phi');

phiLES=zeros(20,1);
phiLES(1) = -interp3(m,m,m,EyV,10,19,11);
for i=2:ny
    phiLES(i) = phiLES(i-1) - interp3(m,m,m,EyV,10,i-1,11);
end
figure(3);
plot((0:19)+0.5,phiLES,'.-b');
hold on;
plot(0:19,interp3(m,m,m,EyV,10,0:19,11),'-.b');
plot(0:19,interp3(m,m,m,phiV,10,0:19,11),'.-r');
plot((0:19)+0.5,interp3(m,m,m,EyV2,10,0:19,11),'-.r');
hold off;
xlabel('y');
ylabel('phi');

phiLES=zeros(nx,ny,nz);
phivec=zeros(nx*ny*nz,1);
i=1;
for iz=1:nz
    for iy=1:ny
        phiLES(1,iy,iz) = -interp3(m,m,m,ExV,nx-1,iy,iz) - Exdip;
        phivec(i) = phiLES(1,iy,iz);
        i = i + 1;
        for ix=2:nx
            phiLES(ix,iy,iz) = phiLES(ix-1,iy,iz) - interp3(m,m,m,ExV,ix-1,iy,iz) ...
                - Exdip;
            phivec(i) = phiLES(ix,iy,iz);
            i = i + 1;
        end
    end
end
