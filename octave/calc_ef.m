clear;
L=20;
ccelec=332.0716;
r=4.0;
sigma = 1.5/sqrt(2);
rcut = 4.5;
xyz=[-r/2 + L/2, L/2, L/2;
     r/2 + L/2,  L/2, L/2];

ExLES0=load('../r4/ExLES0.txt');
EyLES0=load('../r4/EyLES0.txt');
EzLES0=load('../r4/EzLES0.txt');

N=round(length(ExLES0)^(1/3));
h = L/N;

ExLES0 = reshape(ExLES0,N,N,N);
EyLES0 = reshape(EyLES0,N,N,N);
EzLES0 = reshape(EzLES0,N,N,N);

ExLESdx=load('../r4/ExLESdx.txt');
EyLESdx=load('../r4/EyLESdx.txt');
EzLESdx=load('../r4/EzLESdx.txt');
ExLESdx = reshape(ExLESdx,N,N,N);
EyLESdx = reshape(EyLESdx,N,N,N);
EzLESdx = reshape(EzLESdx,N,N,N);

ExLESdy=load('../r4/ExLESdy.txt');
EyLESdy=load('../r4/EyLESdy.txt');
EzLESdy=load('../r4/EzLESdy.txt');
ExLESdy = reshape(ExLESdy,N,N,N);
EyLESdy = reshape(EyLESdy,N,N,N);
EzLESdy = reshape(EzLESdy,N,N,N);

ExLESdz=load('../r4/ExLESdz.txt');
EyLESdz=load('../r4/EyLESdz.txt');
EzLESdz=load('../r4/EzLESdz.txt');
ExLESdz = reshape(ExLESdz,N,N,N);
EyLESdz = reshape(EyLESdz,N,N,N);
EzLESdz = reshape(EzLESdz,N,N,N);

S1 = markSupport(xyz(1,:), sigma, rcut, h, N);
S2 = markSupport(xyz(2,:), sigma, rcut, h, N);

xm = zeros(N*N,1);
ym = zeros(N*N,1);
i=1;
for y=0:N-1
    for x=0:N-1
        xm(i) = (x+0.5)*h;
        ym(i) = (y+0.5)*h;
        i = i + 1;
    end
end
plot(xyz(1,1),xyz(1,2),'*k');
hold on;
plot(xyz(2,1),xyz(2,2),'*k');
f=5;
quiver(xm,ym,vec(ExLES0(:,:,N/2+1))*f,vec(EyLES0(:,:,N/2+1))*f,0,'b');
hold off;

xlabel('x');
ylabel('y');
axis([5 15 5 15]);

% $$$ [FnX FnY FnZ] = calcNumForce(0.001, ExLES0, EyLES0, EzLES0, ExLESdx, ...
% $$$                              ExLESdy, ExLESdz, EyLESdx, EyLESdy, ...
% $$$                              EyLESdz, EzLESdx, EzLESdy, EzLESdz);
% $$$ FnX=FnX*ccelec;
% $$$ FnY=FnY*ccelec;
% $$$ FnZ=Fnz*ccelec;

[Ex Ey Ez] = calcElectricField(xyz, sigma, rcut, h, ExLES0, EyLES0, ...
                               EzLES0);
Ex=Ex*ccelec;
Ey=Ey*ccelec;
Ez=Ez*ccelec;

% $$$ [Fx Fy Fz] = calcForce(xyz, sigma, rcut, h, ExLES, EyLES, ...
% $$$                        EzLES);
% $$$ Fx=Fx*ccelec;
% $$$ Fy=Fy*ccelec;
% $$$ Fz=Fz*ccelec;

return;

