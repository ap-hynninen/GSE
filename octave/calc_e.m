L=20;
nx=20;
ny=20;
nz=20;
ccelec=332.0716;

ExLES=load('ExLES.txt');
EyLES=load('EyLES.txt');
EzLES=load('EzLES.txt');

ExLESd=load('ExLESd.txt');
EyLESd=load('EyLESd.txt');
EzLESd=load('EzLESd.txt');

ExLES = reshape(ExLES,20,20,20);
EyLES = reshape(EyLES,20,20,20);
EzLES = reshape(EzLES,20,20,20);

ExLESd = reshape(ExLESd,20,20,20);
EyLESd = reshape(EyLESd,20,20,20);
EzLESd = reshape(EzLESd,20,20,20);

U = calcu(nx,ny,nz,ExLES,EyLES,EzLES);
Ud = calcu(nx,ny,nz,ExLESd,EyLESd,EzLESd);

U = U*ccelec/(8*pi);
U

Ud = Ud*ccelec/(8*pi);
Ud

E=load("E_LES.txt");
ix=zeros(nx*ny*nz,1);
iy=zeros(nx*ny*nz,1);
iz=zeros(nx*ny*nz,1);
i=1;
for z=1:nz
    for y=1:ny
        for x=1:nx
            ix(i) = x-1;
            iy(i) = y-1;
            iz(i) = z-1;
            i = i + 1;
        end
    end
end

ix=E(1:368,1);
iy=E(1:368,2);
iz=E(1:368,3);

dU=calcdu(nx,ny,nz,ExLES,EyLES,EzLES,ix,iy,iz);
dU = dU*ccelec/(8*pi);
dU

(5/6*ExLES(10,12,12) + 1/12*ExLES(9,12,12) + 1/12*ExLES(11,12,12))*ccelec