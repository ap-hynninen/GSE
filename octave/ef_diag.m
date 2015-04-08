clear;
L=20;
ccelec=332.0716;
r=1.0;
sigma = 1.5/sqrt(2);
rcut = 4.5;
a = r/(2.0*sqrt(3.0));
xyz=[-a + L/2, -a + L/2, -a + L/2;
     a + L/2,  a + L/2, a + L/2];

ExLES=load('../ExLES.txt');
EyLES=load('../EyLES.txt');
EzLES=load('../EzLES.txt');

N=round(length(ExLES)^(1/3));
h = L/N;

ExLES = reshape(ExLES,N,N,N);
EyLES = reshape(EyLES,N,N,N);
EzLES = reshape(EzLES,N,N,N);

rcut = 6.0;
[Ex Ey Ez] = calcElectricField(xyz, sigma, rcut, h, ExLES, EyLES, ...
                               EzLES);
Ex=Ex*ccelec;
Ey=Ey*ccelec;
Ez=Ez*ccelec;
