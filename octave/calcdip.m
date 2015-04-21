clear;
r=5.0;
a=r/(2.0*sqrt(3.0));
L=20;
M=20;
h=1;
rcut = 6.364;
rcutsq = rcut*rcut;
sigma = 2.12132034355964/sqrt(2.0);
x=zeros(2,1);
y=zeros(2,1);
z=zeros(2,1);
q=zeros(2,1);

x(1) = -a + 0.5*L;
y(1) = -a + 0.5*L;
z(1) = -a + 0.5*L;
x(2) = a + 0.5*L;
y(2) = a + 0.5*L;
z(2) = a + 0.5*L;
q(1) = -1.0;
q(2) = 1.0;

dipx=0;
dipy=0;
dipz=0;
for k=0:M-1
    dz = z-k*h;
    wz = exp(-dz.^2/(2*sigma^2))*(2*pi*sigma^2)^(-3/2);
    for j=0:M-1
        dy = y-j*h;
        wy = exp(-dy.^2/(2*sigma^2));
        for i=0:M-1
            dx = x-i*h;
            wx = exp(-dx.^2/(2*sigma^2));
            r = sqrt(dx.^2 + dy.^2 + dz.^2);
            if (r(1) < rcutsq && r(2) < rcutsq)
                
            end
            %rho = sum()
            %dipx = dipx + sum(q.*wx.*wy.*wz*i*h);
        end
    end
end