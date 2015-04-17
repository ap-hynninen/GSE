clear;
ccelec = 332.0716;

ExLES=load('../ExLES.txt');
EyLES=load('../EyLES.txt');
EzLES=load('../EzLES.txt');

sigma = 2.12132034355964/sqrt(2.0);
L = 20;
M = 20;
h = 1;

Ex=reshape(ExLES,M,M,M);
Ey=reshape(EyLES,M,M,M);
Ez=reshape(EzLES,M,M,M);

ix=7;
iy=10;
iz=10;

dx=0;
dy=0;
dz=0;

nx=7;
ny=7;
nz=7;

Jx = zeros(2*nx+1, 2*ny+1, 2*nz+1);
Jy = zeros(2*nx+1, 2*ny+1, 2*nz+1);
Jz = zeros(2*nx+1, 2*ny+1, 2*nz+1);

[wx, wy, wz, jlinkx, jlinky, jlinkz] = calcWeights(nx, ny, nz, dx, dy, dz, sigma, h, 0);
for tz=-nz:nz
    for ty=-ny:ny
        for tx=-nx:nx
            Jx(tx+nx+1,ty+ny+1,tz+nz+1) = jlinkx(tx+nx+1)*wy(ty+ny+1)*wz(tz+nz+1);
            Jy(tx+nx+1,ty+ny+1,tz+nz+1) = jlinky(ty+ny+1)*wx(tx+nx+1)*wz(tz+nz+1);
            Jz(tx+nx+1,ty+ny+1,tz+nz+1) = jlinkz(tz+nz+1)*wy(ty+ny+1)*wx(tx+nx+1);
        end
    end
end

%ExDip = 5.214182e-01/ccelec;
ExDip = 3.128039/ccelec;

AK = (97.0/120.0);
BK = (12.0/120.0);
CK = (-1.0/240.0);
fx = 0;
fy = 0;
fz = 0;
for tz=-nz:nz
    oz = mod(iz+tz+M, M);
    for ty=-ny:ny
        oy = mod(iy+ty+M, M);
        for tx=-nx:nx
            ox = mod(ix+tx+M, M);

            Exm2 = Ex(mod(ox-2+M,M)+1,oy+1,oz+1);
            Exm  = Ex(mod(ox-1+M,M)+1,oy+1,oz+1);
            Exn  = Ex(ox+1,oy+1,oz+1);
            Exp  = Ex(mod(ox+1,M)+1,oy+1,oz+1);
            Exp2 = Ex(mod(ox+2,M)+1,oy+1,oz+1);

            Exm2 = Exm2 + ExDip;
            Exm  = Exm + ExDip;
            Exn  = Exn + ExDip;
            Exp  = Exp + ExDip;
            Exp2 = Exp2 + ExDip;
            
            Eym2 = Ey(ox+1,mod(oy-2+M,M)+1,oz+1);
            Eym  = Ey(ox+1,mod(oy-1+M,M)+1,oz+1);
            Eyn  = Ey(ox+1,oy+1,oz+1);
            Eyp  = Ey(ox+1,mod(oy+1,M)+1,oz+1);
            Eyp2 = Ey(ox+1,mod(oy+2,M)+1,oz+1);

            Ezm2 = Ez(ox+1,oy+1,mod(oz-2+M,M)+1);
            Ezm  = Ez(ox+1,oy+1,mod(oz-1+M,M)+1);
            Ezn  = Ez(ox+1,oy+1,oz+1);
            Ezp  = Ez(ox+1,oy+1,mod(oz+1,M)+1);
            Ezp2 = Ez(ox+1,oy+1,mod(oz+2,M)+1);

            Jxn = Jx(tx+nx+1,ty+ny+1,tz+nz+1);
            Jyn = Jy(tx+nx+1,ty+ny+1,tz+nz+1);
            Jzn = Jz(tx+nx+1,ty+ny+1,tz+nz+1);
            
            Jxm = 0;
            if (tx-1 >= -nx)
                Jxm = Jx(tx+nx-1+1,ty+ny+1,tz+nz+1);
            end
            Jym = 0;
            if (ty-1 >= -ny)
                Jym = Jy(tx+nx+1,ty+ny-1+1,tz+nz+1);
            end
            Jzm = 0;
            if (tz-1 >= -nz)
                Jzm = Jz(tx+nx+1,ty+ny+1,tz+nz-1+1);
            end

            Jxp = 0;
            if (tx+1 <= nx)
                Jxp = Jx(tx+nx+1+1,ty+ny+1,tz+nz+1);
            end
            Jyp = 0;
            if (ty+1 <= ny)
                Jyp = Jy(tx+nx+1,ty+ny+1+1,tz+nz+1);
            end
            Jzp = 0;
            if (tz+1 <= nz)
                Jzp = Jz(tx+nx+1,ty+ny+1,tz+nz+1+1);
            end

            Jxm2 = 0;
            if (tx-2 >= -nx)
                Jxm2 = Jx(tx+nx-2+1,ty+ny+1,tz+nz+1);
            end
            Jym2 = 0;
            if (ty-2 >= -ny)
                Jym2 = Jy(tx+nx+1,ty+ny-2+1,tz+nz+1);
            end
            Jzm2 = 0;
            if (tz-2 >= -nz)
                Jzm2 = Jz(tx+nx+1,ty+ny+1,tz+nz-2+1);
            end

            Jxp2 = 0;
            if (tx+2 <= nx)
                Jxp2 = Jx(tx+nx+2+1,ty+ny+1,tz+nz+1);
            end
            Jyp2 = 0;
            if (ty+2 <= ny)
                Jyp2 = Jy(tx+nx+1,ty+ny+2+1,tz+nz+1);
            end
            Jzp2 = 0;
            if (tz+2 <= nz)
                Jzp2 = Jz(tx+nx+1,ty+ny+1,tz+nz+2+1);
            end

            fx = fx + Jxn*(AK*Exn + 0.5*BK*(Exm + Exp) + 0.5*CK*(Exm2 + Exp2)) + Exn*(0.5*BK*(Jxm + Jxp) + 0.5*CK*(Jxm2 + Jxp2));
            fy = fy + Jyn*(AK*Eyn + 0.5*BK*(Eym + Eyp) + 0.5*CK*(Eym2 + Eyp2)) + Eyn*(0.5*BK*(Jym + Jyp) + 0.5*CK*(Jym2 + Jyp2));
            fz = fz + Jzn*(AK*Ezn + 0.5*BK*(Ezm + Ezp) + 0.5*CK*(Ezm2 + Ezp2)) + Ezn*(0.5*BK*(Jzm + Jzp) + 0.5*CK*(Jzm2 + Jzp2));
        end
    end
end


fx = -fx*ccelec;
fy = -fy*ccelec;
fz = -fz*ccelec;

[fx fy fz]