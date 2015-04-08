function [Ex Ey Ez] = calcElectricField(xyz, sigma, rcut, h, ExM, EyM, EzM)
    N = size(ExM,1);
    rcut2 = rcut^2;
    nk = ceil(rcut/h);
    Ex = zeros(size(xyz,1),1);
    Ey = zeros(size(xyz,1),1);
    Ez = zeros(size(xyz,1),1);
    for i=1:size(xyz,1)
        x = xyz(i,1);
        y = xyz(i,2);
        z = xyz(i,3);
        ix = round(x/h-0.5);
        iy = round(y/h-0.5);
        iz = round(z/h-0.5);
        printf('x,y,z = %f %f %f ix,iy,iz = %d %d %d\n',x,y,z,ix,iy,iz);
        for lz=iz-nk:iz+nk
            for ly=iy-nk:iy+nk
                for lx=ix-nk:ix+nk
                    r2 = ((lx+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    if (r2 < rcut2)
                        fac = exp(-r2/(2*sigma*sigma));
                        nx = mod(lx+N,N);
                        ny = mod(ly+N,N);
                        nz = mod(lz+N,N);
                        Ex(i) = Ex(i) + ExM(nx+1,ny+1,nz+1)*fac;
                        Ey(i) = Ey(i) + EyM(nx+1,ny+1,nz+1)*fac;
                        Ez(i) = Ez(i) + EzM(nx+1,ny+1,nz+1)*fac;
                    end
                end
            end
        end
        Ex(i) = Ex(i)*(2*pi*sigma*sigma)^(-3/2)*h*h*h;
        Ey(i) = Ey(i)*(2*pi*sigma*sigma)^(-3/2)*h*h*h;
        Ez(i) = Ez(i)*(2*pi*sigma*sigma)^(-3/2)*h*h*h;
    end
end
