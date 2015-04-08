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
        ix = round(x/h);
        iy = round(y/h);
        iz = round(z/h);
        printf('x,y,z = %f %f %f ix,iy,iz = %d %d %d\n',x,y,z,ix,iy,iz);
        for lz=iz-nk:iz+nk
            for ly=iy-nk:iy+nk
                for lx=ix-nk:ix+nk
                    r2_mx = ((lx-1+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    r2_my = ((lx+0.5)*h-x)^2 + ((ly-1+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    r2_mz = ((lx+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz-1+0.5)*h-z)^2;
                    r2_nn = ((lx+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    r2_px = ((lx+1+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    r2_py = ((lx+0.5)*h-x)^2 + ((ly+1+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    r2_pz = ((lx+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz+1+0.5)*h-z)^2;
                    
                    w_mx = exp(-r2_mx/(2*sigma*sigma));
                    w_my = exp(-r2_my/(2*sigma*sigma));
                    w_mz = exp(-r2_mz/(2*sigma*sigma));
                    w_nn = exp(-r2_nn/(2*sigma*sigma));
                    w_px = exp(-r2_px/(2*sigma*sigma));
                    w_py = exp(-r2_py/(2*sigma*sigma));
                    w_pz = exp(-r2_pz/(2*sigma*sigma));
                    
                    mx = mod(lx-1+N,N);
                    my = mod(ly-1+N,N);
                    mz = mod(lz-1+N,N);
                    nx = mod(lx+N,N);
                    ny = mod(ly+N,N);
                    nz = mod(lz+N,N);
                    px = mod(lx+1+N,N);
                    py = mod(ly+1+N,N);
                    pz = mod(lz+1+N,N);
                    
                    Ex_mx = ExM(mx+1,ny+1,nz+1);
                    Ey_my = EyM(nx+1,my+1,nz+1);
                    Ez_mz = EzM(nx+1,ny+1,mz+1);
                    Ex_nn = ExM(nx+1,ny+1,nz+1);
                    Ey_nn = EyM(nx+1,ny+1,nz+1);
                    Ez_nn = EzM(nx+1,ny+1,nz+1);
                    Ex_px = ExM(px+1,ny+1,nz+1);
                    Ey_py = EyM(nx+1,py+1,nz+1);
                    Ez_pz = EzM(nx+1,ny+1,pz+1);
                    
                    Ex(i) = Ex(i) + w_nn*(5/6*Ex_nn + 1/24*Ex_mx + 1/24*Ex_px) + 1/24*(w_mx+w_px)*Ex_nn;
                    Ey(i) = Ey(i) + w_nn*(5/6*Ey_nn + 1/24*Ey_my + 1/24*Ey_py) + 1/24*(w_my+w_py)*Ey_nn;
                    Ez(i) = Ez(i) + w_nn*(5/6*Ez_nn + 1/24*Ez_mz + 1/24*Ez_pz) + 1/24*(w_mz+w_pz)*Ez_nn;
                end
            end
        end
        Ex(i) = Ex(i)*(2*pi*sigma*sigma)^(-3/2)*h*h*h;
        Ey(i) = Ey(i)*(2*pi*sigma*sigma)^(-3/2)*h*h*h;
        Ez(i) = Ez(i)*(2*pi*sigma*sigma)^(-3/2)*h*h*h;
    end
end
