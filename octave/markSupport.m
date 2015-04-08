function [S] = markSupport(xyz, sigma, rcut, h, N)
    rcut2 = rcut^2;
    nk = ceil(rcut/h);
    S = zeros(N,N,N);
    for i=1:size(xyz,1)
        x = xyz(1,1);
        y = xyz(1,2);
        z = xyz(1,3);
        ix = round(x/h-0.5);
        iy = round(y/h-0.5);
        iz = round(z/h-0.5);
        for lz=iz-nk:iz+nk
            for ly=iy-nk:iy+nk
                for lx=ix-nk:ix+nk
                    r2 = ((lx+0.5)*h-x)^2 + ((ly+0.5)*h-y)^2 + ((lz+0.5)*h-z)^2;
                    if (r2 < rcut2)
                        fac = exp(-r2/(2*sigma*sigma));
                        nx = mod(lx+N,N);
                        ny = mod(ly+N,N);
                        nz = mod(lz+N,N);
                        S(nx+1,ny+1,nz+1) = fac;
                    end
                end
            end
        end
    end
end
