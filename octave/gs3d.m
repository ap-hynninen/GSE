function [xx] = gs3d(b, x0, iters, nx, ny, nz)
    nxy = nx*ny;
    xx = x0;
    for k = 1:iters
        for z = 0:nz-1
            zp = mod(z+1,nz);
            zm = mod(z-1+nz,nz);
            for y = 0:ny-1
                yp = mod(y+1,ny);
                ym = mod(y-1+ny,ny);
                for x = 0:nx-1
                    xp = mod(x+1,nx);
                    xm = mod(x-1+nx,nx);
                    xx(x+y*nx+z*nxy+1) = (xx(xm+y*nx+z*nxy+1) + ...
                                          xx(xp+y*nx+z*nxy+1) + ...
                                          xx(x+ym*nx+z*nxy+1) + ...
                                          xx(x+yp*nx+z*nxy+1) + ...
                                          xx(x+y*nx+zm*nxy+1) + ...
                                          xx(x+y*nx+zp*nxy+1) - ...
                                          b(x+y*nx+z*nxy+1))/6.0;
                end
            end
        end
    end
end
