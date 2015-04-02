function [xx] = gs2d(b, x0, iters, nx, ny)
    xx = x0;
    for k = 1:iters
        for y = 0:ny-1
            yp = mod(y+1,ny);
            ym = mod(y-1+ny,ny);
            for x = 0:nx-1
                xp = mod(x+1,nx);
                xm = mod(x-1+nx,nx);
                xx(x+y*nx+1) = (xx(xm+y*nx+1) + ...
                                xx(xp+y*nx+1) + ...
                                xx(x+ym*nx+1) + ...
                                xx(x+yp*nx+1) - ...
                                b(x+y*nx+1))/4.0;
            end
        end
    end
end
