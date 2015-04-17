function [wx, wy, wz, jlinkx, jlinky, jlinkz] = calcWeights(nx, ny, nz, dx, dy, dz, sigma, h, s);

    pref = (2*pi*sigma^2)^(-1/2);
    wx = zeros(2*nx+1, 1);
    dwx = zeros(2*nx+1, 1);
    wy = zeros(2*ny+1, 1);
    dwy = zeros(2*ny+1, 1);
    wy = zeros(2*nz+1, 1);
    dwy = zeros(2*nz+1, 1);
    for j=-nx:nx
        wx(j+nx+1) = pref*exp(-(dx+s+j)^2*h^2/(2*sigma^2));
        dwx(j+nx+1) = ((-dx-j)*h^2/(sigma^2))*pref*exp(-(dx+j)^2*h^2/(2*sigma^2));
    end
    for j=-ny:ny
        wy(j+ny+1) = pref*exp(-(dy+s+j)^2*h^2/(2*sigma^2));
        dwy(j+ny+1) = ((-dy-j)*h^2/(sigma^2))*pref*exp(-(dy+j)^2*h^2/(2*sigma^2));
    end
    for j=-nz:nz
        wz(j+nz+1) = pref*exp(-(dz+s+j)^2*h^2/(2*sigma^2));
        dwz(j+nz+1) = ((-dz-j)*h^2/(sigma^2))*pref*exp(-(dz+j)^2*h^2/(2*sigma^2));
    end
    
    jlinkx = zeros(2*nx+1, 1);
    jlinky = zeros(2*ny+1, 1);
    jlinkz = zeros(2*nz+1, 1);
    jlinkx(1) = dwx(1);
    for j=-nx+1:nx
        jlinkx(j+nx+1) = jlinkx(j-1+nx+1) + dwx(j+nx+1);
    end
    jlinky(1) = dwy(1);
    for j=-ny+1:ny
        jlinky(j+ny+1) = jlinky(j-1+ny+1) + dwy(j+ny+1);
    end
    jlinkz(1) = dwz(1);
    for j=-nz+1:nz
        jlinkz(j+nz+1) = jlinkz(j-1+nz+1) + dwz(j+nz+1);
    end

end
