clear;

h=1;
for i=1:4
    nx=7/h;
    ny=7/h;
    nz=7/h;
    dx=0;
    dy=0;
    dz=0;
    sigma = 2.12132034355964/sqrt(2.0);
    
    [wx, wy, wz, jlinkx, jlinky, jlinkz] = calcWeights(nx, ny, nz, dx, dy, dz, sigma, h, 0);
    
    if (h==1)
        c = '.-b';
    elseif (h==0.5)
        c = '.-r';
    elseif (h==0.25)
        c = '.-g';
    else
        c = '.-k';
    end
    %figure;
    %plot((-nx:nx)*h,wx,c);
    plot((-nx:nx)*h,jlinkx,c);
    hold on;
    h = h/2;
end
hold off;

