function U=calcu(nx,ny,nz,Ex,Ey,Ez)

    U = 0;
    for z=0:nz-1
        for y=0:ny-1
            for x=0:nx-1
                Ex_m1 = Ex(mod(x-1+nx,nx)+1,y+1,z+1);
                Ex_nn = Ex(x+1,y+1,z+1);
                Ex_p1 = Ex(mod(x+1,nx)+1,y+1,z+1);
                Ux = 5/6*Ex_nn^2 + 1/12*Ex_nn*(Ex_m1 + Ex_p1);
                
                Ey_m1 = Ey(x+1,mod(y-1+ny,ny)+1,z+1);
                Ey_nn = Ey(x+1,y+1,z+1);
                Ey_p1 = Ey(x+1,mod(y+1,ny)+1,z+1);
                Uy = 5/6*Ey_nn^2 + 1/12*Ey_nn*(Ey_m1 + Ey_p1);
                
                Ez_m1 = Ez(x+1,y+1,mod(z-1+nz,nz)+1);
                Ez_nn = Ez(x+1,y+1,z+1);
                Ez_p1 = Ez(x+1,y+1,mod(z+1,nz)+1);
                Uz = 5/6*Ez_nn^2 + 1/12*Ez_nn*(Ez_m1 + Ez_p1);
                
                U = U + Ux + Uy + Uz;
            end
        end
    end

end
