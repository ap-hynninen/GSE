function U=calcusmp(nx,ny,nz,Ex,Ey,Ez)

    U = 0;
    for z=0:nz-1
        for y=0:ny-1
            for x=0:nx-1
                Ex_nn = Ex(x+1,y+1,z+1);
                Ux = Ex_nn;
                
                %Ey_nn = Ey(x+1,y+1,z+1);
                %Uy = Ey_nn^2;
                
                %Ez_nn = Ez(x+1,y+1,z+1);
                %Uz = Ez_nn^2;
                
                U = U + Ux;% + Uy + Uz;
            end
        end
    end

end
