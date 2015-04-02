function U=calcdusmp(nx,ny,nz,Ex,Ey,Ez,ix,iy,iz)
    U = 0;
    for i=1:length(ix)
        x = ix(i);
        y = iy(i);
        z = iz(i);
        
        Ex_m1 = Ex(mod(x-1+nx,nx)+1,y+1,z+1);
        Ex_nn = Ex(x+1,y+1,z+1);
        Ex_p1 = Ex(mod(x+1,nx)+1,y+1,z+1);
        Ex_p2 = Ex(mod(x+2,nx)+1,y+1,z+1);
        %dEx_nn = Ex_p1-Ex_nn;
        dEx_nn = 5/6*(Ex_p1-Ex_nn) + 1/12*(Ex_p2-Ex_p1) + 1/12*(Ex_nn-Ex_m1);
        %Ux = 2*Ex_nn*dEx_nn;
        Ux = dEx_nn;
        
        Ey_m1 = Ey(mod(x-1+nx,nx)+1,y+1,z+1);
        Ey_nn = Ey(x+1,y+1,z+1);
        Ey_p1 = Ey(mod(x+1,nx)+1,y+1,z+1);
        Ey_p2 = Ey(mod(x+2,nx)+1,y+1,z+1);
        %dEy_nn = Ey_p1-Ey_nn;
        dEy_nn = 5/6*(Ey_p1-Ey_nn) + 1/12*(Ey_p2-Ey_p1) + 1/12*(Ey_nn-Ey_m1);
        Uy = 2*Ey_nn*dEy_nn;
        
        Ez_m1 = Ez(mod(x-1+nx,nx)+1,y+1,z+1);
        Ez_nn = Ez(x+1,y+1,z+1);
        Ez_p1 = Ez(mod(x+1,nx)+1,y+1,z+1);
        Ez_p2 = Ez(mod(x+2,nx)+1,y+1,z+1);
        %dEz_nn = Ez_p1-Ez_nn;
        dEz_nn = 5/6*(Ez_p1-Ez_nn) + 1/12*(Ez_p2-Ez_p1) + 1/12*(Ez_nn-Ez_m1);
        Uz = 2*Ez_nn*dEz_nn;

        U = U + Ux;% + Uy + Uz;
    end

end
