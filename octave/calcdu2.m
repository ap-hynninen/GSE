function U=calcdu2(nx,ny,nz,Ex,Ey,Ez,Exd,Eyd,Ezd,ix,iy,iz)
    U = 0;
    for i=1:length(ix)
        x = ix(i);
        y = iy(i);
        z = iz(i);
        
        Ex_m2 = Ex(mod(x-2+nx,nx)+1,y+1,z+1);
        Ex_m1 = Ex(mod(x-1+nx,nx)+1,y+1,z+1);
        Ex_nn = Ex(x+1,y+1,z+1);
        Ex_p1 = Ex(mod(x+1,nx)+1,y+1,z+1);
        Ex_p2 = Ex(mod(x+2,nx)+1,y+1,z+1);
        Ex_p3 = Ex(mod(x+3,nx)+1,y+1,z+1);
        
        dEx_m1 = Ex_nn-Ex_m1;
        dEx_nn = Ex_p1-Ex_nn;
        dEx_p1 = Ex_p2-Ex_p1;

        Ux = 10/6*Ex_nn*dEx_nn + 1/12*(Ex_nn*dEx_m1 + Ex_m1*dEx_nn ...
                                       + Ex_nn*dEx_p1 + Ex_p1*dEx_nn);
                
        Ey_m2 = Ey(mod(x-2+nx,nx)+1,y+1,z+1);
        Ey_m1 = Ey(mod(x-1+nx,nx)+1,y+1,z+1);
        Ey_nn = Ey(x+1,y+1,z+1);
        Ey_p1 = Ey(mod(x+1,nx)+1,y+1,z+1);
        Ey_p2 = Ey(mod(x+2,nx)+1,y+1,z+1);
        Ey_p3 = Ey(mod(x+3,nx)+1,y+1,z+1);

        dEy_m1 = Ey_nn-Ey_m1;
        dEy_nn = Ey_p1-Ey_nn;
        dEy_p1 = Ey_p2-Ey_p1;

        Uy = 10/6*Ey_nn*dEy_nn + 1/12*(Ey_nn*dEy_m1 + Ey_m1*dEy_nn ...
                                       + Ey_nn*dEy_p1 + Ey_p1*dEy_nn);
        
        Ez_m2 = Ez(mod(x-2+nx,nx)+1,y+1,z+1);
        Ez_m1 = Ez(mod(x-1+nx,nx)+1,y+1,z+1);
        Ez_nn = Ez(x+1,y+1,z+1);
        Ez_p1 = Ez(mod(x+1,nx)+1,y+1,z+1);
        Ez_p2 = Ez(mod(x+2,nx)+1,y+1,z+1);
        Ez_p3 = Ez(mod(x+3,nx)+1,y+1,z+1);

        dEz_m1 = Ez_nn-Ez_m1;
        dEz_nn = Ez_p1-Ez_nn;
        dEz_p1 = Ez_p2-Ez_p1;

        Uz = 10/6*Ez_nn*dEz_nn + 1/12*(Ez_nn*dEz_m1 + Ez_m1*dEz_nn ...
                                       + Ez_nn*dEz_p1 + Ez_p1*dEz_nn);


        U = U + Ux + Uy + Uz;
    end

end
