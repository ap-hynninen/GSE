function [] = calcNumForce(d, Ex, Ey, Ez, Exdx, Exdy, Exdz, Eydx, Eydy, Eydz, Ezdx, Ezdy, Ezdz);

    dEx_dxi=(Exdx-Ex)/d;
    dEx_dyi=(Exdy-Ex)/d;
    dEx_dzi=(Exdz-Ex)/d;

    dEy_dxi=(Eydx-Ey)/d;
    dEy_dyi=(Eydy-Ey)/d;
    dEy_dzi=(Eydz-Ey)/d;

    dEz_dxi=(Ezdx-Ez)/d;
    dEz_dyi=(Ezdy-Ez)/d;
    dEz_dzi=(Ezdz-Ez)/d;

    F=[0 0 0];
    for z=0:N-1
        for y=0:N-1
            for x=0:N-1
                mx = mod(x-1+N,N);
                my = mod(y-1+N,N);
                mz = mod(z-1+N,N);
                px = mod(x+1,N);
                py = mod(y+1,N);
                pz = mod(z+1,N);

                Ex_m = Ex(mx+1,y+1,z+1);
                Ex_n = Ex(x+1,y+1,z+1);
                Ex_p = Ex(px+1,y+1,z+1);

                Ey_m = Ey(x+1,my+1,z+1);
                Ey_n = Ey(x+1,y+1,z+1);
                Ey_p = Ey(x+1,py+1,z+1);

                Ez_m = Ez(x+1,y+1,mz+1);
                Ez_n = Ez(x+1,y+1,z+1);
                Ez_p = Ez(x+1,y+1,pz+1);

                DEx_m=[dEx_dxi(mx+1,y+1,z+1) dEx_dyi(mx+1,y+1,z+1) dEx_dzi(mx+1,y+1,z+1)];
                DEx_n=[dEx_dxi(x+1,y+1,z+1)  dEx_dyi(x+1,y+1,z+1)  dEx_dzi(x+1,y+1,z+1)];
                DEx_p=[dEx_dxi(px+1,y+1,z+1) dEx_dyi(px+1,y+1,z+1) dEx_dzi(px+1,y+1,z+1)];
                
                DEy_m=[dEy_dxi(x+1,my+1,z+1) dEy_dyi(x+1,my+1,z+1) dEy_dzi(x+1,my+1,z+1)];
                DEy_n=[dEy_dxi(x+1,y+1,z+1)  dEy_dyi(x+1,y+1,z+1)  dEy_dzi(x+1,y+1,z+1)];
                DEy_p=[dEy_dxi(x+1,py+1,z+1) dEy_dyi(x+1,py+1,z+1) dEy_dzi(x+1,py+1,z+1)];

                DEz_m=[dEz_dxi(x+1,y+1,mz+1) dEz_dyi(x+1,y+1,mz+1) dEz_dzi(x+1,y+1,mz+1)];
                DEz_n=[dEz_dxi(x+1,y+1,z+1)  dEz_dyi(x+1,y+1,z+1)  dEz_dzi(x+1,y+1,z+1)];
                DEz_p=[dEz_dxi(x+1,y+1,pz+1) dEz_dyi(x+1,y+1,pz+1) dEz_dzi(x+1,y+1,pz+1)];

                t = ...
                    DEx_n*(5/6*Ex_n + 1/12*Ex_m + 1/12*Ex_p) + ...
                    Ex_n*(5/6*DEx_n + 1/12*DEx_m + 1/12*DEx_p) + ...
                    DEy_n*(5/6*Ey_n + 1/12*Ey_m + 1/12*Ey_p) + ...
                    Ey_n*(5/6*DEy_n + 1/12*DEy_m + 1/12*DEy_p) + ...
                    DEz_n*(5/6*Ez_n + 1/12*Ez_m + 1/12*Ez_p) + ...
                    Ez_n*(5/6*DEz_n + 1/12*DEz_m + 1/12*DEz_p);

                F = F + t;
            end
        end
    end
    F = -F/(8*pi);

end
