rho=load('rho.txt');
t=load('rhoK1.txt');
rhoK1ref = complex(t(1:2:length(t)), t(2:2:length(t)));

rho_3d=reshape(rho, 64, 64, 64);
rhoK1ref_3d=reshape(rhoK1ref, 33, 64, 64);

rhoK1_3d=fftn(rho_3d);

max(max(max(abs(rhoK1ref_3d-rhoK1_3d(1:33,:,:)))))

