A = load('L.txt');
b = load('rho.txt');

D = diag(-6*ones(1,length(b)));
R = A - D;

%x = gauss_seidel(A, b, b*0, 100);

x_gs = gs3d(b, b*0, 100, 4, 4, 4);
printf('gs3d error %e\n',max(abs(A*x_gs-b)));

x_rbgs = rbgs3d(b, b*0, 100, 4, 4, 4);
printf('rbgs3d error %e\n',max(abs(A*x_rbgs-b)));

return;

x = b*0;
for i=1:100000
    x=inv(D)*(b - R*x);
end