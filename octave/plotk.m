L=10;
k=2*pi*(1:4)/L;
sig=2.12;
sig1=1.34081;
b=2.0;
N=4;

f=exp(-0.5*k.^2*(sig^2 - 2*sig1^2))*4*pi./k;

g=0;
for j=0:N-1
    g=g+exp(-0.5*k.^2*(sig^2*b^(2*j) - 2*sig1^2));
end
g=g*2*log(b);

plot(k,f);
hold on;
plot(k,g,'r');
hold off;