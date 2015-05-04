function [y] = dispgrid(d, x);
    n = round(length(x)^(1/3));
    t = reshape(x, n, n, n);
    y = zeros(n, n, n);
    for zs=0:n-1
        zd = mod(zs+d+n,n);
        for ys=0:n-1
            yd = mod(ys+d+n,n);
            for xs=0:n-1
                xd = mod(xs+d+n,n);
                y(xd+1,yd+1,zd+1) = t(xs+1,ys+1,zs+1);
            end
        end
    end
    y = vec(y);
end
