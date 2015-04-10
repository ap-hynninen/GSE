function [y] = Gk(qx, qy, qz);
    q = sqrt(qx.^2 + qy.^2 + qz.^2);
    y = q.^(-2) + q.^(-4).*(qx.^4 + qy.^4 + qz.^4)/12;
end

