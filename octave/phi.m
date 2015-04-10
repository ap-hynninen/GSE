function [y] = phi(qx, qy, qz);

    y = (1/2 + 1/6*(cos(qx) + cos(qy) + cos(qz)))./...
    (4 - 2/3*(cos(qx) + cos(qy) + cos(qz)) - 2/3*(cos(qx).*cos(qy) + cos(qx).*cos(qz) + cos(qy).*cos(qz)));
    
end
