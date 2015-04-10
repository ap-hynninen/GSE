function [y] = Ginv(qx, qy, qz);
    y = 6 - 2*cos(qx) - 2*cos(qy) - 2*cos(qz);
end
