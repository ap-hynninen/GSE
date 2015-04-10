function [y] = G3(qx, qy, qz);

    y = 1./(4*sin(qx)./(5/6+1/6*cos(qx)) + 4*sin(qy)./(5/6+1/6*cos(qy)) ...
            + 4*sin(qz)./(5/6+1/6*cos(qz)));
end
