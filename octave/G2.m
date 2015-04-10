function [y] = G2(qx, qy, qz);

    y = 1./(12*( (1-cos(qx))./(5+cos(qx)) + (1-cos(qy))./(5+cos(qy)) ...
                 + (1-cos(qz))./(5+cos(qz)) ));
    
end
