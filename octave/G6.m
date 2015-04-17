function [y] = G6(qx, qy, qz);

    a = 97/120;
    b = 24/120;
    c = -1/120;
    y = 0.5./( (1-cos(qx))./(a + b*cos(qx) + c*cos(2*qx)) + ...
               (1-cos(qy))./(a + b*cos(qy) + c*cos(2*qy)) + ...
               (1-cos(qz))./(a + b*cos(qz) + c*cos(2*qz)));
    
end
