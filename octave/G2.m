function [y] = G2(qx, qy, qz);

    a = 5/6;
    b = 1/6;
    y = 0.5./( (1-cos(qx))./(a+b*cos(qx)) + ...
             (1-cos(qy))./(a+b*cos(qy)) + ...
             (1-cos(qz))./(a+b*cos(qz)) );
    
end
