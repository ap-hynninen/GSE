function [y] = G5(qx, qy, qz);

    a=13/30;
    b=16/30;
    c=1/30;
    y = 2./( (1-cos(2*qx))./(a+b*cos(qx)+c*cos(2*qx)) +...
             (1-cos(2*qy))./(a+b*cos(qy)+c*cos(2*qy)) +...
             (1-cos(2*qz))./(a+b*cos(qz)+c*cos(2*qz)) );
    
end
