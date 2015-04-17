function [y] = G7(qx,qy,qz);

    a = 7033/3780;
    b = -1741/2520;
    c = -223/1260;
    d = 11/1512;
    y = 1./( (13/2-8*cos(qx)+3/2*cos(2*qx))./(a+b*cos(qx)+c*cos(2*qx)+d*cos(3*qx)) + ...
             (13/2-8*cos(qy)+3/2*cos(2*qy))./(a+b*cos(qy)+c*cos(2*qy)+d*cos(3*qy)) + ...
             (13/2-8*cos(qz)+3/2*cos(2*qz))./(a+b*cos(qz)+c*cos(2*qz)+d*cos(3*qz)));
    
end
