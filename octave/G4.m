function [y] = G4(qx,qy,qz);

%a = 5/6;
    a = 0.763007592378452;
    b = 1.07032574095488;
    %y = (b+(1-b)/3*cos(qx)+(1-b)/3*cos(qy)+(1-b)/3*cos(qz))./(2*( (1-cos(qx))./(a+(1-a)*cos(qx)) + ...
    %                                                  (1-cos(qy))./(a+(1-a)*cos(qy)) + ...
    %                                                  (1-cos(qz))./(a+(1-a)*cos(qz)) ));

    y = (1/2)*(b+(1-b)/3*cos(qx)+(1-b)/3*cos(qy)+(1-b)/3*cos(qz)).* ...
        (a+(1-a)/3*(cos(qx)+cos(qy)+cos(qz)))./(3-(cos(qx)+cos(qy)+cos(qz)));

end
