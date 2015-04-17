function [y] = G3(qx, qy, qz);

    a = 5/6;
    b = 1;
    y = (1/2)*(b+(1-b)/3*cos(qx)+(1-b)/3*cos(qy)+(1-b)/3*cos(qz)).* ...
        (a+(1-a)/3*(cos(qx)+cos(qy)+cos(qz)))./(3-(cos(qx)+cos(qy)+cos(qz)));

end
