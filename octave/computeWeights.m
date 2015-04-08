function [w dw] = computeWeights(sigma, dx, nx);

    sigmasq = sigma*sigma;
    w = zeros(2*nx+1,1);
    dw = zeros(2*nx+1,1);

    for i=-nx:nx
        w(i+nx+1) = exp(-(dx+i)*(dx+i)/(2*sigmasq))/sqrt(2*pi*sigmasq);
    end

    qx = 0;
    dqx = 0;
    for i=-nx:nx
        qx = qx + w(i+nx+1);
        dqx = dqx - w(i+nx+1)*(dx+i)/sigmasq;
    end
    
    for i=-nx:nx
        dw(i+nx+1) = w(i+nx+1)*((-dx-i)/sigmasq - 0*dqx/qx);
    end

end
