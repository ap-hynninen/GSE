function res=G(r, sigma)

    res = exp(-r.^2/(2*sigma^2))/sigma; %*(2*pi*sigma^2)^(-3/2);
    
end
