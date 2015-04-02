function res=Usigma(r, N, b, sigma)
    res = 0;
    for j=0:N-1
        sigmaj = sigma*b^j;
        res = res + G(r, sigmaj);
    end
    res = res*2*log(b);
end
