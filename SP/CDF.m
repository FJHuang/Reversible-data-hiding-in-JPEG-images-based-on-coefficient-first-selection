function y = CDF(lambda,x)
    y = 0.5+0.5*sign(x)*(1-exp(-lambda*abs(x)));
end

