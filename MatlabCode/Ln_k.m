function [poly] = Ln_k(n, k, x)
    poly = 0;
    for m = 0:n
        poly = poly + (((-1)^m).*factorial(n+k)./(factorial(n-m).*factorial(k+m).*factorial(m))).*x.^m;
    end
end
