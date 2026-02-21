function f = ITR1D(a,nu)
n = length(a);
aux = (n/pi)^(1/2) * a;
aux = dct(aux);
aux = aux .* sqrt(nu);
f = dct(aux,Type=3);