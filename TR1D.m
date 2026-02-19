function a = TR1D(f,nu)
n = length(f);
aux = dct(f);
aux = aux ./ sqrt(nu);
aux = dct(aux,Type=3);
a = (pi/n)^2 * aux;