function [y,e] = qmidtread(x, xsc, n)
% Resumen de la función qmidtread
% x : vector de muestras de entrada (longitud arbitraria)
% xsc : valor de sobrecarga del cuantificador
% n : número de bits
%
% y : Vector con valores cuantificados
% e : Error de cuantificación, x-y 

delta=2*xsc/(2^n-1);
L=2^n-1; %niveles de cuantificacion
U=delta*floor(L/2);

for n=1:length(x)
    k=round(abs(x(n))/delta);
    if abs(x(n))>=xsc
        y(n)=sign(x(n))*U;
    else
        y(n)=sign(x(n))*delta*k;
    end 
    e(n)=x(n)-y(n); 
end

