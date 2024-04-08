function [y,e] = qmidriser(x, xsc, n)
% Resumen de la función qmidriser
% x : vector de muestras de entrada (longitud arbitraria)
% xsc : valor de sobrecarga del cuantificador
% n : número de bits
%
% y : Vector con valores cuantificados
% e : Error de cuantificación, x-y 

delta=2*xsc/2^n;
L=2^n;
deltamed=delta/2;
U=delta*((L/2)-1)+deltamed; %xsc-deltamed


for n=1:length(x)
    k=floor(abs(x(n))/delta);
    if abs(x(n))>=xsc
        y(n)=sign(x(n))*U;
    elseif x(n)==0 %Preguntar
        y(n)=deltamed;
    else
       y(n)=sign(x(n))*delta*(k+0.5);
    end 
    e(n)=x(n)-y(n);
end



