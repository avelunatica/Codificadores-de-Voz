function SNRq = SNR(x,xq)
% Resumen de la función SNRq
% x= señal original
% xq= señal cuantificada
%SNRq= relación señal ruido
x=x.';
SNRq=10*log10(sum(x.^2)./sum((x-xq).^2));

end

