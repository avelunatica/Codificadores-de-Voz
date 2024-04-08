function [SNRseg,SNRm,m] = SNRS(x,xq,L)
% Resumen de la función SNRq
% x= señal original
% xq= señal cuantificada
%L= longitud del segmento
%SNRm = relación señal ruido por tramas (vector)
%SNRseg= relación señal a ruido segmental
%m=vector de referencia en tiempo discreto para S
start_slice=1;
stop_slice=L;
x=x.';
for n=1:floor(length(x)/L) %utilizar floor para descartar el último segmento
    xi=x(start_slice:stop_slice);
    xiq=xq(start_slice:stop_slice);
    SNRm(n)=10*log10(sum(xi.^2)./sum((xi-xiq).^2));
    m(n)=start_slice;
    start_slice=start_slice+L;
    stop_slice=stop_slice+L;
end
SNRseg=mean(SNRm);
end

