function [syn d]=syslpc(s, Lframe, p)
% Parámetros adicionales:
% s: señal de voz
% Lframe: longitud de cada trama (segmento de voz)
% p: orden de predicción

% syn: señal reconstruída
% d: secuencia d[n] (error de predicción) e(n) = s(n) − sb(n) 

start_slice=1;
stop_slice=Lframe;
z1=[];
z2=[];
LPC=[];
Ep=[];
RC=[];
d=[];
syn=[];
xn=s.'; 
i=1;

for n=1:floor(length(xn)/Lframe) %utilizar floor para descartar el último segmento
    frame=xn(start_slice:stop_slice);
    
    Rcorr=xcorr(frame,frame); 
    Rcorr=Rcorr(length(frame):length(frame)+p);
    [ak,ep,rc] =levinson(Rcorr,p);
    
    LPC(i, :) = ak;
    Ep(i, :)=ep; 
    RC(i,:) = rc;
    
    LAR(i,:) = rc2lar(rc);
    LSF(i,:) = poly2lsf(ak);

    start_slice=start_slice+Lframe;
    stop_slice=stop_slice+Lframe;   

    [di z1]=filter(LPC(i,:),1,frame,z1);
    [syni z2]=filter(1, LPC(i,:), di,z2);

    d=[d di];
    syn=[syn syni];
end



end