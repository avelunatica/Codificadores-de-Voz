function [syn d]=syslpc_vqq(s, Lframe, p, VQ, par_string)
% Parámetros adicionales:
% s: señal de voz
% Lframe: longitud de cada trama (segmento de voz)
% p: orden de predicción
% VQ: matriz con los centroides por filas del VQ a aplicar.
% par_string: indica los coeficientes a cuantificar (’LPC’,’RC’,’LSF’ o ’LAR’)

% syn: señal reconstruída
% d: secuencia d[n] (error de predicción) e(n) = s(n) − sb(n) 
start_slice=1;
ntramas = floor(length(s)/Lframe);
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

for n=1:ntramas %utilizar floor para descartar el último segmento
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
end

start_slice=1;
stop_slice=Lframe;

if strcmp ('LPC', par_string)
        LPC=LPC(:,2:end);
        [LPCVQ,dist]=VQuantize(LPC,VQ);
        LPCVQ= [1 LPCVQ];
        LPC= [1 LPC]
            for n=1:ntramas %utilizar floor para descartar el último segmento
            frame=xn(start_slice:stop_slice);
            start_slice=start_slice+Lframe;
            stop_slice=stop_slice+Lframe;   
            [di z1]=filter(LPCVQ(i,:),1,frame,z1);
            [syni z2]=filter(1, LPCVQ(i,:), di,z2);
            d=[d di];
            syn=[syn syni];
            end
             if not((ntramas*Lframe)==length(xn))
            frame=xn(stop_slice:length(xn));
            [di,z1] = filter(LPCVQ(i,:),1,frame,z1);
            [syni,z2] = filter(1,LPCVQ(i,:),frame,z2);
             d=[d di];
             syn=[syn syni];
            end
    
            
elseif strcmp (par_string, 'RC')
        [RCVQ,dist]=VQuantize(RC,VQ);
        LPCVQ=rc2poly(RCVQ)
            for n=1:ntramas %utilizar floor para descartar el último segmento
            frame=xn(start_slice:stop_slice);
            start_slice=start_slice+Lframe;
            stop_slice=stop_slice+Lframe;   
            [di z1]=filter(LPCVQ(i,:),1,frame,z1);
            [syni z2]=filter(1, LPCVQ(i,:), di,z2);
            d=[d di];
            syn=[syn syni];
            end   
            if not((ntramas*Lframe)==length(xn))
            frame=xn(stop_slice:length(xn));
            [di,z1] = filter(LPCVQ(i,:),1,frame,z1);
            [syni,z2] = filter(1,LPCVQ(i,:),frame,z2);
             d=[d di];
             syn=[syn syni];
            end
    
        
elseif strcmp (par_string,'LAR')
         [LARVQ,dist]=VQuantize(LAR,VQ);
         LPCVQ=rc2poly(lar2rc(LARVQ));
             for n=1:ntramas %utilizar floor para descartar el último segmento
            frame=xn(start_slice:stop_slice);
            start_slice=start_slice+Lframe;
            stop_slice=stop_slice+Lframe;   
            [di z1]=filter(LPCVQ(i,:),1,frame,z1);
            [syni z2]=filter(1, LPCVQ(i,:), di,z2);
            d=[d di];
            syn=[syn syni];
             end  
            if not((ntramas*Lframe)==length(xn))
            frame=xn(stop_slice:length(xn));
            [di,z1] = filter(LPCVQ(i,:),1,frame,z1);
            [syni,z2] = filter(1,LPCVQ(i,:),frame,z2);
             d=[d di];
             syn=[syn syni];
            end
    
        
elseif strcmp (par_string,'LSF')
         [LSFVQ,dist]=VQuantize(LSF,VQ);
         LPCVQ=lsf2poly(LSFVQ);
            for n=1:ntramas %utilizar floor para descartar el último segmento
            frame=xn(start_slice:stop_slice);
            start_slice=start_slice+Lframe;
            stop_slice=stop_slice+Lframe;   
            [di z1]=filter(LPCVQ(i,:),1,frame,z1);
            [syni z2]=filter(1, LPCVQ(i,:), di,z2);
            d=[d di];
            syn=[syn syni];
            end 
            if not((ntramas*Lframe)==length(xn))
            frame=xn(stop_slice:length(xn));
            [di,z1] = filter(LPCVQ(i,:),1,frame,z1);
            [syni,z2] = filter(1,LPCVQ(i,:),frame,z2);
             d=[d di];
             syn=[syn syni];
            end
    
end
