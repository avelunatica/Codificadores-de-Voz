function [LPC, Ep, RC, LSF, LAR]=speech2lpc(s, p, window, wshift)

% Extracción de parámmetros LPC y equivalentes de una señal.

% s: señal de voz o sonora.
% p: orden de predicción
% window, wshift: ventana a utilizar y su desplazamiento en muestras.

% Salidas: parámetros correspondientes a cada trama por filas.
% LPC, Ep: LPCs (matriz) y energía del error de predicción (vector).
% RC: coeficientes de reflexión (matriz)
% LSF: Line Spectral Frequencies (matriz)
% LAR: Log Area Ratios (matriz)

w_start = 1;
w_end = length(window);

LPC = [];
Ep = [];
RC = [];
LSF = [];
LAR = [];

for i = 1: floor(length(s)/wshift)
    
frame = s(w_start:w_end);
frame = frame.*window;

Rcorr=xcorr(frame,frame);
Rcorr=Rcorr(length(frame):length(frame)+p);
[ak,ep,rc] =levinson(Rcorr,p);

LPC(i,:) = ak;
Ep(i,:) = ep;
RC(i,:) = rc;

LAR(i,:) = rc2lar(rc);
LSF(i,:) = poly2lsf(ak);

w_start = w_start + wshift;
w_end = w_end + wshift;
end 

end

