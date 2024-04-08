function [VQ Dist no_asig] = Kmeans(data, nbits, VQini, umbral, display)

% VQ: Matriz con los centroides resultantes por filas. 
% vDist: vector con la distorsi�n resultante en cada iteraci�n 
% del algoritmo 
% no_asig: n�mero de centroides no asignados (contienen NaN)
%
% data : matriz con vectores de entrenamiento por fila. 
%		       Dimensi�n VQ = n�mero de columnas de data.
% nbits: numero de bits del VQ deseado.
% VQini: centroides iniciales por filas. Puede ser una matriz vac�a.
% umbral: umbral para la variaci�n relativa de la distorsi�n (MSE) 
%         entre iteraciones utilizado como criterio de parada (ej. 1e-2).
% display: si vale 1 y la dimensi�n de los vectores de entrenamiento es 2,
%          muestra alguna gr�fica ilustrativa.
% Posibles colores en gr�ficas: vectores de entrenamiento (verde), centroides
%          en cada iteraci�n (rojo), centroides resultantes (azul).
%
% Ejemplo:  load testdata; [VQ,MSE]=Kmeans(testdata,4,[],1e-3,1);



sdata=size(data);
N=sdata(1);
L=sdata(2);
M=2^nbits; 

MSE=[];

if(N<M) error('N�mero de datos insuficiente');
end

if(display && (L==2)) display=1;
else display=0;
end


%Inicializaci�n aleatoria si la inicializaci�n del VQ no est� especificada.
if(length(VQini)<M)
    rind=randperm(N);
    VQ=data(rind(1:M),:);
else VQ=VQini;
end

distM=[];
celda=zeros(N,1);

pdistorsion=1;
distorsion=1;

it=0;

if (display)
    figure(1);
    clf;
    plot(data(:,1),data(:,2),'g.')
    hold on;
end


%Comienzan las iteraciones
while(abs(pdistorsion)> umbral) %Criterio de parada

    it=it+1;
   
    %Asignaci�n de vectores a celdas
    parfor i=1:N
       
        difM=repmat(data(i,:),M,1)-VQ;
        distM=vecnorm(difM,2,2);
        
        [mindist minind]=min(distM);
        celda(i)=minind(1);

    end;

    VQ=[];
    % C�lculo del centroide de cada celda.
    for j=1:M
        dat=data(find(celda==j),:);
        VQ(j,:)=mean(dat,1); % vector de NaN si "dat" vac�o.
        if(display) plot(VQ(j,1),VQ(j,2),'r.')
        end
    end
    
    %N�mero de centroides no asignados
    no_asig=length(find(isnan(VQ(:,1))==1));
    
    
    distorsion_prev=distorsion;
    
    distorsion=0;
    
    %C�lculo del MSE (medida de distorsion utilizada)
   
    dif=data-VQ(celda,:);
    distorsion=sum(sum(dif.^2))/N;
    
    MSE(it)=distorsion;
    
    pdistorsion=(distorsion-distorsion_prev)/distorsion;

end

%Normalizaci�n por la dimensi�n del vector
Dist=MSE/L;

if(display)
    plot(VQ(:,1),VQ(:,2),'b.')
    hold off;
end

        
    
  


