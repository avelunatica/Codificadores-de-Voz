function [VQ Dist no_asig] = Kmeans(data, nbits, VQini, umbral, display)

% VQ: Matriz con los centroides resultantes por filas. 
% vDist: vector con la distorsión resultante en cada iteración 
% del algoritmo 
% no_asig: número de centroides no asignados (contienen NaN)
%
% data : matriz con vectores de entrenamiento por fila. 
%		       Dimensión VQ = número de columnas de data.
% nbits: numero de bits del VQ deseado.
% VQini: centroides iniciales por filas. Puede ser una matriz vacía.
% umbral: umbral para la variación relativa de la distorsión (MSE) 
%         entre iteraciones utilizado como criterio de parada (ej. 1e-2).
% display: si vale 1 y la dimensión de los vectores de entrenamiento es 2,
%          muestra alguna gráfica ilustrativa.
% Posibles colores en gráficas: vectores de entrenamiento (verde), centroides
%          en cada iteración (rojo), centroides resultantes (azul).
%
% Ejemplo:  load testdata; [VQ,MSE]=Kmeans(testdata,4,[],1e-3,1);



sdata=size(data);
N=sdata(1);
L=sdata(2);
M=2^nbits; 

MSE=[];

if(N<M) error('Número de datos insuficiente');
end

if(display && (L==2)) display=1;
else display=0;
end


%Inicialización aleatoria si la inicialización del VQ no está especificada.
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
   
    %Asignación de vectores a celdas
    parfor i=1:N
       
        difM=repmat(data(i,:),M,1)-VQ;
        distM=vecnorm(difM,2,2);
        
        [mindist minind]=min(distM);
        celda(i)=minind(1);

    end;

    VQ=[];
    % Cálculo del centroide de cada celda.
    for j=1:M
        dat=data(find(celda==j),:);
        VQ(j,:)=mean(dat,1); % vector de NaN si "dat" vacío.
        if(display) plot(VQ(j,1),VQ(j,2),'r.')
        end
    end
    
    %Número de centroides no asignados
    no_asig=length(find(isnan(VQ(:,1))==1));
    
    
    distorsion_prev=distorsion;
    
    distorsion=0;
    
    %Cálculo del MSE (medida de distorsion utilizada)
   
    dif=data-VQ(celda,:);
    distorsion=sum(sum(dif.^2))/N;
    
    MSE(it)=distorsion;
    
    pdistorsion=(distorsion-distorsion_prev)/distorsion;

end

%Normalización por la dimensión del vector
Dist=MSE/L;

if(display)
    plot(VQ(:,1),VQ(:,2),'b.')
    hold off;
end

        
    
  


