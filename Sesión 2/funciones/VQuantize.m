function [y, Dist] = VQuantize(x,VQ)

% y: matriz con los vectores resultantes de la cuantificación
% Dist: Error cuadrático medio normalizado por la dimensión del vector.
%
% x : matriz con vectores a cuantificar por fila.
% VQ: Matriz con los centroides del VQ a utilizar por filas.

sx=size(x);
sVQ=size(VQ);

y=[];
MSE=0;
M=sVQ(1);

parfor i=1:sx(1)
    
    difM=repmat(x(i,:),M,1)-VQ;
    distM=vecnorm(difM,2,2);
        
    [mindist minind]=min(distM);
    y=[y; VQ(minind(1),:)];
    MSE=MSE+mindist.^2; 
end
MSE=MSE/sx(1);
Dist=MSE/sx(2);





end

