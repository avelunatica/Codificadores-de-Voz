 function [VQ Dist] = bsVQ(data, nbits, epsilon, umbral, display)
% VQ: Matriz con los centroides resultantes por filas.
% vDist: vector con la distorsión (MSE) en cada iteración del algoritmo Kmeans
% sobre la última partición del binary splitting (última ejecución Kmeans).
%
% data : matriz con vector de datos por fila.
% nbits: numero de bits del VQ deseado.
% epsilon: parametro para dividir cada centroide en 2 (ej: 0.025).
% threshold: umbral de parada para K-means (ej:0.01)
% display: si vale 1 y la dimensión es 2 entonces representa gráficas ilustrativas.
%
% Ejemplo: load traindata; [VQ, vDist]=bsVQ(traindata,4,0.025,1e-3,1);

VQini = [];
Dist = [];

%Bucle duplicar tamaño VQ y estimación centroides k-means

for i= 0:nbits %la primera es con VQini vacía porque va desde 0 (todo el conjunto)
    %Obtén un VQ con un único centroide, el correspondiente a todo el conjunto de entrenamiento.
    [VQ vDist] = Kmeans(data, i, VQini, umbral, display);
    %Duplica el tamaño del VQ, creando a partir de cada centroide actual, y , los vectores.
    Ymas = VQ * (1 + epsilon);
    Ymenos = VQ * (1 - epsilon);
    %Matriz con centroides
    VQini=[Ymas;Ymenos];
    
    Dist=[Dist vDist];
    
end

if display == 1
    figure;

% plot(vDist);
% xlabel ('Iteración'); ylabel ('Distorsión');
% title 'Distorsión por cada iteración';

end


end
    
    
    
    
