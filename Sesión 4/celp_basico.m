function [sh,B,G,AK,Tv,indv, bitsPorMuestra] = celp_basico(s,Ltrama,Lsubtrama,p)


Lframe = Ltrama;
Lsubframe = Lsubtrama;

% inicializo os slices para as tramas
start_slice_frame=1;
stop_slice_frame=Lframe;

%CALCULO No BITS
indBiblEst = 9 * (Ltrama/Lsubtrama); %una vez cada subtrama
indBiblAdapt = 7 * (Ltrama/Lsubtrama); %una vez cada subtrama
LPCbitsPorTrama = p*3; %3 bits por LPC, p totales en la trama
ganancias = 8 * (Ltrama/Lsubtrama) * 2; %8 bits por ganancia por la cantidad de subtramas, x2 porque hay 2 (est + adapt)

bitsPorTrama = LPCbitsPorTrama + ganancias + indBiblEst + indBiblAdapt;
bitsPorMuestra = bitsPorTrama / Ltrama;
 
% Inicializamos as variables de saida
sh=[];             
B=[];             
G=[];             
AK=[];            
Tv=[];             
indv=[];         

% Inicializamos a libraria estocastica
M = 512; N=Lsubtrama; rng(1); v = randn(M,N);
 
% Inicializamos a libraría adaptativa
libraria_adaptativa=zeros(1,3*Lsubtrama); %debemos ter en conta que agora o noso valor de subtrama pode ser de 40 ou 60.

% Inicializamos as condicións iniciais de cada filtro (cada un ten as súas)
z_adap_inicial=[];
z_est_inicial=[];

% Outras variables
u=[];
u0=[];
U=[];
UF=[];
B=[];
Bposibles=[];
G=[];
Gposibles=[];
AK=[];
Y1=[];
Y2=[];

for i = 1:floor(length(s)/Lframe) % Para cada trama
    
    frame=s(start_slice_frame:stop_slice_frame);
    Rcorr=xcorr(frame,frame);
    Rcorr=Rcorr(length(frame):length(frame)+p);
    [ak,ep,rc] =levinson(Rcorr,p);
    AK=[AK ak];
    
    for n = 1:4 % Para cada subtrama
        
        subtrama_actual = frame(Lsubframe*(n-1)+1:Lsubframe*n);
        
        % Parte adaptativa -----------------------------------
        for t = 1 : length(libraria_adaptativa) - Lsubtrama + 1 
       
            d20 = libraria_adaptativa(t:t + Lsubtrama - 1);                 

            [y20, z_adap_final] = filter(1,ak,d20, z_adap_inicial); %LPC síntese
            b=(subtrama_actual*y20.')/(y20*y20'+ eps);  %tendo en conta n=1:4 
            y2=y20*b;
            u0 = subtrama_actual - y2;
            Et=sum(u0.*u0);

            % variables: gardoas pero non sei se algunha non a necesitarei
            Y20(t,:)=y20;
            Y2(t,:)=y2; 
            ET(t) = Et; % Teño que buscar o erro minimo 
            Bposibles(t) = b; % Sabendo o erro minimo, consigo a ganancia
            matriz_z_libraria_adaptativa(t,:) = z_adap_final';
            U(t,:) =  u0;   
            % pause;
        
        end %--------------------------------------------------

        [valor_ET, Indice_t]=min(ET);
        b_que_minimiza = Bposibles(Indice_t);
        B = [B  b_que_minimiza]; % vector coas ganancias que minimizan a enerxia do erro en cada subtrama
        T = 3*length(Lsubframe) - Indice_t + 1;
        Tv=[Tv T]; %vector candidato (ou indice t)
        z_adap_inicial= matriz_z_libraria_adaptativa(Indice_t,:); 
        d2 = libraria_adaptativa(Indice_t:Indice_t+Lsubtrama-1)*b_que_minimiza;
        U0=U(Indice_t,:);
        
        if(i<2) 
             figure   
             subplot(511), plot(0:Lsubtrama-1, subtrama_actual);xlabel('n');ylabel('s'),grid,
             title("subtrama " + n + " de la trama " + i)
             subplot(512), plot(-3*Lsubtrama:-1, libraria_adaptativa),xlabel('n');ylabel('adapt.'),grid,
             subplot(513), plot(0:Lsubtrama-1, d2),xlabel('n');ylabel('d2o*b'),grid,
             subplot(514), plot(0:Lsubtrama-1, Y2(Indice_t,:)),xlabel('n');ylabel('y2'),grid,
             subplot(515), plot(0:Lsubtrama-1, U(Indice_t,:)),xlabel('n');ylabel('u0'),grid,
        end     

        % Parte estocástica --------------------------------------

        for j = 1:M % Seguindo as diapos debería ser unha l, pero non diferencio 'l' de '1'
            
            [y10, z_est_final_c] = filter(1,ak,v(j,:), z_est_inicial); %LPC síntese
            g=(U0*y10.')/(y10*y10'+eps);
            y1=y10*g;
            uf = U0 - y1;
            El=sum(uf.*uf);

            % Variables: gardoas pero non sei se algunha non a necesitarei
            Gposibles(j)=g;
            UF(j,:)=uf;
            EL(j) = El;
            z_est_final=z_est_final_c';
            matriz_z_libraria_estocastica(j,:) = z_est_final;
            Y1(j,:)=y1;
            Y10(j,:)=y10;
            
        end %------------------------------------------------------------------------------------------------

         [valor_EL, Indice_indv]=min(EL);
         g_que_minimiza=Gposibles(Indice_indv); %es el error con la ganancia
         G=[G g_que_minimiza];
         indv = [indv  Indice_indv];
         z_est_inicial = matriz_z_libraria_estocastica(Indice_indv,:);
         d1 = v(Indice_indv,:)*g_que_minimiza; 
         
         d = d2 + d1;  
         
         % Actualizo la librería adaptativa
         libraria_adaptativa_pasada = libraria_adaptativa (Lsubframe+1:3*Lsubframe);
         libraria_adaptativa = [libraria_adaptativa_pasada d];                    

         % sh: señal reconstruida
         shi = Y2(Indice_t,:) + Y1(Indice_indv,:);
         sh =[sh shi];
         
         if(i<2) 
             figure   
             subplot(511), plot(0:Lsubtrama-1, subtrama_actual);xlabel('n');ylabel('s'),grid,
             title("subtrama " + n + " de la trama " + i)
             subplot(512), plot(0:Lsubtrama-1, d1),xlabel('n');ylabel('v*g'),grid,
             subplot(513), plot(0:Lsubtrama-1, Y1(Indice_indv,:)),xlabel('n');ylabel('y1'),grid,
             subplot(514), plot(0:Lsubtrama-1, UF(Indice_indv,:)),xlabel('n');ylabel('uf'),grid,
             subplot(514), plot(0:Lsubtrama-1, shi),xlabel('n');ylabel('sh'),grid,
         end
         
    end
    
    % actualizo os slices
    start_slice_frame = start_slice_frame+Lframe;
    stop_slice_frame = stop_slice_frame+Lframe;
    
end
end







