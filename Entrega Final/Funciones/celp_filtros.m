function [sh,B,G,AK,Tv,indv, bitsPorMuestra] = celpFiltrosAlternativos(s,Ltrama,Lsubtrama,p)

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
z_w=[];
z_inv_W_iniciales=[];
z_Hf_iniciales=[];


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

EntradaNula = zeros(1,Lsubtrama);

for i = 1:floor(length(s)/Lframe) % Para cada trama
    
    frame=s(start_slice_frame:stop_slice_frame);
    Rcorr=xcorr(frame,frame);
    Rcorr=Rcorr(length(frame):length(frame)+p);
    [ak,ep,rc] =levinson(Rcorr,p);
    AK=[AK ak];
    akPond = ak.*(0.8.^(0:p));
    
    for n = 1:4 % Para cada subtrama
        
        subtrama_actual = frame(Lsubframe*(n-1)+1:Lsubframe*n);
        
        % Paso la subtrama actual por el filtro W(Z) y actualizo las condiciones para la siguiente iteración   
        [u_n, z_w] = filter(ak,akPond,subtrama_actual, z_w);
        
        
        [y3, z_Hf_finales] = filter(1,akPond,EntradaNula,z_Hf_iniciales);
        u3 = u_n - y3;
        
        % Parte adaptativa -----------------------------------
        for t = 1 : length(libraria_adaptativa) - Lsubtrama + 1 
       
            d20 = libraria_adaptativa(t:t + Lsubtrama - 1);                 

            y20 = filter(1,akPond,d20); %LPC síntese
            b=(u_n*y20.')/(y20*y20'+ eps);  %tendo en conta n=1:4 
            y2=y20*b;
            u0 = u3 - y2; % en la cuatro era subtrama actual (sin filtrar) - y2
            Et=sum(u0.*u0);

            % variables: gardoas pero non sei se algunha non a necesitarei
            Y20(t,:)=y20;
            Y2(t,:)=y2; 
            ET(t) = Et; % Teño que buscar o erro minimo 
            Bposibles(t) = b; % Sabendo o erro minimo, consigo a ganancia
            U(t,:) =  u0;   
        
        end %--------------------------------------------------

        [valor_ET, Indice_t]=min(ET);
        b_que_minimiza = Bposibles(Indice_t);
        B = [B  b_que_minimiza]; % vector coas ganancias que minimizan a enerxia do erro en cada subtrama
        T = 3*length(Lsubframe) - Indice_t + 1;
        Tv=[Tv T]; %vector candidato (ou indice t)
        
        d2 = libraria_adaptativa(Indice_t:Indice_t+Lsubtrama-1)*b_que_minimiza;
        [y2, z_adap] = filter(1,akPond,d2);
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
            
            y10 = filter(1,akPond,v(j,:)); %LPC síntese
            g=(U0*y10.')/(y10*y10'+eps);
            y1=y10*g;
            uf = U0 - y1;
            El=sum(uf.*uf);

            % Variables: gardoas pero non sei se algunha non a necesitarei
            Gposibles(j)=g;
            UF(j,:)=uf;
            EL(j) = El;
            % matriz_z_libraria_estocastica(j,:) = z_est_final';
            Y1(j,:)=y1;
            Y10(j,:)=y10;
            
        end %------------------------------------------------------------------------------------------------

         [valor_EL, Indice_indv]=min(EL);
         g_que_minimiza=Gposibles(Indice_indv); %es el error con la ganancia
         G=[G g_que_minimiza];
         indv = [indv  Indice_indv];
         
         d1 = v(Indice_indv,:)*g_que_minimiza;
         [y1,z_est]=filter(1,akPond,d1);

         d = d2 + d1;  
         
         % Actualizo la librería adaptativa
         libraria_adaptativa_pasada = libraria_adaptativa (Lsubframe+1:3*Lsubframe);
         libraria_adaptativa = [libraria_adaptativa_pasada d];
         
         z_Hf_iniciales = z_Hf_finales + z_est + z_adap;

         % sh: señal reconstruida
         shi = y2 + y1 + y3;
         % Ahora para reconstruir la señal de voz será necesario pasar y2(t)+y1(l) por 1/W(Z)
         [shiFiltrada, z_inv_W_finales] = filter(akPond,ak,shi, z_inv_W_iniciales);
         z_inv_W_iniciales = z_inv_W_finales; % actualizo condiciones del filtro para la sig. iteración
         sh =[sh shiFiltrada];
         
         if(i<2) 
             figure   
             subplot(511), plot(0:Lsubtrama-1, subtrama_actual);xlabel('n');ylabel('s'),grid,
             title("subtrama " + n + " de la trama " + i)
             subplot(512), plot(0:Lsubtrama-1, d1),xlabel('n');ylabel('v*g'),grid,
             subplot(513), plot(0:Lsubtrama-1, Y1(Indice_indv,:)),xlabel('n');ylabel('y1'),grid,
             subplot(514), plot(0:Lsubtrama-1, UF(Indice_indv,:)),xlabel('n');ylabel('uf'),grid,
             subplot(514), plot(0:Lsubtrama-1, shiFiltrada),xlabel('n');ylabel('sh'),grid,
         end
         
    end
    
    % actualizo os slices
    start_slice_frame = start_slice_frame+Lframe;
    stop_slice_frame = stop_slice_frame+Lframe;
    
end

end

