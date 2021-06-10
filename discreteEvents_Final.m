clc 
clear 
close all


DIFS = 10e-3; % Espacio de intertrama distribuido
SIFS = 5e-3; % Espacio de intertrama corto
durRTS = 11e-3; %paquetes de control
durCTS = 11e-3; %paquetes de control
durACK = 11e-3; %paquetes de control
durDATA = 43e-3; %paquetes de datos
sigma = 1e-3; %duracion de cada miniranura
I = 7; % Grado 
k = 15; %Tamaño del buffer
Rd = 18; %Numero de ranuras dormir (E rara)
N = [5,10,15,20]; %Numero de nodos por grado
n = 2; %indice de nodos
W = [16,32,64,128,256]; %Maximo numero de miniranuras de tiempo
w = 2; %indice de miniranuras
lanmbda = [0.0005,0.001, 0.005, 0.03]; % Tasa de generacion de paquetes/seg
l = 1; %indice lanmbda

% Duracion de la ranura
T = sigma.*W(w)+ DIFS + 3*SIFS + durRTS+ durCTS+ durDATA + durACK;

%Ciclo de trabajo 
Tc = (2+Rd).*T;

%Matriz del sistema
Nodos = zeros(N(n),I);

%Estados de los buffers
hashMap = containers.Map('KeyType','char','ValueType','any');

collisions = zeros(1, 7);
fullbuffer = zeros(1, 7);
nodeTimes = zeros(1, 7);
countersuccess = zeros(1, 7);
counterontransit = zeros(1, 7);
piepackages = zeros(1, 3);
transpergrade = zeros(1, 7);

for f = 1: N(n)
   for c = 1: I
       hashMap(strcat(int2str(f),int2str(c))) = zeros(1,15);
   end
end

%Garantiza que hay almenos un arribo al inicar la simulacion 
ta = -1;
pg=0;
tsim = 0;
final = 0;
counter = 0;

contendientes = zeros(1, N(n));
numPaquete = 0;

for t = 1:500000*Tc
     %Generacion de paquete
    if(ta < tsim)
        lanmbda2 = lanmbda(l)*N(n)*I;
        U = (1e6*rand())/1e6;
        nuevot = -(1/lanmbda2)*log(1-U);
        nodo = randi([1 N(n)],1);
        grado = randi([1 I],1);
        
        if(Nodos(nodo,grado) < k)
            
            Nodos(nodo,grado) = Nodos(nodo,grado) + 1;
            filas = int2str(nodo);
            columnas = int2str(grado);
            llave = strcat(filas,columnas);
            Aux = hashMap(llave);
            index = find(Aux==0);
            if ~isempty(index)
                pg = pg+1;
                Aux(index(1)) = pg;
                hashMap(llave) = Aux;
                Aux = 0;
            end
            
            %Conteo paquetes
            %Pa=[Id,nodo,grado,ta,tsink]
            Pa(pg,:)=zeros(1,5);
            Pa(pg,1) = pg;
            Pa(pg,2) = nodo;
            Pa(pg,3) = grado;
            Pa(pg,4) = ta;
            
        end
        ta = tsim+nuevot;
    end       
    
   
     if mod(t, Tc) == 0  %Transmitir cada ciclo de trabajo.
        for grado = I:-1:1 %Proceso de transmisión
            for row = 1: N(n)
                if Nodos(row,grado) ~= 0
                    contendientes(row) = randi([1,W(w)],1);
                else
                    contendientes(row) = NaN;
                end
            end
            %Colisiones

            nodoGanador = min(contendientes);
            colisiones = find(contendientes == nodoGanador);
            if length(colisiones) > 1
                %Hay colision
                for index_colision = 1: length(colisiones)
                    %Eliminamos los paquetes que colisionaron
                    %Actualizar HASH
                     llaveHash = strcat(int2str(colisiones(index_colision)), int2str(grado));
                     aux =  hashMap(llaveHash);
                     Pa(aux(1), 5) = -1; %Actualizar Pa
                     aux(1) = [];
                     aux(15) = 0;
                     hashMap(llaveHash) = aux;
                    %Actualizar NODOS
                    Nodos(colisiones(index_colision), grado) = Nodos(colisiones(index_colision), grado) - 1; 
                end
            else
%                 %Podemos transmitir
%                 %Revisar si hay algo que transmitir al menos un paquete en el grado.
                if ~isnan(nodoGanador) 
                    %Restar al buffer del hashMap paquete del nodo (Grado I)
                    indexNodo = find(contendientes == nodoGanador);
                    llaveHash = strcat(int2str(indexNodo), int2str(grado));
                    aux =  hashMap(llaveHash);
                    numPaquete = aux(1);
                    aux(1) = [];
                    aux(15) = 0;
                    hashMap(llaveHash) = aux;
                    Nodos(indexNodo, grado) = Nodos(indexNodo, grado) - 1;%Actualizar matriz Nodos (grado I)

                    if (grado - 1 )~= 0 %Sumar al buffer del hashMap paquete del nodo (Grado I -1)
                        llaveHash = strcat(int2str(indexNodo), int2str(grado-1));
                        aux =  hashMap(llaveHash);
                        aux(find(aux == 0,1)) = numPaquete;
                        if aux(15) == 0
                            transpergrade(grado) = transpergrade(grado) + 1;
                            hashMap(llaveHash) =  aux;
                            Nodos(indexNodo, grado - 1) = Nodos(indexNodo, grado -1) + 1; %Actualizar matriz Nodos (grado I-1)
                        else
                            Pa(numPaquete, 5) = -2;
                        end
                    else
                        Pa(numPaquete, 5) = tsim;
                        transpergrade(grado) = transpergrade(grado) + 1;
                    end
                end
            end
        end
        
        
        counter = counter +1;
     end
     %Termina Transmitir cada ciclo de trabajo.
     tsim = tsim + T;
end

%Estadística de paquetes
for pck = 1: pg
    if Pa(pck, 5) == -1
        collisions(1, Pa(pck, 3)) = collisions(1, Pa(pck, 3)) + 1;
    elseif Pa(pck, 5) == -2
        fullbuffer(1, Pa(pck, 3)) = fullbuffer(1, Pa(pck, 3)) + 1;
    end
end

% aquí
Pa(1,4) = 0;
for pck = 1:pg
   if Pa(pck,5) ~= -1 && Pa(pck, 5) ~= -2 && Pa(pck, 5) ~= 0
      nodeTimes(1, Pa(pck, 3)) =  nodeTimes(1, Pa(pck, 3)) + (Pa(pck, 5) - Pa(pck, 4));
   end
end
    
for pck = 1:pg
   if Pa(pck,5) ~= -1 && Pa(pck, 5) ~= -2 
      switch Pa(pck, 3)
          case 1
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
          case 2
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
          case 3
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
          case 4
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
          case 5
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
          case 6
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
          case 7
              if Pa(pck,5) == 0
                  counterontransit(1, Pa(pck, 3)) = counterontransit(1, Pa(pck, 3)) + 1;
              end
              countersuccess(1, Pa(pck, 3)) = countersuccess(1, Pa(pck, 3)) + 1;
      end
   end
end
   
piepackages(1, 1) = sum(collisions);
piepackages(1, 2) = sum(fullbuffer);
piepackages(1, 3) = sum(countersuccess);

totaltrans = sum(transpergrade);
ptx = transpergrade ./ totaltrans;

endtoend = sum(nodeTimes) / tsim;

figure(1)
stem(collisions, 'LineWidth',2)
xlim([0 8])
title('Lost packages per grade(Collisions)', 'FontSize', 30)
ylabel('# lost packages', 'FontSize', 30)
xlabel('grade', 'FontSize', 30)
grid on

figure(2)
stem(fullbuffer, 'LineWidth',2)
xlim([0 8])
title('Lost packages per grade (Full buffer)', 'FontSize', 30)
ylabel('# lost packages', 'FontSize', 30)
xlabel('grade', 'FontSize', 30)
grid on

figure(3)
plot(nodeTimes ./ tsim, 'LineWidth',2)
xlim([0 8])
title('Delay source-to-sink', 'FontSize', 30)
ylabel('Time[s]', 'FontSize', 30)
xlabel('grade', 'FontSize', 30)
grid on

figure(4)
stem(transpergrade, 'LineWidth',2)
xlim([0 8])
title('Throughput', 'FontSize', 30)
ylabel('[Packages/Tc]', 'FontSize', 30)
xlabel('grade', 'FontSize', 30)
grid on

figure(5)
labels = {'Collisions', 'Full Buffer', 'Successfull Transmission'};
pie(piepackages, '%.2f%%')
title('Generated packages', 'FontSize', 30)
% Create legend
lgd = legend(labels);

figure(6)
labels = {'Grade 1', 'Grade 2', 'Grade 3', 'Grade 4', 'Grade 5', 'Grade 6', 'Grade 7'};
pie(collisions, '%.2f%%')
title('Lost packages by collisions', 'FontSize', 30)
% Create legend
lgd = legend(labels);

figure(7)
labels = {'Grade 1', 'Grade 2', 'Grade 3', 'Grade 4', 'Grade 5', 'Grade 6', 'Grade 7'};
pie(fullbuffer, '%.2f%%')
title('Lost packages by full buffer', 'FontSize', 30)
% Create legend
lgd = legend(labels);

figure(8)
labels = {'Grade 1', 'Grade 2', 'Grade 3', 'Grade 4', 'Grade 5', 'Grade 6', 'Grade 7'};
pie(fullbuffer + collisions, '%.2f%%')
title('Lost packages', 'FontSize', 30)
% Create legend
lgd = legend(labels);