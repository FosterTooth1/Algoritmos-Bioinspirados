clear, clc, close all;

% Solicitar al usuario el número de veces que se debe repetir el algoritmo
Num_iteraciones = input('Ingrese el número de veces que se repetirá el algoritmo: ');

%Creamos un arreglo donde se almacenaran los individuos con mejor aptitud (minimzacion) por cada iteracion
Resultados_Generales= zeros(1,Num_iteraciones);

for iteracion = 1:Num_iteraciones
    
    Num_pob=100;
    Num_Gen=200;
    Pm=0.15;
    Num_var=11;
    
    Distancias = [
        0 3091 927 1876 2704 94 2999 1641 3471 1838 3013;
        3091 0 2542 1681 375 2994 138 1442 389 1407 290;
        927 2542 0 1337 2169 930 2464 1100 2935 1168 2465;
        1876 1681 1337 0 1308 1778 1603 240 2075 163 1604;
        2704 375 2169 1308 0 2603 366 1069 767 1034 296;
        94 2994 930 1778 2603 0 2898 1543 3369 1740 2915;
        2999 138 2464 1603 366 2898 0 1364 531 1329 338;
        1641 1442 1100 240 1069 1543 1364 0 1836 201 1365;
        3471 389 2935 2075 767 3369 531 1836 0 1801 607;
        1838 1407 1168 163 1034 1740 1329 201 1801 0 1330;
        3013 290 2465 1604 296 2915 338 1365 607 1330 0
    ];
    
    %Poblacion inicial
    Pob=zeros(Num_pob, Num_var);
    for i=1:Num_pob
        Pob(i,:)= randperm(Num_var);
    end
    
    %Evaluación
    Aptitud_Pob=zeros(Num_pob, 1);
    for i=1:Num_pob
        Aptitud_Pob(i)= Costo(Pob(i,:), Distancias, Num_var);
    end
    
    % Inicializar mejor histórico
    [Mejor_Aptitud_Historico, idx] = min(Aptitud_Pob);
    Mejor_Individuo_Historico = Pob(idx, :);
    
    iter= 1;
    while iter <= Num_Gen
        % Cruzamiento y evaluación
        Hijos = zeros(Num_pob / 2, Num_var);
        Aptitud_Hijos = zeros(Num_pob / 2, 1);
        for i = 1:2:Num_pob
            % Generar hijo
            hijo = EdgeRecombination(Pob(i, :), Pob(i+1, :), Num_var);
            Hijos((i+1)/2, :) = hijo;
            % Evaluar hijo
            Aptitud_Hijos((i+1)/2) = Costo(hijo, Distancias, Num_var);
        end
    
        % Sustitución
        Nueva_Pob = zeros(Num_pob, Num_var);
        Nueva_Aptitud = zeros(Num_pob, 1);
        for i = 1:2:Num_pob
            familia = [Pob(i:i+1, :); Hijos((i+1)/2, :)];
            aptitudes = [Aptitud_Pob(i:i+1); Aptitud_Hijos((i+1)/2)];
            [~, indices] = sort(aptitudes);
            Nueva_Pob(i:i+1, :) = familia(indices(1:2), :);
            Nueva_Aptitud(i:i+1) = aptitudes(indices(1:2));
        end
    
        % Mezclar la nueva población
        indices = randperm(Num_pob);
        Pob = Nueva_Pob(indices, :);
        Aptitud_Pob = Nueva_Aptitud(indices);
    
        % Mutación
        for i = 1:Num_pob
            if rand <= Pm
                Pob(i, :) = randperm(Num_var);
                Aptitud_Pob(i) = Costo(Pob(i, :), Distancias, Num_var);
            end
        end
        
        % Actualizar mejor histórico si es necesario
        [Mejor_Aptitud_Generacion, idx] = min(Aptitud_Pob);
        if Mejor_Aptitud_Generacion < Mejor_Aptitud_Historico
            Mejor_Aptitud_Historico = Mejor_Aptitud_Generacion;
            Mejor_Individuo_Historico = Pob(idx, :);
        end
    
        % Información de salida de generación actual
        %fprintf('Iteración: %d, Mejor aptitud de la generación: %f\n', iter, Mejor_Aptitud_Generacion);
        iter = iter + 1;
    end
    
    % Mostrar resultados finales
    [Mejor_Aptitud_Final, idx] = min(Aptitud_Pob);
    %fprintf('Mejor solución final de la última generación: %s, Costo: %d\n', mat2str(Pob(idx, :)), Mejor_Aptitud_Final);
    fprintf('Solucion %d : %s, Costo: %d\n', iteracion, mat2str(Mejor_Individuo_Historico), Mejor_Aptitud_Historico);
    
    Resultados_Generales(iteracion) = Mejor_Aptitud_Historico;

    iteracion=iteracion+1;
end

%Se obtienen las estadisticas de los mejores reultados de cada iteracion
%disp(Mejor_Aptitud);
mejor = min(Resultados_Generales);
media = mean(Resultados_Generales);
peor = max(Resultados_Generales);
desviacion_estandar = std(Resultados_Generales);

disp('Los resultados generales del algoritmo genetico son: ');
disp(['Mejor: ', num2str(mejor)]);
disp(['Media: ', num2str(media)]);
disp(['Peor: ', num2str(peor)]);
disp(['Desviacion estandar: ', num2str(desviacion_estandar)]);

function Ciudades_Vecinas = CrearLista(padre1,  padre2, noCiudades)
%arreglo de celdas
Ciudades_Vecinas = cell(1,noCiudades); %arreglo de celdas
for i = 1:noCiudades
    ciudad = padre1(i);
    indx1 = [i-1, i+1];
    indx2 = find(padre2 == ciudad);

    %casos en los que indx estan fuera del limite
    indx1(indx1 == 0) = noCiudades;
    indx1(indx1 > noCiudades) = 1;
    indx2(indx2 == 0) = noCiudades;
    indx2(indx2 > noCiudades) = 1;

    %ciudades no repetidas
    vecinos = unique([padre1(indx1), padre2(indx2)]);
    Ciudades_Vecinas{ciudad} = vecinos;
end
end

function hijo = EdgeRecombination(padre1, padre2, noCiudades)
hijo = zeros(1, noCiudades);
Ciudades_Vecinas = CrearLista(padre1, padre2, noCiudades);
rd = round(rand);
city_actual = rd*padre1(1) + (1-rd)*padre2(1);
hijo(1) = city_actual;

for i = 2:noCiudades
    for j = 1:length(Ciudades_Vecinas)
        indx = Ciudades_Vecinas{j} == city_actual;
        Ciudades_Vecinas{j}(indx) = [];
    end

    vecinos_actuales = Ciudades_Vecinas{city_actual};
    if isempty(vecinos_actuales)
        cities = setdiff(1:noCiudades, hijo(1:i-1));
        city_actual = randsample(cities, 1);
    else
        noConex = [];
        for j = 1:length(vecinos_actuales)
            noConex(end+1) = length(Ciudades_Vecinas{vecinos_actuales(j)});
        end

        [~, Minindx] = find(noConex == min(noConex));
        idx = randsample(Minindx, 1);  % Elegir uno de los mínimos al azar
        city_actual = vecinos_actuales(idx);
    end
    hijo(i) = city_actual;
end
end

function costo = Costo(recorrido, Distancias, noCiudades)
    costo = 0;
    for i = 1:noCiudades-1
        costo = costo + Distancias(recorrido(i), recorrido(i+1));
    end
    costo = costo + Distancias(recorrido(noCiudades), recorrido(1)); % Retorno a la ciudad inicial
end


