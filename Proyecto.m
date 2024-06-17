clear, clc, close all;

% Solicitar al usuario el número de veces que se debe repetir el algoritmo
Num_iteraciones = input('Ingrese el número de veces que se repetirá el algoritmo: ');

% Creamos un arreglo donde se almacenaran los individuos con mejor aptitud (minimización) por cada iteración
Resultados_Generales= zeros(1,Num_iteraciones);

% Leer el archivo CSV
data = readmatrix('distancias.csv');

% Eliminar la primera columna, que contiene los encabezados laterales
Distancias = data(:, 2:end);

Permutacion_ciudades= zeros(Num_iteraciones, size(Distancias, 1));

for iteracion = 1:Num_iteraciones
    
    Num_pob=100;
    Num_Gen=200;
    Pm=0.15;
    Num_var=size(Distancias, 1);  % Número de ciudades basado en la matriz de distancias
    
    % Población inicial
    Pob=zeros(Num_pob, Num_var);
    for i=1:Num_pob
        Pob(i,:)= randperm(Num_var);
    end
    
    % Evaluación
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
            Nueva_Aptitud(i:i+1) = [aptitudes(indices(1:2))];
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
    Permutacion_ciudades(iteracion,:) = Mejor_Individuo_Historico;


    iteracion=iteracion+1;
end

% Se obtienen las estadísticas de los mejores resultados de cada iteración
mejor = min(Resultados_Generales);
media = mean(Resultados_Generales);
peor = max(Resultados_Generales);
desviacion_estandar = std(Resultados_Generales);

disp('Los resultados generales del algoritmo genético son: ');
disp(['Mejor: ', num2str(mejor)]);
disp(['Media: ', num2str(media)]);
disp(['Peor: ', num2str(peor)]);
disp(['Desviación estándar: ', num2str(desviacion_estandar)]);

% Leer el archivo CSV
fid = fopen('distancias.csv');
nombres_ciudades_cell = textscan(fid, '%s', 'Delimiter', ',');
fclose(fid);

% Convertir a matriz de caracteres (string)
nombres_ciudades = char(nombres_ciudades_cell{1}(2:end)); % Convertir y omitir la primera fila

% Obtener solo el número de ciudades necesarias
num_ciudades = size(Distancias, 1);  % Número de ciudades basado en la matriz de distancias
nombres_ciudades = nombres_ciudades(1:num_ciudades, :); % Desde la segunda celda hasta el número deseado

% Imprimir el orden de ciudades recorridas por la mejor solución encontrada
[~, idx] = min(Resultados_Generales);
mejor_solucion = Permutacion_ciudades(idx,:);

fprintf('Orden de ciudades recorridas:\n');
for i = 1:length(mejor_solucion)
    nombre_ciudad = nombres_ciudades(mejor_solucion(i), :);
    fprintf('%d. %s\n', i, nombre_ciudad);
end




function Ciudades_Vecinas = CrearLista(padre1,  padre2, noCiudades)
% arreglo de celdas
Ciudades_Vecinas = cell(1,noCiudades); %arreglo de celdas
for i = 1:noCiudades
    ciudad = padre1(i);
    indx1 = [i-1, i+1];
    indx2 = find(padre2 == ciudad);

    % casos en los que indx estan fuera del limite
    indx1(indx1 == 0) = noCiudades;
    indx1(indx1 > noCiudades) = 1;
    indx2(indx2 == 0) = noCiudades;
    indx2(indx2 > noCiudades) = 1;

    % ciudades no repetidas
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
