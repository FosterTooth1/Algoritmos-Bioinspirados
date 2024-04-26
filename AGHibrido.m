function []= AGHibrido()
clear, clc, close all;

Num_pob=100;
Num_Gen=500;
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

Soluciones=[];

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

iter= 1;
while iter <= Num_Gen
    %Seleccion
    idx = randi(Num_pob, [1 Num_pob])';
    Padres = Pob(idx, :);
    Aptitud_Padres = Aptitud_Pob(idx);

    %Cruzamiento
    Hijos = [];
    for i = 1:2:Num_pob
        Hijo = EdgeRecombination(Padres(i,:), Padres(i+1,:), Num_var);
        Hijos = [Hijos; Hijo];
    end

    %Evaluacion
    Aptitud_Hijos = zeros(Num_pob,1);
    for i = 1:Num_pob/2
        Aptitud_Hijos(i) = Costo(Hijos(i,:), Distancias, Num_var);
    end

    %Sustitucion
    Pob = [];
    Aptitud_Pob = [];
    j = 1;
    for i = 1:2:Num_pob
        %Familia
        Familia = [Padres(i:i+1, :); Hijos(j,:)];
        %Ordenamiento de la familia de mejor a peor
        [v, idx] = sort([Aptitud_Padres(i:i+1); Aptitud_Hijos(j)]);
        %Ganadores de la familia
        Pob(i:i+1,:) = Familia(idx(1:2),:);
        %Aptitudes
        Aptitud_Pob(i:i+1,1) = v(1:2);
        j = j + 1;
    end  

    %Mutacion
    for i = 1:Num_pob
        if rand <= Pm
            ind = randi(Num_pob);
            Pob(ind,:) = randperm(Num_var);
            Aptitud_Pob(ind) = Costo(Pob(ind,:), Distancias, Num_var);
        end
    end

    %Elite
    [~, ind] = sort(Aptitud_Pob);
    p_Elite(iter,:) = Pob(ind(1),:);
    Aptitud_Elite(iter) = Aptitud_Pob(ind(1));
    
    % Salida
    fprintf('Iteración: %d, p_Elite: %s, Aptitud_Elite: %f\n', iter, mat2str(p_Elite(iter,:)), Aptitud_Elite(iter));
    iter = iter + 1;
end

Mejor_Historico = find(Aptitud_Elite == min(Aptitud_Elite));
Pos_Mejor = Mejor_Historico(1);

Solucion = [p_Elite(Num_Gen,:), Aptitud_Elite(Num_Gen)];

fprintf('Mejor ultima generacion de recorrido: %s, Costo: %d\n', mat2str(Solucion(1:end-1)), Solucion(end));
fprintf('Mejor recorrido Historico: %s, Costo: %d\n', mat2str(p_Elite(Pos_Mejor,:)), Aptitud_Elite(Pos_Mejor));

end

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


