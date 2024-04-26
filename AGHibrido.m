function []= AGHibrido()
clear, clc, close all;

Np=50;
Ng=200;
Pm=0.2;
Nvar=11;
m=3;

d = [
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

runs=1;
Soluciones=[];

%Poblacion inicial
p=zeros(Np, Nvar);
for i=1:Np
    p(i,:)= randperm(Nvar);
end

%Evaluación
VFO_p=zeros(Np, 1);
for i=1:Np
    VFO_p(i)= FO_TSP(p(i,:), d, Nvar);
end

iter= 1;
while iter<= Ng
    %Seleccion
    idx=randi(Np, [1 Np])';
    Padres=p(idx, :);
    VFO_padres=VFO_p(idx);
    
    %Cruzamiento
    Hijos=[];
    for i=1:2:Np
        Hijo = EdgeRecombination(Padres(i,:), Padres(i+1,:), Nvar);
        Hijos=[Hijos; Hijo];
    end

    %Evaluacion
    VFO_h=zeros(Np,1);
    for i=1:Np/2
        VFO_h(i)=FO_TSP(Hijos(i,:),d,Nvar);
    end
    
    %Sustitucion
    p=[];
    VFO_p=[];
    j = 1;
    for i=1:2:Np
            %Familia
            Familia=[Padres(i:i+1, :); Hijos(j,:)];
            %Ordenamiemto de la familia de mejor a peor
            [v,idx]=sort([VFO_padres(i:i+1); VFO_h(j)]);
            %Ganadores de la familia
            p(i:i+1,:)=Familia(idx(1:2),:);
            %Aptitudes
            VFO_p(i:i+1,1)=v(1:2);
            j=j+1;
    end
    %Mutacion
    for i=1:Np
        if rand<=Pm
            ind=randi(Np);
            p(ind,:)=randperm(Nvar);
            VFO_p(ind)=FO_TSP(p(ind,:), d,Nvar);
        end
    end
    %Elite
    [~,ind]= sort(VFO_p);
    p_Elite=p(ind(1),:);
    VFO_Elite=VFO_p(ind(1));
    %Salida
    iter=iter+1;

end

Soluciones = [Soluciones; p_Elite, VFO_Elite];  % Añadir permutación y costo

% Antes de mostrar la salida
disp('Permutaciones y sus costos:');
for i = 1:size(Soluciones, 1)
    fprintf('Permutación: %s, Costo: %d\n', mat2str(Soluciones(i, 1:end-1)), Soluciones(i, end));
end

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

function costo = FO_TSP(recorrido, d, noCiudades)
    costo = 0;
    for k = 1:noCiudades-1
        costo = costo + d(recorrido(k), recorrido(k+1));
    end
    costo = costo + d(recorrido(noCiudades), recorrido(1)); % Retorno a la ciudad inicial
end


