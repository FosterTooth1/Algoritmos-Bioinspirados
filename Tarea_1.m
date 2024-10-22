%%Flores Lara Alberto 5BV1
%Tarea 2

clc

% Solicitar al usuario el número de veces que se debe repetir el algoritmo
Num_iteraciones = input('Ingrese el número de veces que se repetirá el algoritmo: ');

% Solicitar al usuario el número de generaciones en el algoritmo
Num_Generaciones = input('Ingrese el número de generaciones en el algoritmo: ');

%%Establecemos nuestra poblacion y cantidad de variables con las que vamos a trabajar
Num_pob = input('Ingrese el numero de individuos dentro de la población: ');
Num_var = input('Ingrese la cantidad de variables que tiene cada individuo: ');

%%Establecemos los limites inferiores y superiores de las variables en un
%%ciclo for y los guardamos en un arreglo

Limite_Inferior = zeros(1, Num_var);
for i = 1:Num_var
    valor = input(['Ingrese el limite inferior de la variable ', num2str(i), ': ']);
    Limite_Inferior(i) = valor;
end


Limite_Superior = zeros(1, Num_var);
for i = 1:Num_var
    % Solicita al usuario que ingrese un valor
    valor = input(['Ingrese el limite superior de la variable ', num2str(i), ': ']);
    Limite_Superior(i) = valor;
end


%Precicion del algoritmo
Precision=input('Ingrese el numero de precision del algoritmo: ');

%Probabilidad de cruza
Pc=input('Ingrese la probabilidad de cruza del algoritmo: ');

%Probabilidad de mutación
Pm=input('Ingrese la probabilidad de mutacion del algoritmo: ');

Num_bits= zeros(1,Num_var);
    
%Calculamos la cantidad de bits para cada variable 
for i=1:Num_var 
    Num_bits(i)=fix(log2((Limite_Superior(i)-Limite_Inferior(i)) * 10 ^ Precision) + 0.9); 
end

%Creamos un arreglo donde se almacenaran los individuos con mejor aptitud (minimzacion) por cada iteracion
Resultados_Generales= zeros(1,Num_iteraciones);

for iteracion = 1:Num_iteraciones

    %Generamos la población en una matriz m*n donde m es la poblacion y n son
    %la suma de los bits de todas las variables
    Poblacion = randi([0 1],[Num_pob sum(Num_bits)]);
    
    %Creamos un arreglo donde se almacenaran los individuos con mejor aptitud (minimzacion) por cada generacion
    Mejor_Aptitud = zeros(1,Num_Generaciones);

    %Aqui se regulan la cantidad de generaciones del algoritmo
    for generacion = 1:Num_Generaciones
        %disp(['Iteración número: ', num2str(generacion)]);
        %disp(Poblacion)

        % Convertir la población binaria a población real
        Poblacion_real = zeros(Num_pob, Num_var);

        %Iteramos dentro de cada variable en todo individuo
        for i = 1:Num_pob
            bit_inicio = 1;
            for j = 1:Num_var
                %Buscamos el rango de los bits a los que pertenece cada variable
                bit_final = bit_inicio + Num_bits(j) - 1;
                binario = Poblacion(i, bit_inicio:bit_final);
                %Lo convertimos de binario a decimal
                %disp('Binario')
                %disp(binario)
                valor_binario = binario_a_decimal(binario);
                %disp('Decimal')
                %disp(valor_binario)
                %Hacemos la conversion a reales
                Poblacion_real(i, j) = Limite_Inferior(j) + ((valor_binario * (Limite_Superior(j) - Limite_Inferior(j) )) / ((2 ^ Num_bits(j))-1));
                %disp("Valor real")
                %disp(Poblacion_real(i,j))
                bit_inicio = bit_final + 1;
            end
        end
        
        aptitud= zeros(1,Num_pob);
        % Evaluación de la población en la función objetivo
        for i = 1:Num_pob
            %aptitud(i) = (1-Poblacion_real(i,1))^2 + (100-Poblacion_real(i,2))^2;
            aptitud(i) = 20 + (Poblacion_real(i,1)^2 - 10 * cos(2 * pi * Poblacion_real(i,1))) + (Poblacion_real(i,2)^2 - 10 * cos(2 * pi * Poblacion_real(i,2)));
        end
        
        % Gráfica de la población en el espacio 2D para la generación actual
        figure(iteracion); % Crear una nueva figura para cada iteración
        scatter(Poblacion(:, 1), Poblacion(:, 2), 50, 'filled'); % Graficar la población actual
        title(['Iteración ', num2str(iteracion), ' - Generación ', num2str(generacion)]);
        xlabel('Variable X');
        ylabel('Variable Y');
        xlim([Limite_Inferior(1), Limite_Superior(1)]);
        ylim([Limite_Inferior(2), Limite_Superior(2)]);
        grid on;
        drawnow; % Actualizar la gráfica en cada generación

        %{
        disp(['Esta es la generacion numero: ', num2str(generacion)]);
        for i = 1:Num_pob
            disp(['Esta es la aptitud del individuo ', num2str(i), ': ']);
            disp(aptitud(i));
        end
        %}

        %Obtener la posicion de la poblacion con mejor aptitud (minimzacion)
        posiciones_menor = find(aptitud == min(aptitud));
        Pos_Mejor = posiciones_menor(1);
        %disp(Pos_Mejor)
        %disp(['El individuo en la posicion ', num2str(Pos_Mejor), ' tiene la mejor aptitud con: ', num2str(aptitud(Pos_Mejor))]);
        %Almacenar la mejor aptitud por generacion
        Mejor_Aptitud(generacion)= aptitud(Pos_Mejor);

        %Codigo del Torneo
        Padres = zeros(Num_pob,sum(Num_bits));
        Torneo = [randperm(Num_pob); randperm(Num_pob)]';

        for i = 1:Num_pob
            if aptitud(Torneo(i,1))<aptitud(Torneo(i,2))
                Padres(i,:) = Poblacion(Torneo(i,1),:);

            else
                Padres(i,:) = Poblacion(Torneo(i,2),:);
            end
        end

        %{
        for i = 1:Num_pob
            disp(['Esta es la codificacion en bits del Padre ', num2str(i), ':']);
            disp(Padres(i,:));
        end
        %}

        %Cruzamiento
        Hijos=zeros(Num_pob, sum(Num_bits));
        for i =1:2:Num_pob-1
            %disp(['Padres :', num2str(i), num2str(i+1)]);
            Num_random= rand;

            if Num_random <= Pc
                puntos=randperm(sum(Num_bits)-1, 2);
                puntos=sort(puntos);
                %disp('Los puntos de corte son:')
                %disp(puntos);
                hijo_1= [Padres(i,1:puntos(1)), Padres(i+1,puntos(1)+1:puntos(2)), Padres(i,puntos(2)+1:end)];
                hijo_2= [Padres(i+1,1:puntos(1)), Padres(i,puntos(1)+1:puntos(2)), Padres(i+1,puntos(2)+1:end)];
            else
                %disp('Los dos hijos son los padres');
                hijo_1 = Padres(i,:);
                hijo_2 = Padres(i+1,:);
            end
            Hijos(i,:) = hijo_1;
            Hijos(i+1,:) = hijo_2;

        end

        if mod(Num_pob, 2) ~= 0
            Hijos(Num_pob,:)=Padres(Num_pob,:);
        end
        
        %{
        disp('Estos son los Hijos:')
        for i = 1:Num_pob
            disp(['Esta es la codificacion en bits del Hijo ', num2str(i), ':']);
            disp(Hijos(i,:));
        end
        %}

        %Mutacion
        for i=1:Num_pob
            Num_random= rand;
            if Num_random <= Pm
                x = randi([1, sum(Num_bits)]);
                Hijos(i, x) = ~Hijos(i, x);
            end
        end
        
        %{
        disp('Estos son los Hijos tras la mutacion:')
        for i = 1:Num_pob
            disp(['Esta es la codificacion en bits del Hijo ', num2str(i), ':']);
            disp(Hijos(i,:));
        end
        %}

        %Sustitucion
        Pos_aleatoria = randi([1, Num_pob]);
        Hijos(Pos_aleatoria,:) = Poblacion(Pos_Mejor,:);
        
        %Convertimos a los hijos en la poblacion de la siguiente generacion
        if generacion < Num_Generaciones
            Poblacion = Hijos;
        end
    end

    
    %Se obtienen las estadisticas de los resultados obtenidos en cada generacion
    %disp(Mejor_Aptitud);
    mejor = min(Mejor_Aptitud);
    %{
    media = mean(Mejor_Aptitud);
    peor = max(Mejor_Aptitud);
    desviacion_estandar = std(Mejor_Aptitud);

    disp(['Los resultados de la iteracion numero ',num2str(iteracion),' del algoritmo genetico son: ']);
    disp(['Mejor: ', num2str(mejor)]);
    disp(['Media: ', num2str(media)]);
    disp(['Peor: ', num2str(peor)]);
    disp(['Desviacion estandar: ', num2str(desviacion_estandar)]);
    disp('');
    %}
    Resultados_Generales(iteracion) = mejor;
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


function decimal = binario_a_decimal(binario)
    % Inicializa el valor decimal
    decimal = 0;
    % Itera sobre cada dígito binario en la cadena
    for i = 1:length(binario)
        % Suma el valor del dígito binario (0 o 1) multiplicado por su posición en la cadena
        decimal = decimal + binario(i)*2^(length(binario)-i);
    end
end