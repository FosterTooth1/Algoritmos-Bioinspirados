clc
% Parámetros de entrada 
%Inicializar variables
Np=5; %Número de población
Nvar=2; %Número de variables
li=[-2 -2]; %limite izquierdo de x=-2 & y=-2
ls=[2 2]; %limite izquierdo de x=2 & y=2
Pr=2; %Cantidad de decimales

%Paso 1. Calcular el número de bits para representar los valores de -2 a 2
%con presición 2
for i=1:Nvar 
    Nbits(i)=fix(log2(ls(i)-li(i))*10^Pr+0.9); 
end

%Paso 2. Generar la población
Poblacion = randi([0 1],[Np sum(Nbits)]); %Esta línea genera una matriz de números
%aleatorios entre 0 y 1 con diensiones Np por la suma de los elementos en
%Nbits

%Paso 3. Decodificar los individuos para evaluar en FO

%aleloX = Poblacion(:,1:Nbits(1)); %substraemos los codigos que pertenecen al alelo x
%aleloY = Poblacion(:,Nbits(1)+1:end);%substraemos los codigos de pertenecen al alelo y

%Sacar la parte entera de X & Y
Xentero = round(bi2de(Poblacion(:,1:Nbits(1))));
Yentero = round(bi2de(Poblacion(:,Nbits(1)+1:end)));
%disp(Xentero);
%disp(Yentero);

%Aplicamos la segunda formula
for i=1:Np

    Xreal(i) = li(1) + ((Xentero(i) *  (ls(1) - li(1)))/(2^Nbits(1)-1));
    Yreal(i) = li(2) + ((Yentero(i) *  (ls(2) - li(2)))/(2^Nbits(2)-1));
end

disp(Xreal);
disp(Yreal);

%Evaluamos en la FO f(x,y) = (1-x)^2 + (100 - y)^2

fXY = (1-Xreal).^2 + (100-Yreal).^2 ;

disp(fXY);