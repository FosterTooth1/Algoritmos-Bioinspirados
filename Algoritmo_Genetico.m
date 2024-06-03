clc;
clear;
warning off all;

Ng = 50; % Número de generaciones
Ncells = 101; % Número de células (debe ser impar para tener un centro)

% Estado inicial
cells = zeros(Ng, Ncells);
cells(1, ceil(Ncells/2)) = 1; % Iniciar con una célula viva en el centro

% Definir la regla 150
regla = [0 1 1 0 1 0 0 1];

% Función para obtener el nuevo estado de la célula
nuevo_estado = @(vec) regla(bin2dec(num2str(vec)) + 1);

% Simulación del autómata celular
for gen = 2:Ng
    % Actualizar células con una convolución
    for i = 2:Ncells-1
        cells(gen, i) = nuevo_estado(cells(gen-1, i-1:i+1));
    end
    cells(gen, 1) = nuevo_estado([0, cells(gen-1, 1), cells(gen-1, 2)]);
    cells(gen, Ncells) = nuevo_estado([cells(gen-1, Ncells-1), cells(gen-1, Ncells), 0]);
end

% Visualizar el autómata celular
figure;
imagesc(cells);
colormap(abyss);
title('Autómata celular - Regla Mod 2');
xlabel('Células');
ylabel('Generaciones');
axis equal tight;