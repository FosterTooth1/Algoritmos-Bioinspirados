% Parámetros
T = 100;

% Cargar la imagen
img = imread('elmalilla.jpg');

% Convertir a escala de grises
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Obtener las dimensiones de la imagen
[rows, cols] = size(img);

% Convertir a matriz de bits (los bits de la primera fila representan los pixeles de la primera fila...)
img_bits = zeros(rows, 8 * cols);
for i = 1:rows
    img_bits(i, :) = reshape(dec2bin(img(i, :), 8)', 1, cols * 8) - '0';
end



% Cifrar la imagen
img_cifrada = cifrar_imagen(img_bits, T);

% Convertir la imagen cifrada a valores de píxeles
img_cifrada_pix = zeros(rows, cols, 'uint8');
for i = 1:rows
    for j = 1:cols
        start_idx = (j - 1) * 8 + 1;
        end_idx = j * 8;
        binary_segment = img_cifrada(i, start_idx:end_idx);
        decimal_value = sum(binary_segment .* (2.^(7:-1:0)));
        img_cifrada_pix(i, j) = uint8(decimal_value);
    end
end

figure(1)
imshow(img_cifrada_pix);
title('Imagen cifrada- Autómata celular');


% Descifrar la imagen con animación
img_descifrada = descifrar_imagen(img_cifrada, T);
img_descifrada_pix = img_cifrada_pix; % Inicializar la imagen descifrada en píxeles

[rows1, cols1] = size(img_descifrada);
for i = 1:2:rows1-1
    % Convertir la fila descifrada actual a valores de píxeles para mostrar la animación
    for k = 1:cols
        start_idx = (k - 1) * 8 + 1;
        end_idx = k * 8;
        binary_segment = img_descifrada(i:i+1, start_idx:end_idx);
        decimal_value = sum(binary_segment .* (2.^(7:-1:0)), 2);
        img_descifrada_pix(i:i+1, k) = uint8(decimal_value);
    end

    % Mostrar la imagen descifrada parcialmente
    imshow(img_descifrada_pix);
    title('Descifrando Imagen - Autómata celular');
    drawnow;
    pause(10e-8); % Ajusta el valor para controlar la velocidad de la animación
end

% Función para cifrar la imagen
function img_bits = cifrar_imagen(img_bits, T)
    [rows, cols] = size(img_bits);
    for i = 1:2:rows-1
        for j = 1:T
            for k = 2:cols-1
                % Función para cifrar la imagen (regla 90R)
                temp(1,k-1) = mod(img_bits(i+1, k-1) + img_bits(i+1, k+1) + img_bits(i, k), 2);
            end
            img_bits(i, :) = img_bits(i+1, :);
            img_bits(i+1, 2:end-1) = temp;
        end
    end
end

% Función para descifrar la imagen
function img_descifrada = descifrar_imagen(img_cifrada, T)
    [rows, cols] = size(img_cifrada);
    img_descifrada = img_cifrada; % Inicializar la imagen descifrada

    for i = 1:2:rows-1
        for j = 1:T
            temp = zeros(1, cols - 2);
            for k = 2:cols-1
                % Función para descifrar la imagen (regla 90R)
                temp(k-1) = mod(img_descifrada(i, k-1) + img_descifrada(i, k+1) + img_descifrada(i+1, k), 2);
            end
            img_descifrada(i+1, :) = img_descifrada(i, :);
            img_descifrada(i, 2:end-1) = temp;
        end
    end
end