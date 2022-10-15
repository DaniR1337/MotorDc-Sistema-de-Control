clc; clear figures; 

% Extracción de las mediciones
M = readmatrix("delta_85a115.csv");
[t, u, y] = deal(M(:, 1), M(:, 2), M(:, 3));

% Figura con respuesta del sistema
figure(1)
plot(t, y, 'm')
grid on

% Figura con señal del controlador
figure(2)
plot(t, u, 'r')
grid on