clc; clear variables, clear figures; 
% Integrantes: Daniel Alejandro Rodríguez Alvarado, Nataly Delgado Huertas y Sylvia Fonseca Cruz
% Carné: C06575, C02583 y C03039

% Obtención del IAE para el modelo obtenido con Toolbox

% Extracción de las mediciones
M = readmatrix("delta_85a115.csv");
[t, u, yr] = deal(M(:, 1), M(:, 2), M(:, 3));

% Se identifica el valor inicial y final de la entrada y respuesta
yi = mean(yr(1:51)); % Media de 'yr' enre 0s y 1s
yf = mean(yr(386:546)); % Media de 'yr' entre 4.5s y 5.5s
ui = mean(u(1:51));
uf = mean(u(386:546));

s = tf('s'); 

% Función de transferencia
P = 1.1453/(0.54856*s+1);
tz = (0:0.001:6);

% Se interpola la respuesta real
y = interp1(t, yr, tz);

% Se grafica la respuesta real con el modelo del Toolbox
figure(1)
entrada = heaviside(tz-3)*(uf-ui); 
[yz, tz] = lsim(P, entrada, tz);

plot(tz, y, 'b', tz, yz+yi, 'r--', 'linewidth', 2)

% Índice integral de error absoluto IAE
e = abs(y.'-(yz+yi));
JIAE = trapz(tz(3001:6001), e(3001:6001)); % Índice de error de 3s a 6s
disp(JIAE)