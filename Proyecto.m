 % Integrantes: 
%   Daniel Alejandro Rodríguez Alvarado C06575
%   Nataly Delgado Huertas C02583
%   Sylvia Fonseca Cruz C03039

clc; clear variables, clear figures;

%% Extracción de las mediciones
M = readmatrix("delta_85a115.csv");
[t, u, yr] = deal(M(:, 1), M(:, 2), M(:, 3));

%% Gráficas iniciales
% Figura con respuesta del sistema
figure(1)
plot(t, yr, 'm', 'linewidth', 2)
grid on

% Figura con señal del controlador
figure(2)
plot(t, u, 'r', 'linewidth', 2)
grid on

%% Modelo POMTM Alfaro (123c)
% Se identifica el valor inicial y final de la entrada y respuesta
yi = mean(yr(1:51)); % Media de 'yr' enre 0s y 1s
yf = mean(yr(386:546)); % Media de 'yr' entre 4.5s y 5.5s
ui = mean(u(1:51));
uf = mean(u(386:546));

% Se halla la ganancia K
K = (yf-yi)/(uf-ui);

% Se halla p1 (25%)
[~, ix25] = min(abs((yr(1:500)-yi) - 0.25*(yf-yi))); 
[y25, t25] = deal(yr(ix25), t(ix25));

% Se halla p2 (75%)
[~, ix75] = min(abs((yr(1:500)-yi) - 0.75*(yf-yi))); 
[y75, t75] = deal(yr(ix75), t(ix75));

% Se hallan otros tiempos de interés
t0 = 3; 
t1 = t25-t0;
t2 = t75-t0;

% Modelo 
a = 0.910;
b = 1.262;
tau = a*(t2-t1);
L = b*t1 + (1-b)*t2;
s = tf('s');
P = K/(tau*s+1);
tz = (0:0.001:10);

% Se interpola la respuesta real
y = interp1(t, yr, tz);

% Se grafica la respuesta real con el modelo POMTM
figure(3)
entrada = heaviside(tz-3)*(uf-ui).*heaviside(6-tz) + heaviside(tz-9)*(uf-ui);
plot(tz, y, 'b', 'linewidth', 2)
hold on
[yz, tz] = lsim(P, entrada, tz);
plot(tz, yz+yi, 'r--', 'linewidth', 2)
plot(t, u, 'm', 'linewidth', 2)
legend('Real', 'POMTM Alfaro (123c)', 'Entrada')
grid on
% Índice integral de error absoluto
e = abs(y.'-(yz+yi));
JIAE1 = trapz(tz(2001:6001), e(2001:6001)); % Índice de error de 3s a 6s

%% Modelo PDMTM Alfaro (123c)
% Modelo 
a2 = 0.5776;
b2 = 1.5552;
tau2 = a2*(t2-t1);
L2 = b2*t1 + (1-b2)*t2;
P2 = K/(tau2*s+1)^2;
tz = (0:0.001:10);

% Se grafica la respuesta real con el modelo PDMTM
figure(4)
plot(tz, y, 'b', 'linewidth', 2)
hold on
[yz2, tz] = lsim(P2, entrada, tz);
plot(tz, yz2+yi, 'r--', 'linewidth', 2)
plot(t, u, 'm', 'linewidth', 2)
legend('Real', 'PDMTM Alfaro (123c))', 'Entrada')
grid on

% Índice integral de error absoluto
e2 = abs(y.'-(yz2+yi));
JIAE2 = trapz(tz(2001:6001), e2(2001:6001)); % Índice de error de 3s a 6s

