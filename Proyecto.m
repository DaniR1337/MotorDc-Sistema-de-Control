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

%% Modelo Alfaro (123c)
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

% Se grafica la respuesta real con el modelo 123c de Alfaro
figure(3)
entrada = heaviside(tz-3)*(uf-ui).*heaviside(6-tz) + heaviside(tz-9)*(uf-ui);
plot(tz, y, 'b', 'linewidth', 2)
hold on
[yz, tz] = lsim(P, entrada, tz);
plot(tz, yz+yi, 'r--', 'linewidth', 2)
plot(t, u, 'm', 'linewidth', 2)
legend('Real', 'Alfaro (123c)', 'Entrada')
grid on
% Índice integral de error absoluto
e = abs(y.'-(yz+yi));
JIAE1 = trapz(tz(3001:6001), e(3001:6001)); % Índice de error de 3s a 6s

%% Modelo Ho
% Se halla p1 (35%)
[~, ix35] = min(abs((yr(1:500)-yi) - 0.35*(yf-yi))); 
[y35, t35] = deal(yr(ix35), t(ix35));

% Se halla p2 (85%)
[~, ix85] = min(abs((yr(1:500)-yi) - 0.85*(yf-yi))); 
[y85, t85] = deal(yr(ix85), t(ix85));

% Se hallan otros tiempos de interés 
t3 = t35-t0;
t4 = t85-t0;

% Modelo 
aHo = 0.670;
bHo = 1.290;
tauHo = aHo*(t4-t3);
LHo = bHo*t3 + (1-b)*t4;
P2 = K/(tauHo*s+1);
tz = (0:0.001:10);

% Se grafica la respuesta real con el modelo Ho
figure(4)
plot(tz, y, 'b', 'linewidth', 2)
hold on
[yz2, tz] = lsim(P2, entrada, tz);
plot(tz, yz2+yi, 'r--', 'linewidth', 2)
plot(t, u, 'm', 'linewidth', 2)
legend('Real', 'Ho', 'Entrada')
grid on

% Índice integral de error absoluto
e2 = abs(y.'-(yz2+yi));
JIAE2 = trapz(tz(3001:6001), e2(3001:6001)); % Índice de error de 3s a 6s

