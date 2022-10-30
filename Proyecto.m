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
disp(JIAE1)

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
disp(JIAE2)

%% Modelo System Toolbox 

% Función de transferencia
P3 = 1.1453/(0.54856*s+1);
tz = (0:0.001:10);

% Se interpola la respuesta real
y = interp1(t, yr, tz);

% Se grafica la respuesta real con el modelo del Toolbox
figure(5)
entrada = heaviside(tz-3)*(uf-ui)-heaviside(tz-6)*(uf-ui)+heaviside(tz-9)*(uf-ui); 
[yz3, tz] = lsim(P3, entrada, tz);

plot(tz, y, 'b', tz, yz3+yi, 'r--', 'linewidth', 2)

% Índice integral de error absoluto IAE
e = abs(y.'-(yz3+yi));
JIAE3 = trapz(tz3(3001:6001), e(3001:6001)); % Índice de error de 3s a 6s
disp(JIAE3)

%% Gráficas de los modelos juntas con base a los métodos de identificación

figure(6)

plot(t, u, 'm--', 'linewidth', 1.5)
hold on
entrada = heaviside(tz-3)*(uf-ui).*heaviside(6-tz) + heaviside(tz-9)*(uf-ui);
[yz, tz] = lsim(P, entrada, tz);
plot(tz, yz+yi, 'linewidth', 1.5)
[yz2, tz] = lsim(P2, entrada, tz);
plot(tz, yz2+yi, 'linewidth', 1.5)
entrada = heaviside(tz-3)*(uf-ui)-heaviside(tz-6)*(uf-ui)+heaviside(tz-9)*(uf-ui); 
[yz3, tz] = lsim(P3, entrada, tz);
plot(tz, yz3+yi, 'linewidth', 1.5)
grid on
title("Resuesta de los modelos ante la señal deseada")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("Valor de referencia", "Alfaro 123c POMTM", "Alfaro 123c PDMTM","Toolbox Matlab")

%%  Método de LGR, controlador PI para POMTM, Respuesta del controlador a lazo cerrado
% Kp=0.86058, K=1.1619, T=0.1847
s=tf('s');
L_LGR = (0.86058*1.1619)/(0.1847*s);
Myr_LGR = L_LGR/(1+L_LGR);

% Gráfica
figure(7)
hold on
plot(t, u, 'm--', 'linewidth', 1.5)
[yz4, ~] = lsim(Myr_LGR, entrada, tz);
plot(tz, yz4+yi, 'b', 'linewidth', 2)
grid on
title("Respuesta del controlador con base al método del LGR")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("Valor de referencia", "Señal del LGR")
