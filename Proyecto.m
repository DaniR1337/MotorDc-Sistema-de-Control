clc; clear variables, clear figures; 
% Integrantes: Daniel Alejandro Rodríguez Alvarado
% Carné: C06575

%% Extracción de las mediciones
M = readmatrix("delta_85a115.csv");
[t, u, yr] = deal(M(:, 1), M(:, 2)*100/255, M(:, 3));

%% Gráficas iniciales
% Figura con respuesta del sistema
figure(1)
plot(t, yr, 'm', 'linewidth', 2)
grid on

% Figura con señal del controlador
figure(2)
plot(t, u, 'r', 'linewidth', 2)
grid on

%% Alfaro (123c)
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
legend('Real', 'Alfaro 123c')
grid on

% Índice integral de error absoluto
e = abs(y.'-(yz+yi));
JIAE = trapz(tz(3001:6001), e(3001:6001)); % Índice de error de 3s a 6s
