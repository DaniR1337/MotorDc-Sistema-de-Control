clc; clear variables, clear figures; 
% Integrantes: Daniel Alejandro Rodríguez Alvarado
% Carné: C06575

%% Extracción de las mediciones
M = readmatrix("delta_85a115.csv");
[t, u, y] = deal(M(:, 1), M(:, 2)*100/255, M(:, 3));

%% Gráficas iniciales
% Figura con respuesta del sistema
figure(1)
plot(t, y, 'm', 'linewidth', 2)
grid on

% Figura con señal del controlador
figure(2)
plot(t, u, 'r', 'linewidth', 2)
grid on

%% Alfaro (123c)
% Se identifica el valor inicial y final de la entrada y respuesta
yi = mean(y(1:51)); % Media de 'y' enre 0s y 1s
yf = mean(y(386:546)); % Media de 'y' entre 4.5s y 5.5s
ui = mean(u(1:51));
uf = mean(u(386:546));

% Se halla la ganancia K
K = (yf-yi)/(uf-ui);

% Se halla p1 (25%)
[~, ix25] = min(abs((y(1:500)-yi) - 0.25*(yf-yi))); 
[y25, t25] = deal(y(ix25), t(ix25));

% Se halla p2 (75%)
[~, ix75] = min(abs((y(1:500)-yi) - 0.75*(yf-yi))); 
[y75, t75] = deal(y(ix75), t(ix75));

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
tz = (0:0.01:10);


figure(3)
entrada = ui + heaviside(tz-3)*(uf-ui).*heaviside(6-tz) + heaviside(tz-9)*(uf-ui);
[yz, tz] = lsim(P, entrada, tz);
plot(tz, yz, 'linewidth', 2)
hold on
plot(t, y, 'linewidth', 2)
legend('Alfaro 123c', 'Real')
grid on


