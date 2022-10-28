clc; clear variables, clear figures; 

%% Extracción de las mediciones
M = readmatrix("delta_85a115.csv");
[t, u, yr] = deal(M(:, 1), M(:, 2), M(:, 3));

% Se identifica el valor inicial y final de la entrada y respuesta
yi = mean(yr(1:51)); % Media de 'yr' enre 0s y 1s
yf = mean(yr(386:546)); % Media de 'yr' entre 4.5s y 5.5s
ui = mean(u(1:51));
uf = mean(u(386:546));

s = tf('s'); 

P = 1.1598/(0.58415*s+1);
tz = (0:0.001:6);

% Se interpola la respuesta real
y = interp1(t, yr, tz);

% Se grafica la respuesta real con el modelo 123c de Alfaro
figure(3)
entrada = heaviside(tz-3)*(uf-ui); 
[yz, tz] = lsim(P, entrada, tz);

plot(tz, y, 'b', tz, yz+yi, 'r--', 'linewidth', 2)

% Índice integral de error absoluto
e = abs(y.'-(yz+yi));
JIAE = trapz(tz(3001:6001), e(3001:6001)); % Índice de error de 3s a 6s
disp(JIAE)

%% Prueba Alfaro
%Kp=0.86058, K=1.162, T=0.1847
s=tf('s')
Lalfaro = (0.86058*1.162)/(0.1847*s)
%sisotool(Lalfaro)
%% Sintesis analitica servo
tau=6.38

Kp=1/(tau*0.1847);
ti=0.1847;
td=0;
K=1.162;
T=0.1847;

Ls=(K*Kp*(ti*td*s^2+ti*s+1))/(ti*s*(T*s+1))
sisotool(Ls)

%% Método de Klein, controlador PI para POMTM
Tm = 0.1847;
Km = 1.162;
tau_m = 0.038;
Kc = 0.28*Tm/(Km*(tau_m+0.1*Tm));
Tc = 0.53*Tm;
Cc = Kc*(1+1/(Tc*s));
P=1.162/(0.1847*s+1);
%sisotool(Cc)