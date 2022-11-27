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
P3 = 1.1444/(0.16619*s+1);
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
JIAE3 = trapz(tz(3001:6001), e(3001:6001)); % Índice de error de 3s a 6s
disp(JIAE3)

%% Gráficas de los modelos juntas con base a los métodos de identificación
figure(6)

plot(t, u, 'm--', 'linewidth', 1.5)
hold on
entrada = heaviside(tz-3)*(uf-ui).*heaviside(6-tz) + heaviside(tz-9)*(uf-ui)+heaviside(tz-20)*(uf-ui);
[yz, tz] = lsim(P, entrada, tz);
plot(tz, yz+85, 'linewidth', 1.5)
[yz2, tz] = lsim(P2, entrada, tz);
plot(tz, yz2+85, 'linewidth', 1.5)
entrada = heaviside(tz-3)*(uf-ui)-heaviside(tz-6)*(uf-ui)+heaviside(tz-9)*(uf-ui); 
[yz3, tz] = lsim(P3, entrada, tz);
plot(tz, yz3+85, 'linewidth', 1.5)
grid on
title("Respuesta de los modelos ante la señal deseada")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("r(t): Valor de referencia", "Alfaro 123c POMTM", "Alfaro 123c PDMTM","Toolbox Matlab")

%% Gráficas de los controladores para los tres métodos de sintonización

% -Método de LGR, controlador PI para POMTM, Respuesta del controlador a lazo cerrado-
% Kp=0.86058, K=1.1619, T=0.1847
s=tf('s');
L_LGR = (0.86058*1.1619)/(0.1847*s);
Myr_LGR = L_LGR/(1+L_LGR);
% Gráfica
figure(7)
hold on
[yz4, ~] = lsim(Myr_LGR, entrada, tz);
plot(t, u, 'm--', 'linewidth', 1.5);
plot(tz, yz4+85, 'b', 'linewidth', 2);
grid on
title("Respuesta a lazo cerrado con base al método del LGR")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("r(t)", "y(t)")

% -Método de SA, controlador PID estándar con Td=0 para POMTM, Respuesta
% del controlador a lazo cerrado-
% Kp=1/6.38, K=1.1619, T=0.1847
s=tf('s');
K_SA= 1.1619;
T_SA=0.1847;
Kp_SA=1/(6.38*T_SA);
C_SA = Kp_SA*((T_SA*s+1)/(T_SA*s));
P_SA = K_SA/(T_SA*s+1);
L_SA = C_SA * P_SA;
Myr_SA = L_SA/(1+L_SA);
% Gráfica
figure(8)
hold on
[yz5, ~] = lsim(Myr_SA, entrada, tz);
plot(t, u, 'm--', 'linewidth', 1.5);
plot(tz, yz5+85, 'b', 'linewidth', 2);
grid on
title("Respuesta a lazo cerrado con base al método SA")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("r(t)", "y(t)")

% -Método de Klein, controlador PI para POMTM, Respuesta del controlador a
% lazo cerrado-
% Parámetros
s=tf('s');
K_K = 1.1619;
T_K = 0.1847;
tau_K = 0.038;
Kc_K = 0.28*T_K/(K_K*(tau_K+0.1*T_K));
Tc_K = 0.53*T_K;
C_K = Kc_K*((Tc_K*s+1)/(Tc_K*s));
P_K = K_K/(T_K*s+1);
L_K = C_K * P_K;
Myr_K = L_K/(1+L_K);
% Gráfica
figure(9)
hold on
[yz6, ~] = lsim(Myr_K, entrada, tz);
plot(t, u, 'm--', 'linewidth', 1.5);
plot(tz, yz6+85, 'b', 'linewidth', 2);
grid on
title("Respuesta a lazo cerrado con base al método Klein")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("r(t)", "y(t)")

% --Gráfica de los tres métodos de sintonización, Respuesta a lazo
% cerrado--
figure(10)
plot(t, u, 'k--', 'linewidth', 1.5)
hold on
entrada = heaviside(tz-3)*(uf-ui).*heaviside(6-tz) + heaviside(tz-9)*(uf-ui)+heaviside(tz-20)*(uf-ui);
[yz4, ~] = lsim(Myr_LGR, entrada, tz);
plot(tz, yz4+85, 'm','linewidth', 2)
[yz5, ~] = lsim(Myr_SA, entrada, tz);
plot(tz, yz5+85, 'b--','linewidth', 2)
[yz6, ~] = lsim(Myr_K, entrada, tz);
plot(tz, yz6+85, 'g','linewidth', 2)
grid on
title("Respuesta a lazo cerrado de los métodos de sintonización")
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
legend("r(t)","y(t) LGR", "y(t) SA", "y(t) Klein")


%% Métodos de Sintonización-Optimización para robustez y esfuerzo de control
% Datos importantes de cada controlador
s=tf('s');
P=1.1619/(0.1847*s+1);
C_LGR = 0.86058*(1+(1/(0.1847*s)));
C_SA = (1/(6.38*0.1847))*((0.1847*s+1)/(0.1847*s));
Kc_K = 0.28*0.1847/(1.1619*(0.038+0.1*0.1847));
C_K = Kc_K*(((0.53*0.1847)*s+1)/((0.53*0.1847)*s));
H_LGR = 1 + P*C_LGR;
H_SA = 1 + P*C_SA;
H_K = 1 + P*C_K;
%margin(H);
S_LGR = 1/(H_LGR);
S_SA = 1/(H_SA);
S_K = 1/(H_K);
%bode(G)

%-------------------
% Robustez del controlador diseñado, sensibilidad máxima
gpeak_LGR = getPeakGain(S_LGR)
gpeak_SA = getPeakGain(S_SA)
gpeak_K = getPeakGain(S_K)

%-------------------
% Esfuerzo de control
Mur_LGR = C_LGR / (1+C_LGR*P);
Mur_SA= C_SA / (1+C_SA*P);
Mur_K = C_K / (1+C_K*P);

% Extracción de los datos obtenidos en las pruebas reales del controlador
M_LGR = readmatrix("Datos_LGR_Grupo02_07.txt");
[u_LGR, m_LGR, y_LGR] = deal(M_LGR(:, 3), M_LGR(:, 5), M_LGR(:, 7));
deltat = 3.3/389;
t = transpose(0:deltat:3.3-deltat);
M_SA = readmatrix("Datos_SA_Grupo02_07.txt");
[u_SA, m_SA, y_SA] = deal(M_SA(:, 3), M_SA(:, 5), M_SA(:, 7));
M_K = readmatrix("Datos_Klein_Grupo02_07.txt");
[u_K, m_K, y_K] = deal(M_K(:, 3), M_K(:, 5), M_K(:, 7));

r = 30*heaviside(t-1); 
ur_LGR = lsim(Mur_LGR,r,t);
ur_SA = lsim(Mur_SA,r,t);
ur_K = lsim(Mur_K,r,t);

% Valores del esfuerzo de control
TVur = sum(abs(diff(ur_LGR)))
TVur = sum(abs(diff(ur_SA)))
TVur = sum(abs(diff(ur_K)))


%-------------------
% Información del desempeño sistema de cada controlador 
stepinfo(y_LGR(295:683),t,115)
stepinfo(y_SA(295:683),t,115)
stepinfo(y_K(295:683),t,115)

%-------------------
% Error permamente de cada controlador
error_LGR=(abs(mean(y_LGR(583:683)-115))/115)*100
error_SA=(abs(mean(y_SA(583:683)-115))/115)*100
error_K=(abs(mean(y_K(583:683)-115))/115)*100

