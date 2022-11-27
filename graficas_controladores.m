clc; clear variables, clear figures; 
% Integrantes: Daniel Alejandro Rodríguez Alvarado, Nataly Delgado Huertas y Sylvia Fonseca Cruz
% Carné: C06575, C02583 y C03039

% Método de sintonización LGR

M_LGR = readmatrix("Datos_LGR_Grupo02_07.txt");
[u_LGR, m_LGR, y_LGR] = deal(M_LGR(:, 3), M_LGR(:, 5), M_LGR(:, 7));
deltat_LGR = (6)/length(u_LGR);
t_LGR = transpose(0:deltat_LGR:6-deltat_LGR);


% Método de sintonización Síntesis Analítica (SA)

% Extracción de los datos controlador SA
M_SA = readmatrix("Datos_SA_Grupo02_07.txt");
[u_SA, m_SA, y_SA] = deal(M_SA(:, 3), M_SA(:, 5), M_SA(:, 7));
deltat_SA = 6/length(u_SA);
t_SA = transpose(0:deltat_SA:6-deltat_SA);


% Método de sintonización Klein

% Extracción de los datos controlador Klein
M_K = readmatrix("Datos_Klein_Grupo02_07.txt");
[u_K, m_K, y_K] = deal(M_K(:, 3), M_K(:, 5), M_K(:, 7));
deltat_K = 6/length(u_K);
t_K = transpose(0:deltat_K:6-deltat_K);


% Gráfica respuesta del controlador LGR
figure(1)
plot(t_LGR, y_LGR, t_LGR, m_LGR, t_LGR, u_LGR, 'LineWidth',1.5)
legend('Señal de salida', 'Señal de control', 'Señal de entrada')
grid on
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
title("Respuesta del controlador a partir del método de sintonización de LGR")


% Gráfica respuesta del controlador SA
figure(2)
plot(t_SA, y_SA, t_SA, m_SA, t_SA, u_SA, 'LineWidth',1.5)
legend('Señal de salida', 'Señal de control', 'Señal de entrada')
grid on
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
title("Respuesta del controlador a partir del método de sintonización de SA")


% Gráfica respuesta del controlador Klein
figure(3)
plot(t_K(1:end-2), y_K(1:end-2), t_K, m_K, t_K, u_K, 'LineWidth',1.5)
legend('Señal de salida', 'Señal de control', 'Señal de entrada')
grid on
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
title("Respuesta del controlador a partir del método de sintonización de Klein")


% Gráfica de la señal de salida de los 3 controladores
figure(4)
hold on
plot(t_LGR, u_LGR, 'LineWidth',1.5)
plot(t_SA, u_SA, 'LineWidth',1.5)
plot(t_K, u_K, 'LineWidth',1.5)
plot(t_LGR, y_LGR, 'LineWidth',1.5)
plot(t_SA, y_SA, 'LineWidth',1.5)
plot(t_K(1:end-2), y_K(1:end-2), 'LineWidth',1.5)
legend('r(t) LGR','r(t) SA','r(t) Klein', 'y(t) LGR', 'y(t) SA', 'y(t) Klein')
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
title("Señal de salida de de los controladores ante un valor de referencia")
hold off

% Gráfica de la señal de control de los 3 controladores
figure(5)
hold on
plot(t_LGR, u_LGR, 'LineWidth',1.5)
plot(t_SA, u_SA, 'LineWidth',1.5)
plot(t_K, u_K, 'LineWidth',1.5)
plot(t_LGR, m_LGR, 'LineWidth',1.5)
plot(t_SA, m_SA, 'LineWidth',1.5)
plot(t_K(7:end-2), m_K(7:end-2), 'LineWidth',1.5)
legend('r(t) LGR', 'r(t) SA', 'r(t) Klein','Señal de control LGR','Señal de control SA','Señal de control Klein')
xlabel("Tiempo (s)")
ylabel("Revoluciones por minuto (RPM)")
title("Señal de control de de los controladores ante un valor de referencia")
hold off
