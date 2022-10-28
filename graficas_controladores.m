%% Extracción de las mediciones datos controlador LGR
M_LGR = readmatrix("Datos_LGR_Grupo02_07.txt");
[u_LGR, m_LGR, y_LGR] = deal(M_LGR(:, 3), M_LGR(:, 5), M_LGR(:, 7));
deltat_LGR = 6/length(u_LGR);
t_LGR = transpose(0:deltat_LGR:6-deltat_LGR);
figure(1)
plot(t_LGR, y_LGR, t_LGR, u_LGR, t_LGR, m_LGR)

%% Extracción de las mediciones datos controlador SA
M_SA = readmatrix("Datos_SA_Grupo02_07.txt");
[u_SA, m_SA, y_SA] = deal(M_SA(:, 3), M_SA(:, 5), M_SA(:, 7));
deltat_SA = 6/length(u_SA);
t_SA = transpose(0:deltat_SA:6-deltat_SA);
figure(2)
plot(t_SA, y_SA, t_SA, m_SA, t_SA, u_SA)

%% Extracción de las mediciones datos controlador Klein
M_K = readmatrix("Datos_Klein_Grupo02_07.txt");
[u_K, m_K, y_K] = deal(M_K(:, 3), M_K(:, 5), M_K(:, 7));
deltat_K = 6/length(u_K);
t_K = transpose(0:deltat_K:6-deltat_K);
figure(3)
plot(t_K, y_K, t_K, m_K, t_K, u_K)

