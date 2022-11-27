clc; clear variables, clear figures; 
% Integrantes: Daniel Alejandro Rodríguez Alvarado, Nataly Delgado Huertas y Sylvia Fonseca Cruz
% Carné: C06575, C02583 y C03039

% Pruebas de los distintos modelos de sintonización 

%%  Método de LGR, controlador PI para POMTM
% Kp=0.86058, K=1.162, T=0.1847
s=tf('s');
L_LGR = (0.86058*1.162)/(0.1847*s);
sisotool(L_LGR)

%%  Método de Síntesis Analitica para servo
s=tf('s');
tau=6.38;
Kp=1/(tau*0.1847);
ti=0.1847;
td=0;
K=1.162;
T=0.1847;

Ls=(K*Kp*(ti*td*s^2+ti*s+1))/(ti*s*(T*s+1));
sisotool(Ls)

%% Método de Klein, controlador PI para POMTM
s=tf('s');
Tm = 0.1847;
Km = 1.162;
tau_m = 0.038;
Kc = 0.28*Tm/(Km*(tau_m+0.1*Tm));
Tc = 0.53*Tm;
Cc = Kc*(1+1/(Tc*s));
P=1.162/(0.1847*s+1);
L_K = Cc * P;
sisotool(L_K)