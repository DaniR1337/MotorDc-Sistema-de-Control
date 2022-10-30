clc; clear all; close all;
% Identificación Toolbox 
% Integrantes: Daniel Alejandro Rodríguez Alvarado, Nataly Delgado Huertas y Sylvia Fonseca Cruz
% Carné: C06575, C02583 y C03039

datos = readmatrix('delta_85a115.csv');
t=datos(:,1);
u=datos(:,2);
y=datos(:,3);
figure(1)
%plot(t,u,t,y)
u_1=u-85;
y_1=y-107.2;
figure(2)
plot(t,u_1,t,y_1)

systemIdentification



