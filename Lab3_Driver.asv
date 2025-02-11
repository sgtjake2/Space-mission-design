close all

%Constants for Model

G = 6.674e-11; %kg^-1m^3s^-2
Ms = 1988410e24; %kg
Mp = 5.97219e24; %kg

Vel_s = [4.508e-3 -1.335e-2 -5.769e-3]*1000; %m/s
Vel_p = [2.749e1 -9.479 -4.106]*1000; %m/s

Pos_s = [-1.313e6 -1.803e5 -4.319e4]*1000; %m
Pos_p = [-5.328e7 -1.309e8 -5.671e7]*1000; %m

%Earth Simulink Outputs
XE = out.simout1.Data(:,1);
YE = out.simout1.Data(:,2);
ZE = out.simout1.Data(:,3);

%Sun Simulink Outputs
XS = out.simout.Data(:,1);
YS = out.simout.Data(:,2);
ZS = out.simout.Data(:,3);

%NASA Refrence Data
data = xlsread("Sun-Earth-2023-2024-equotoarl(1).csv")*1000;
x = data(:,3);
y = data(:,4);
z = data(:,5);

figure(1)
plot3(XE,YE,ZE)
hold on
grid on
plot3(XS,YS,ZS,"*")
plot3(x,y,z)
daspect([1 1 1])
