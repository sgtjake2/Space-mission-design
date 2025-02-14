close all

%Constants for Model

G = 6.674e-11; %kg^-1m^3s^-2
Ms = 1988410e24; %kg
Mp = 5.97219e24; %kg
Mm = 6.4171e23;

Vel_s = [4.508e-3 -1.335e-2 -5.769e-3]*1000; %m/s
Vel_p = [2.749848392692773E+01 -1.032646088959717E+01 1.535581598741675E-03]*1000; %m/s
Vel_v = [1 1 1]*1000; %m/s
Vel_m = [-8.699321075487822E+00 -2.018860873085168E+01 -2.092884078165804E-01]*1000;

Pos_s = [-1.313e6 -1.803e5 -4.319e4]*1000; %m
Pos_p = [-5.328792362541436E+07 -1.426707827423932E+08 3.981048565176874E+04]*1000; %m
Pos_m = [-2.299762717931310E+08 9.863481268283156E+07 7.712126942588501E+06]*1000;

%% out

%Earth Simulink Outputs
XE = out.simout1.Data(:,1);
YE = out.simout1.Data(:,2);
ZE = out.simout1.Data(:,3);

%Sun Simulink Outputs
XS = out.simout.Data(:,1);
YS = out.simout.Data(:,2);
ZS = out.simout.Data(:,3);

%OVs Simulink Outputs
XVm = out.simout3.Data(:,1);
YVm = out.simout3.Data(:,2);
ZVm = out.simout3.Data(:,3);

%OVp Simulink Outputs
XVp = out.simout2.Data(:,1);
YVp = out.simout2.Data(:,2);
ZVp = out.simout2.Data(:,3);

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
%plot3(x,y,z)
plot3(XVm,YVm,ZVm)
%plot3(XVp,YVp,ZVp)
daspect([1 1 1])
