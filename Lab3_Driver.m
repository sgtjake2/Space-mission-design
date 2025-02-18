close all
clear

%Constants for Model

G = 6.674e-11; %kg^-1m^3s^-2
Ms = 1988410e24; %kg
Mp = 5.97219e24; %kg
Mm = 6.4171e23;

Vel_s = [4.508e-3 -1.335e-2 -5.769e-3]*1000; %m/s
Vel_p = [2.749848392692773E+01 -1.032646088959717E+01 1.535581598741675E-03]*1000; %m/s
Vel_v = [-2.625497707637767E+01 1.649830245531765E+01 5.627457721235540E-02]*1000; %m/s
Vel_m = [-8.699321075487822E+00 -2.018860873085168E+01 -2.092884078165804E-01]*1000;

Pos_s = [-1.313e6 -1.803e5 -4.319e4]*1000; %m
Pos_p = [-5.328792362541436E+07 -1.426707827423932E+08 3.981048565176874E+04]*1000; %m
Pos_m = [-2.299762717931310E+08 9.863481268283156E+07 7.712126942588501E+06]*1000;
Pos_v = [-5.362193933254258E+07 -1.428620826284412E+08 4.195466200186312E+04]*1000;

out = sim("Lab3_Sim.slx");

%% out

%Earth Simulink Outputs
XE = out.OP.Data(:,1);
YE = out.OP.Data(:,2);
ZE = out.OP.Data(:,3);

%Sun Simulink Outputs
XS = out.OS.Data(:,1);
YS = out.OS.Data(:,2);
ZS = out.OS.Data(:,3);

%OM Simulink Outputs
XM = out.OM.Data(:,1);
YM = out.OM.Data(:,2);
ZM = out.OM.Data(:,3);

%OV Simulink Outputs
XV = out.OV.Data(:,1);
YV = out.OV.Data(:,2);
ZV = out.OV.Data(:,3);

%NASA Refrence Data
%data = xlsread("Sun-Earth-2023-2024-equotoarl(1).csv")*1000;
%data_mars = xlsread("Sun-Mars-2023-2025-ecliptic.csv")*1000;
%x = data_mars(:,3);
%y = data_mars(:,4);
%z = data_mars(:,5);


figure(1)
hold on
E = plot3(XE,YE,ZE);
grid on
S = plot3(XS,YS,ZS,"*");
%plot3(x,y,z)
M = plot3(XM,YM,ZM);
%V = plot3(XV,YV,ZV);

%p = [E,S,M];
%l = ["Earth","Sun","Mars"];
%legend(p,l)

daspect([1 1 1])
view(3)


h = animatedline("Color","red","LineStyle","--","Marker","o");
h1 = animatedline("Color","red","LineStyle","--");
for i = 1:5:length(XE)
    addpoints(h,XV(i),YV(i),ZV(i))
    addpoints(h1,XV(i),YV(i),ZV(i))
    pause(0.001)

    clearpoints(h)
end

%{
for i = 1:length(XE)
    h =  animatedline;
    hE = addpoints(h,XE(i),YE(i),ZE(i),"Marker","o","Color","blue");
    hM = addpoints(h,XM(i),YM(i),ZM(i),"Marker","o","Color","red");
    hV = addpoints(h,XV(i),YV(i),ZV(i),"Marker","o","Color","green");
    pause(0.01)
end
%}

