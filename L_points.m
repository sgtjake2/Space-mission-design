close all
clear
tic
G = 6.674e-11; %kg^-1m^3s^-2

%{
i = 1.727149683568981E-01* pi/180;
e = 0.055;
a = 0.384e9;
P = 0.363e9;
w = 1.452529123375567E+02*(pi/180);
f = 1:1:360;
f = f' * (pi/180);
omega = 7.110013372809593E+01*(pi/180);
tau = 2.460062781279559E+06;
T = 27.321582;
%}

a = sqrt(((-2.4671e8)^2)+((3.1187e8)^2));

mE = 5.97219e24;
mM = 7.349e22;

alpha = mM/(mM+mE);

R = 0.384e9;

ohm = sqrt((G*(mM+mE))/(R^3));
V = [0 0 1 0;
    0 0 0 1;
    9*ohm^2 0 0 2*ohm;
    0 -3*ohm^2 2*ohm 0];

vs = sqrt(G*1000/35000000);

L1 = a*(1-((alpha/3)^(1/3)));
L2 = a*(1+((alpha/3)^(1/3)));
L3 = -a*(1-((5*alpha/12)^(1/3)));
L4 = (a/2)*((mE-mM)/(mE+mM));
L5 = (a/2)*((mE-mM)/(mE+mM));

L2 = 4.475451e8;

theta = -0.7375;
L2x = L2*sin(theta);
L2y = L2 *cos(theta);
L2z = L2 * sin(5.145*pi/180);

%Vectors for Earth, Moon and Satalite
posE = [3.034529801433376E+03 -3.340359653273315E+03 -1.939037192252903E+03]*1000;
posM = [-2.467089971416006E+05 3.118709928882247E+05 3.661070919281326E+04]*1000;
%posS = [L2x L2y 1.576448255327776E+08];

%posE = [3.034529801433376E+03 -3.340359653273315E+03]*1000;
%posM = [-2.467089971416006E+05 3.118709928882247E+05]*1000;
%posS = [L2x L2y L2z];
posS = [-2.066118376482230E+05 4.155370802455372E+05 4.252060750636854E+04+3500]*1000;
%posS = [-2.588076140146283E+05 3.018745235389227E+05]*1000;

velE = [9.245683416082787E-03 6.608382937807280E-03 2.954691314042687E-03]*1000;
velM = [-7.516793153236040E-01 -5.884845299522508E-01 -6.683923719039309E-03]*1000;
%velS = [-1 -1 -1]*642;

%velE = [9.245683416082787E-03 6.608382937807280E-03]*1000;
%velM = [-7.516793153236040E-01 -5.884845299522508E-01]*1000;
%velS = [-1.277 -1]*700;
%velS = (velM.*posS)./posM
%velS = [1.257*velM(1) 1.094*velM(2) 1*velM(3)]+50;
%velS = [-1.012041739961583E+00 -4.800076965180470E-01 1.499994544138000E-02]*1000;
%velS = 0;
velS = velM*([-2.86e8 3.616e8 6.742e7]/posM);
%velS = [-7.253275256058060E-01 -2.281755498821680E-01]*1000;

%Simulink Output
out = sim("Earth_Moon.slx");
OE = out.OE;
OM = out.OM;
OS = out.OS;

data = xlsread("Earth-Moon-2023 (2).csv");

figure(2)
hold on
grid on
OMp = plot3(OM.Data(:,1),OM.Data(:,2),OM.Data(:,3));
OEp = plot3(OE.Data(:,1),OE.Data(:,2),OE.Data(:,3), "Marker","*");
p = plot3(OM.Data(1,1),OM.Data(1,2),OM.Data(1,3),"Marker","*");
OMs = plot3(OS.Data(:,1),OS.Data(:,2),OS.Data(:,3));
%L2p = plot(L2x,L2y,"Marker","*");
dat = plot3(data(:,3)*1000,data(:,4)*1000, data(:,5)*1000);
view(3)

h = animatedline("Color","red","LineStyle","--","Marker","o");
hS = animatedline("Color","red","LineStyle","--");
hM = animatedline("Color","green","LineStyle","--","Marker","o");

daspect([1 1 1])

for i = 1:1000:length(OM.Data)
    addpoints(hM,OM.Data(i,1),OM.Data(i,2),OM.Data(i,3))
    addpoints(h,OS.Data(i,1),OS.Data(i,2),OS.Data(i,3))
    addpoints(hS,OS.Data(i,1),OS.Data(i,2),OS.Data(i,3))

    pause(0.001)

    clearpoints(h)
    clearpoints(hM)
end

toc