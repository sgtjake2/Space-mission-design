close all
clear

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

L1 = a*(1-((alpha/3)^(1/3)));
L2 = a*(1+((alpha/3)^(1/3)));
L3 = -a*(1-((5*alpha/12)^(1/3)));
L4 = (a/2)*((mE-mM)/(mE+mM));
L5 = (a/2)*((mE-mM)/(mE+mM));

theta = -0.7375;
L2x = L2*sin(theta);
L2y = L2 *cos(theta);

%Vectors for Earth, Moon and Satalite
%posE = [3.034529801433376E+03 -3.340359653273315E+03 -1.939037192252903E+03]*1000;
%posM = [-2.467089971416006E+05 3.118709928882247E+05 3.661070919281326E+04]*1000;
%posS = [L2x L2y 1.576448255327776E+08];

posE = [3.034529801433376E+03 -3.340359653273315E+03]*1000;
posM = [-2.467089971416006E+05 3.118709928882247E+05]*1000;
%posS = [L2x L2y];
posS = [-2.588076140146283E+05 3.018745235389227E+05]*1000;

%velE = [9.245683416082787E-03 6.608382937807280E-03 2.954691314042687E-03]*1000;
%velM = [-7.516793153236040E-01 -5.884845299522508E-01 -6.683923719039309E-03]*1000;
%velS = [-1 -1 -1]*642;

velE = [9.245683416082787E-03 6.608382937807280E-03]*1000;
velM = [-7.516793153236040E-01 -5.884845299522508E-01]*1000;
%velS = [-1.277 -1]*700;
%velS = (velM.*posS)./posM
%velS = [1.257*velM(1) 1.094*velM(2)]+100;
velS = [-7.253275256058060E-01 -2.281755498821680E-01]*1000;

%Simulink Output
out = sim("Earth_Moon.slx");
OE = out.OE;
OM = out.OM;
OS = out.OS;

%data = xlsread("Data\Earth-Moon-2023.csv");

%{
L = a.*(1-(e.^2));
r = L./(1+e.*cos(f));

x2 = r.*cos(f+w);
y2 = r.*sin(f+w);

x0 = x2.*cos(omega)-y2.*cos(i).*sin(omega);
y0 = x2.*sin(omega) + y2.*cos(i).*cos(omega);
z0 = y2.*sin(i);
%}

figure(1)
%p1 = plot(x0,y0);
hold on
grid on

%x = P.*cos(w);
%y = P.*sin(w);

%{
x1 = x.*cos(omega)-y.*cos(i).*sin(omega);
y1 = x.*sin(omega) + y.*cos(i).*cos(omega);
z1 = y.*sin(i);
p2 = plot(x1,y1,'*','Color','blue','MarkerSize',10,...
    'MarkerFaceColor','blue');
%}

%{
L1p=plot(L1,0,"Marker","x");
L2p=plot(L2,0,"Marker","x");
L3p=plot(L3,0,"Marker","x");
L4p=plot(L4,((sqrt(3)/2)*a),"Marker","x");
L5p=plot(L5,-((sqrt(3)/2)*a),"Marker","x");
%legend("Moon","L1","L2","L3","L4","L5")
%}

OMp = plot(OM.Data(:,1),OM.Data(:,2));
OEp = plot(OE.Data(:,1),OE.Data(:,2), "Marker","*");
p = plot(OM.Data(1,1),OM.Data(1,2),"Marker","*");
%OMs = plot(OS.Data(:,1),OS.Data(:,2));
%L2p = plot(L2x,L2y,"Marker","*");

daspect([1 1 1])
h = animatedline("Color","red","LineStyle","--","Marker","o");
hS = animatedline("Color","red","LineStyle","--");
hM = animatedline("Color","green","LineStyle","--","Marker","o");
for i = 1:1:length(OM.Data)
    addpoints(hM,OM.Data(i,1),OM.Data(i,2))
    addpoints(h,OS.Data(i,1),OS.Data(i,2))
    addpoints(hS,OS.Data(i,1),OS.Data(i,2))

    pause(0.01)

    clearpoints(h)
    clearpoints(hM)
end

%dat = plot(data(:,3)*1000,data(:,4)*1000,"Marker","+");

%{
figure(2)
hold on
grid on
OMp = plot3(OM.Data(:,1),OM.Data(:,2),OM.Data(:,3));
OEp = plot3(OE.Data(:,1),OE.Data(:,2),OE.Data(:,3), "Marker","*");
p = plot3(OM.Data(1,1),OM.Data(1,2),OM.Data(1,3),"Marker","*");
%OMs = plot3(OS.Data(:,1),OS.Data(:,2),OS.Data(:,3));
%L2p = plot(L2x,L2y,"Marker","*");
view(3)
h = animatedline("Color","red","LineStyle","--","Marker","o");
hS = animatedline("Color","red","LineStyle","--");
hM = animatedline("Color","green","LineStyle","--","Marker","o");
for i = 1:1:length(OM.Data)
    addpoints(hM,OM.Data(i,1),OM.Data(i,2),OM.Data(i,3))
    addpoints(h,OS.Data(i,1),OS.Data(i,2),OS.Data(i,3))
    addpoints(hS,OS.Data(i,1),OS.Data(i,2),OS.Data(i,3))

    pause(1)

    clearpoints(h)
    clearpoints(hM)
end

daspect([1 1 1])
%}


