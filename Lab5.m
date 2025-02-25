clear
close all

xE = 4.628140551638176E+03*1000; %km
xM = 3.762704566571173E+05*1000; %km

D = xE + xM;

mE = 5.97219e24;
mM = 7.349e22;

G = 6.674e-11; %kg^-1m^3s^-2

o = sqrt((G*(mE+mM))/(D^3));

T = 2*pi*sqrt((D^3)/(G*(mE+mM)));

vE = (2*pi*xE)/T;
vM = (2*pi*xM)/T;

vE = o*xE;
vM = o*xM;

posE = [-xE 0 0];
posM = [xM 0 0];

velE = [0 -vE 0];
velM = [0 vM 0];

out = sim("Lab5_Sim");

OE = out.OE.data;
OM = out.OM.data;

figure(1)
hold on
grid on

e = plot3(OE(:,1),OE(:,2),OE(:,3),"Marker","*");
m = plot3(OM(:,1),OM(:,2),OM(:,3));

daspect([1 1 1])

alpha = mM/(mM+mE);

L2 = xM*(1+((alpha/3)^(1/3)));