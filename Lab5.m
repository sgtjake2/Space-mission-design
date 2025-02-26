clear
close all

syms r 

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

mu = mM/mE;

posE = [-xE 0 0];
posM = [xM 0 0];

velE = [0 -vE 0];
velM = [0 vM 0];

alpha = mM/(mM+mE);

%L2 = xM*(1+((alpha/3)^(1/3)));
%L2 = xM*(1+((mM/(3*mE))^(1/3)));
%L2 = D*(1+((alpha/3)^(1/3)));
%L2 = 4.475451e+08;
%r = 3000*1000;

eq = mE/((D+r)^2) + (mM/(r^2)) == (mE/(D^2)) + r*((mM+mE)/(D^3));
sol = solve(eq,r);
sol = vpa(sol);

%xDash = x-(1-(mM/mE)+(hL2-xM));

hL2 = double(sol(3))+xM;

%vSo = sqrt((G*(mE))/L2);
%vSe = o*L2;
vSe = o*(hL2);

rho = D*((mM/mE)^(2/5));

velS = [0 vSe 0];
posS = [hL2 0 0];

%Az = 0.2*xM;
%velS = [0 sqrt(o^2 * Az^2) 0];
%posS = [hL2 0 Az];

out = sim("Lab5_Sim");

OE = out.OE.data;
OM = out.OM.data;
OS = out.OS.data;

figure(1)
hold on
grid on

e = plot3(OE(:,1),OE(:,2),OE(:,3));
m = plot3(OM(:,1),OM(:,2),OM(:,3));
s = plot3(OS(:,1),OS(:,2),OS(:,3));

view(3)
daspect([1 1 1])

hM = animatedline("Color","red","LineStyle","--","Marker","o");
hS = animatedline("Color","green","LineStyle","--","Marker","o");
hE = animatedline("Color","blue","LineStyle","--","Marker","o");
hB = animatedline("Color","black","LineStyle","--");
hC = animatedline("Color","black","LineStyle","-.");
hD = animatedline("Color","black","LineStyle","-.");

for i = 1:100:length(OM)
    addpoints(hM,OM(i,1),OM(i,2),OM(i,3))
    addpoints(hE,OE(i,1),OE(i,2),OE(i,3))
    addpoints(hS,OS(i,1),OS(i,2),OS(i,3))
    b = [OM(i,1) OM(i,2) OM(i,3);
        OE(i,1) OE(i,2) OE(i,3)];
    addpoints(hB, b(:,1),b(:,2),b(:,3))
    c = [OS(i,1) OS(i,2) OS(i,3);
        OE(i,1) OE(i,2) OE(i,3)];
    addpoints(hC, c(:,1),c(:,2),c(:,3))
    d = [OS(i,1) OS(i,2) OS(i,3);
        OM(i,1) OM(i,2) OM(i,3)];
    addpoints(hD, d(:,1),d(:,2),d(:,3))

    pause(0.001)

    if i < length(OM)-100
        clearpoints(hM)
        clearpoints(hE)
        clearpoints(hB)
        clearpoints(hS)
        clearpoints(hC)
        clearpoints(hD)
    end
end

%% huh

%{
clear
close all

k = 0.1;
lambda = sqrt(0.7);
v =0.1;
t=1:100;

Ay = 3000;
Az = 3000;

x = -k*Ay *cos(lambda*t);
y = Ay*sin(lambda*t);
z = Az * sin(v*t+(pi/2));

figure(333)
hold on
grid on
xlabel("x")
ylabel("y")
zlabel("z")

plot3(x,y,z)
view(3)
daspect([1 1 1])

figure(33)
hold on
grid on
xlabel("y")
ylabel("z")

plot(y,z)
daspect([1 1 1])

figure(34)
hold on
grid on
xlabel("t")
ylabel("x,y,z")

plot(t,x)
plot(t,y)
plot(t,z)
%}