clear
close all

syms r 

xE = 4.628140551638176E+03*1000; %m
xM = 3.762704566571173E+05*1000; %m

D = xE + xM;

mE = 5.97219e24; %kg
mM = 7.349e22; %kg
mS = 1000; %kg

G = 6.674e-11; %kg^-1m^3s^-2

o = sqrt((G*(mE+mM))/(D^3)); %rad/s

T = 2*pi*sqrt((D^3)/(G*(mE+mM))); %seconds

vE = o*xE; %m/s
vM = o*xM; %m/s

%mu = mM/mE;

posE = [-xE 0 0];
posM = [xM 0 0];

velE = [0 -vE 0];
velM = [0 vM 0];

%Circular Orbit
%Finding L2 Point
eq = mE/((D+r)^2) + (mM/(r^2)) == (mE/(D^2)) + r*((mM+mE)/(D^3));
sol = solve(eq,r);
sol = vpa(sol);
xS = double(sol(3));

xL2 = xS+xM;

vSe = o*(xL2);

rho = D*((mM/mE)^(2/5));

%Satalite Vectors
velS = [0 vSe 0];
posS = [xL2 0 0];

%Outputting Simulink Data
out = sim("Lab5_Sim");

OE = out.OE.data;
OM = out.OM.data;
OMv = out.OMv.data;
OS = out.OS.data;

%find h and check if h is constant
h = cross(OM,OMv);
h1 = cross(OM(1,:),OMv(1,:));
h2 = cross(OM((47174),:),OMv((47174),:));
h3 = cross(OM((length(OE)),:),OMv((length(OE)),:));

h_hat = zeros(length(OM),3);
r_hat = zeros(length(OM),3);
u_hat = zeros(length(OM),3);

for i = 1:length(OM)
    h_hat_temp = h(i,:)./norm(h(i,:));
    r_hat_temp = OM(i,:)/norm(OM(i,:));
    u_hat_temp = cross(h_hat_temp,r_hat_temp);

    h_hat(i,:) = h_hat_temp;
    r_hat(i,:) = r_hat_temp;
    u_hat(i,:) = u_hat_temp;
end

xr = zeros(length(OM),3);
yr = zeros(length(OM),3);
zr = zeros(length(OM),3);

for j = 1:length(OM)
    xr_temp = dot(OS(j,:),r_hat(j,:));
    yr_temp= dot(OS(j,:),u_hat(j,:));
    zr_temp = dot(OS(j,:),h_hat(j,:));

    xr(j,:) = xr_temp;
    yr(j,:) = yr_temp;
    zr(j,:) = zr_temp;
end

%Plotting Rotational Plane
figure(2)
movegui([1000 300])
subplot(1,2,1)
rotational = plot(yr(:,2),zr(:,3));
hold on
grid on
xlabel("Yr (m)")
ylabel("Zr (m)")
title(["Satalite Path in Rotational Plane",""])
start = plot(yr(1,2),zr(1,3),"Marker","*","Color","green","LineStyle","none");
finish = plot(yr(end,2),zr(end,3),"Marker","*","Color","red","LineStyle","none");
moon = plot(0,0,"MarkerSize",15,"Marker","o","Color","red","LineStyle","none");
hR1 = animatedline("Color","green","LineStyle","none","Marker","o"); %Satalite Marker on rotational frame

xlim([-2000 2000])
ylim([-2000 2000])

plots = [rotational start finish moon hR1];
l = ["Satalite Path","Start Position", "End Posistion","Moon","Satalite"];

legend(plots, l)

hold off

subplot(1,2,2)
rotational = plot(yr(:,2),zr(:,3));
hold on
grid on
xlabel("Yr (m)")
ylabel("Zr (m)")
title(["Satalite Path in Rotational Plane",""])
start = plot(yr(1,2),zr(1,3),"Marker","*","Color","green","LineStyle","none");
finish = plot(yr(end,2),zr(end,3),"Marker","*","Color","red","LineStyle","none");
moon = plot(0,0,"MarkerSize",15,"Marker","o","Color","red","LineStyle","none");
hR2 = animatedline("Color","green","LineStyle","none","Marker","o"); %Satalite Marker on rotational frame

hold off

%Plotting three body simulation
figure(1)
movegui([300 300])
hold on
grid on
xlabel("X (m)")
ylabel("Y (m)")
zlabel("Z (m)")
title("Three Body Simulation of Earth, Moon and Satalite")
mid = int16(length(OS)/2);
e = plot3(OE(:,1),OE(:,2),OE(:,3));
m = plot3(OM(:,1),OM(:,2),OM(:,3));
s = plot3(OS(:,1),OS(:,2),OS(:,3));
start = plot3(OS(1,1),OS(1,2),OS(1,3),"Marker","*","Color","green","LineStyle","none");
finish = plot3(OS(end,1),OS(end,2),OS(end,3),"Marker","*","Color","red","LineStyle","none");

view(3)
%daspect([1 1 1])

hM = animatedline("Color","red","LineStyle","none","Marker","o"); %Moon Marker
hS = animatedline("Color","green","LineStyle","none","Marker","o"); %Satalite Marker
hE = animatedline("Color","blue","LineStyle","--","Marker","o"); %Earth Marker
hB = animatedline("Color","black","LineStyle","--"); %Moon to Earth Line
hC = animatedline("Color","black","LineStyle","-."); %Satalite to Earth Line
hD = animatedline("Color","black","LineStyle","-."); %Moon to Satalite Line

plots = [e m hM s hS start finish];
l = ["Earth","Moon", "","Satalite","","Start Position","End Position"];

legend(plots, l)

hold off

for i = 1:100:length(OM)
    %figure(1)
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

    %figure(2)
    addpoints(hR1,yr(i,2),zr(i,3))
    addpoints(hR2,yr(i,2),zr(i,3))


    pause(0.001)

    if i < length(OM)-100
        clearpoints(hM)
        clearpoints(hE)
        clearpoints(hB)
        clearpoints(hS)
        clearpoints(hC)
        clearpoints(hD)
        clearpoints(hR1)
        clearpoints(hR2)
    end
end

%% hmmmmm

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