%%
% 1 Lunar Orbit  = 2358720 seconds

clear
close all

syms r 

%% Define Parameters
xE = 4.628140551638176E+03*1000; % Earth to barycenter (m)
xM = 3.762704566571173E+05*1000; % Moon to barycenter (m)

D = xE + xM; % Distance from earth to moon (m)

mE = 5.97219e24; % Earth Mass (kg)
mM = 7.349e22; % Moon Mass (kg)

G = 6.674e-11; % Gravitational constant (kg^-1m^3s^-2)

o = sqrt((G*(mE+mM))/(D^3)); % Angular velocity of system (rad/s)

vE = o*xE; % Initial velcity of Earth (m/s)
vM = o*xM; % Initial velcoity of Moon (m/s)

posE = [-xE 0 0]; %Initial position vector of earth
posM = [xM 0 0]; %Initial position vector of moon

velE = [0 -vE 0]; %Initial velocity vector of earth
velM = [0 vM 0]; %Initial velcoity vector of moon


%Finding L2 Point
eq = mE/((D+r)^2) + (mM/(r^2)) == (mE/(D^2)) + r*((mM+mE)/(D^3));
sol = solve(eq,r); %Solve for r
sol = vpa(sol);
xS = double(sol(3)); %Third solution is the requried r

xL2 = xS+xM; %Distance from barycenter to L2 (m)

vSe = o*(xL2); % Initial velocity of spacecraft (m/s)

% Initial Satalite Vectors for halo orbit (radius of 2000km)
% Trial and error method
velS = [-0.00103491*2222.6 vSe+(0.004*2000) 0];
posS = [xL2-(200*2000) 0 1000*2000];

%% Simulation
%Outputting Simulink Data
out = sim("CW_Sim");

%Simulink vectors
OE = out.OE.data;
OM = out.OM.data;
OMv = out.OMv.data;
OS = out.OS.data;

%% Creating rotational frame
%Find angular momentum
%(h sholud be constant throughout whole sim)
h = cross(OM,OMv); 

%Converting from inertial to rotational frame
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

%% Plotting
%Position of middle point of simulation
mid = int16(length(OS)/2);

%Figure 1
%Plotting three body simulation in inertial frame
f1 = figure(1);
movegui([0 550])
hold on
grid on
xlabel("X (km)")
ylabel("Y (km)")
zlabel("Z (km)")
title("Three Body Simulation of Earth, Moon and Satalite")
e = plot3(OE(:,1)/1000,OE(:,2)/1000,OE(:,3)/1000);
m = plot3(OM(:,1)/1000,OM(:,2)/1000,OM(:,3)/1000);
s = plot3(OS(:,1)/1000,OS(:,2)/1000,OS(:,3)/1000);
start = plot3(OS(1,1)/1000,OS(1,2)/1000,OS(1,3)/1000,"Marker","*","Color","green","LineStyle","none");
finish = plot3(OS(end,1)/1000,OS(end,2)/1000,OS(end,3)/1000,"Marker","*","Color","red","LineStyle","none");

view(3)

hM = animatedline("Color","red","LineStyle","none","Marker","o"); %Moon Marker
hS = animatedline("Color","green","LineStyle","none","Marker","o"); %Satalite Marker
hE = animatedline("Color","blue","LineStyle","none","Marker","o"); %Earth Marker
hB = animatedline("Color","black","LineStyle","--"); %Moon to Earth Line
hC = animatedline("Color","black","LineStyle","-."); %Satalite to Earth Line
hD = animatedline("Color","black","LineStyle","-."); %Moon to Satalite Line

plots = [hE m hM s hS start finish];
l = ["Earth","Moon Path", "Moon","Satalite Path","Satalite","Start Position","End Position"];

lgd = legend(plots, l);
lgd.Location = "eastoutside";

hold off

%Figure 2
%Plotting Rotational Plane (Y and Z axes)
f2 = figure(2);
movegui([550 550])
rotational = plot(yr(:,2)/1000,zr(:,3)/1000);
hold on
grid on
xlabel("Yr (km)")
ylabel("Zr (km)")
title(["Satalite Path in Rotational Plane",""])
start = plot(yr(1,2)/1000,zr(1,3)/1000,"Marker","*","Color","green","LineStyle","none");
finish = plot(yr(end,2)/1000,zr(end,3)/1000,"Marker","*","Color","red","LineStyle","none");
moon = plot(0,0,"MarkerSize",15,"Marker","o","Color","red","LineStyle","none");

hR2 = animatedline("Color","green","LineStyle","none","Marker","o"); %Satalite Marker on rotational frame (Front)

ylim([-2500 2100])

daspect([1 1 1])

plots = [rotational start finish moon hR2];
l = ["Satalite Path","Start Position", "End Posistion","Moon","Satalite"];

lgd = legend(plots, l);
lgd.Location = "eastoutside";

hold off


%Figure 3
%Plotting Rotational Plane (X and Z axes)
f3 = figure(3);
movegui([1100 550])
hold on
grid on
side = plot(xr(:,1)/1000,zr(:,3)/1000);
start = plot(xr(1,2)/1000,zr(1,3)/1000,"Marker","*","Color","green","LineStyle","none");
finish = plot(xr(end,2)/1000,zr(end,3)/1000,"Marker","*","Color","red","LineStyle","none");
l2 = plot(xL2/1000,0,"MarkerSize",10,"Marker","o","Color","red","LineStyle","none");
xlabel("Xr (km)")
ylabel("Zr (km)")
title(["Side view of the Satalite Path"; "in the Rotational Frame"])

hR3 = animatedline("Color","green","LineStyle","none","Marker","o");% Satalite Marker on rotational frame (Side)

ylim([-2500 2100])
xlim([4.39e5 4.41e5])

daspect([1 1 1])

plots = [side start finish l2 hR3];
l = ["Satalite Path","Start Position", "End Posistion","L2 Point","Satalite"];

lgd = legend(plots, l);
lgd.Location = "eastoutside";

hold off

%Figure 5
%Unanimated version of Figure 1
f5 = figure(5);
movegui([200 50])
hold on
grid on
xlabel("X (km)")
ylabel("Y (km)")
zlabel("Z (km)")
title("Three Body Simulation of Earth, Moon and Satalite")
e = plot3(OE(1,1)/1000,OE(1,2)/1000,OE(1,3)/1000,Marker="o",LineStyle="none",Color="blue");
m = plot3(OM(:,1)/1000,OM(:,2)/1000,OM(:,3)/1000);
s = plot3(OS(:,1)/1000,OS(:,2)/1000,OS(:,3)/1000);
mM = plot3(OM(mid,1)/1000,OM(mid,2)/1000,OM(mid,3)/1000,Marker="o",LineStyle="none",Color="red");
mS = plot3(OS(mid,1)/1000,OS(mid,2)/1000,OS(mid,3)/1000,Marker="o",LineStyle="none",Color="green");
start = plot3(OS(1,1)/1000,OS(1,2)/1000,OS(1,3)/1000,"Marker","*","Color","green","LineStyle","none");
finish = plot3(OS(end,1)/1000,OS(end,2)/1000,OS(end,3)/1000,"Marker","*","Color","red","LineStyle","none");

view(3)

plots = [e m mM s mS start finish];
l = ["Earth","Moon Path", "Moon","Satalite Path","Satalite","Start Position","End Position"];

lgd = legend(plots, l);
lgd.Location = "eastoutside";

hold off


%Figure 6
%Plotting rotational frame (3D)
f6 = figure(6);
movegui([750 50])
hold on
grid on
xlabel("Xr (km)")
ylabel("Yr (km)")
zlabel("Zr (km)")
title("Three Body Simulation of Earth, Moon and Satalite")
s = plot3(xr(:,1)/1000,yr(:,2)/1000,zr(:,3)/1000);
l2 = plot3(xL2/1000,0,0, LineStyle="none",Marker="o");
start = plot3(xr(1,1)/1000,yr(1,2)/1000,zr(1,3)/1000,"Marker","*","Color","green","LineStyle","none");
finish = plot3(xr(end,1)/1000,yr(end,2)/1000,zr(end,3)/1000,"Marker","*","Color","red","LineStyle","none");

hR4 = animatedline("Color","green","LineStyle","none","Marker","o"); %Satalite Marker on rotational frame

plots = [s l2 start finish hR4];
l = ["Satalite Path","L2 Point", "Start Position","End Position","Satalite"];

lgd = legend(plots, l);
lgd.Location = "eastoutside";

view(3)
daspect([1 1 1])

hold off


%Animation Loop
for i = 1:100:length(OM)
    %Adding points
    %Figure 1
    addpoints(hM,OM(i,1)/1000,OM(i,2)/1000,OM(i,3)/1000)
    addpoints(hE,OE(i,1)/1000,OE(i,2)/1000,OE(i,3)/1000)
    addpoints(hS,OS(i,1)/1000,OS(i,2)/1000,OS(i,3)/1000)
    b = [OM(i,1) OM(i,2) OM(i,3);
        OE(i,1) OE(i,2) OE(i,3)];
    addpoints(hB, b(:,1)/1000,b(:,2)/1000,b(:,3)/1000)
    c = [OS(i,1) OS(i,2) OS(i,3);
        OE(i,1) OE(i,2) OE(i,3)];
    addpoints(hC, c(:,1)/1000,c(:,2)/1000,c(:,3)/1000)
    d = [OS(i,1) OS(i,2) OS(i,3);
        OM(i,1) OM(i,2) OM(i,3)];
    addpoints(hD, d(:,1)/1000,d(:,2)/1000,d(:,3)/1000)

    %Figure 2
    addpoints(hR2,yr(i,2)/1000,zr(i,3)/1000)

    %Figure 3
    addpoints(hR3,xr(i,1)/1000,zr(i,3)/1000)

    %Figure 6
    addpoints(hR4,xr(i,1)/1000,yr(i,2)/1000,zr(i,3)/1000)

    pause(0.001)

    %Deleting Points
    if i < length(OM)-100
        clearpoints(hM)
        clearpoints(hE)
        clearpoints(hB)
        clearpoints(hS)
        clearpoints(hC)
        clearpoints(hD)
        clearpoints(hR2)
        clearpoints(hR3)
        clearpoints(hR4)
    end
end