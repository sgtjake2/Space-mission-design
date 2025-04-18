clear, close all

i = [7 3.4 4.13e-3 1.8]*(pi/180);
e = [0.206 0.007 0.017 0.094];
P = [46 107.5 147.1 206.7];
a = [57.9 108.2 149.6 228];
w = [29.2 54.8 256 287]*(pi/180);
f = 1:1:360;
f = f' * (pi/180);
omega = [48.3 76.7 205 49.5]*(pi/180);
tau = [2460123 2460052 2459946 2460438];
T = [88 224.7 365.2 687];

t = "2002-10-08";
t = datetime(t);

jd = juliandate(t);

M = 2*pi*((jd-tau)./T);
E = [3 3 3 3];

for n = 1:100
    E=M+e.*sin(E);
end

cosTheta =(cos(E)-e)./(1-e.*cos(E));
theta = acos(cosTheta);

L = a.*(1-(e.^2));
r = L./(1+e.*cos(f));

x2 = r.*cos(f+w);
y2 = r.*sin(f+w);

x0 = x2.*cos(omega)-y2.*cos(i).*sin(omega);
y0 = x2.*sin(omega) + y2.*cos(i).*cos(omega);
z0 = y2.*sin(i);

figure(1)
p1 = plot3(x0,y0,z0);
hold on
grid on

[X,Y,Z] = sphere;
circle = 10;
X2 = X * circle;
Y2 = Y * circle;
Z2 = Z * circle;

p3 = surf(X2+5,Y2-5,Z2,FaceColor="yellow");

xlabel("X (km *10^6)")
ylabel("Y (km *10^6)")
zlabel("Z (km *10^6)")
title("Orbital Elements Plot")
daspect([1 1 1])

%f at periapsis
x = P.*cos(w);
y = P.*sin(w);

x1 = x.*cos(omega)-y.*cos(i).*sin(omega);
y1 = x.*sin(omega) + y.*cos(i).*cos(omega);
z1 = y.*sin(i);
p2 = plot3(x1,y1,z1,'*','Color','blue','MarkerSize',10,...
    'MarkerFaceColor','blue');

r3 = L./(1+e.*cos(theta));

x3 = r3.*cos(theta+w);
y3 = r3.*sin(theta+w);

x03 = x3.*cos(omega)-y3.*cos(i).*sin(omega);
y03 = x3.*sin(omega) + y3.*cos(i).*cos(omega);
z03 = y3.*sin(i);

p4 = plot3(x03,y03,z03,'o','Color','red','MarkerSize',10,...
    'MarkerFaceColor','red');

data = xlsread("Sun-Earth-2023-2025(1).csv");
dat = plot3(data(:,3)*1e-6,data(:,4)*1e-6, data(:,5)*1e-6);

l = ["Mercury","Venus","Earth","Mars","Periapsis","Sun","Epoch"];

a = [0 250];
b = [0 0];
c = [0 0];

p5 = plot3(a,b,c,"Color","black","LineWidth",2.5);
text(200,0,10,'Vernal Direction')

h = [p1;p2;p3;p4];

legend(h, l)