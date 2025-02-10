clear, clc, close all
% aero419_orbital_elements_ex01.m
% Aim: to interpret orbital elements as per example in notes.
% More info to follow:
% DJ Walker, 02/11/2022
% Earth Moon Barycentre
a=1.00000011; % semi-major axis (AU)
e=0.01671022; % orbital eccentricity
i=0.00005; % orbital inclination (deg)
%W=-11.26064; % longitude of ascending node (deg)
W=348.73936; % longitude of ascending node (RAAN)(deg)
% Next two are not Keplerian Elements
wbar=102.94719; % longitude of perihelion (deg) 102.94719 (*)
Lambda=100.46435; % mean longitude (deg) 100.46435
w=wbar-W; % argument of perihelion (deg) (A Keplerian element)
if w<0
w=w+360;
end
% calculations:
M=Lambda-wbar; % mean anomaly (deg)
if M<0
M=M+360;
end
% Now solver Kepler's equation for Eccentric Anomaly
E=mean_to_eccentric_anomaly(M*pi/180,e); % (E is in rad)
cosE=cos(E);
sinE=sin(E);
costheta=(cosE-e)/(1-e*cosE);
sintheta=sqrt(1-e^2)*sinE/(1-e*cosE);
theta=atan2(sintheta,costheta);
if theta<0
theta=theta+2*pi;
end
theta=theta*180/pi; % convert true anomaly back to deg
fprintf('M=%8.6frad =%8.3fdeg\n',M*pi/180,M)
fprintf('E=%8.6frad =%8.3fdeg\n',E,E*180/pi)
fprintf('cos(theta0)=%8.7f\n',costheta)
fprintf('sin(theta0)=%8.7f\n',sintheta)
fprintf('theta=%8.6fdeg =%8.3frad\n',theta,theta*pi/180)
fprintf('Argument of perhelion (w)=%8.6fdeg\n',w)
% Generate and plot orbit in inertial coordinates
Npoints=1000;
[ri,opi]=aero419_orbital_elements_to_inertial_xyz(a,e,i,W,w,theta,Npoints);
plot3(ri(1,:),ri(2,:),ri(3,:),'k-')
hold on
plot3([0,opi(1)],[0,opi(2)],[0,opi(3)],'k--') % Draw line from O to Periapsis
daspect([1,1,1])
grid
xlabel('x')
ylabel('y')
zlabel('z')
