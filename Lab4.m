clear, close all

r_min = 147.1e9;
e = [0.5 0.9 1.1 1.25 1.5 2 5];

theta = (-0.5*pi:0.1:(0.5*pi))';
theta2 = (1*pi:0.1:(1.4*pi))';

r = ((1+e)./(1+e.*cos(theta)))*r_min;
r2 = ((1+e)./(1+e.*cos(theta2)))*r_min;

x = (r/r_min).*cos(theta);
y = (r/r_min).*sin(theta);

x2 = (r2/r_min).*cos(theta2);
y2 = (r2/r_min).*sin(theta2);

rmax = 0.5e12;
for i = 1:length(e)
    figure(i)
    polarplot(theta,r(:,i))
    hold on
    polarplot(theta2,r2(:,i))
    rlim([0,rmax])
end

for i = 1:length(e)
    figure(i+10)
    plot(x(:,i),y(:,i))
    hold on
    plot(x2(:,i),y2(:,i))
    grid on
end
