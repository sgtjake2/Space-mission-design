clear, close all

moon = xlsread("Earth-Moon-2023 (2).csv");
jw = xlsread("JamesWebb-2023.csv");

figure(1)
hold on
grid on
m = plot3(moon(:,3),moon(:,4),moon(:,5));
j = plot3(jw(:,3),jw(:,4),jw(:,5));

view(3)
daspect([1 1 1])

