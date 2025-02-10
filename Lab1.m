%% Ecliptic

clear, close all
data = xlsread("Sun-Earth-2023-2025(1).csv");
data_sun = xlsread("Sun-Sun-2023-2024.csv");
data_mars = xlsread("Sun-Mars-2023-2025-ecliptic.csv");

x = data(:,3);
y = data(:,4);
z = data(:,5);

x1 = data_sun(:,3);
y1 = data_sun(:,4);
z1 = data_sun(:,5);

x2 = data_mars(:,3);
y2 = data_mars(:,4);
z2 = data_mars(:,5);

figure(1)
plot3(x,y,z)
hold on
grid on
plot3(x1,y1,z1,'o','Color','yellow','MarkerSize',10,...
    'MarkerFaceColor','yellow')
plot3(x2,y2,z2)
xlabel("X")
ylabel("Y")
zlabel("Z")
title("Ecliptic")
daspect([1 1 1])

t = "2024-12-01";
t = datetime(t);

jd = juliandate(t);

for i = 1:length(data)
    if data(i,1) == jd
        plot3(data(i,3),data(i,4),data(i,5),'o','Color','blue')
        %plot3(data_sun(i,3),data_sun(i,4),data_sun(i,5),'*','Color','yellow')
        plot3(data_mars(i,3),data_mars(i,4),data_mars(i,5),'o','Color','red')
    end
end

figure(3)
plot(x,y)
hold on
grid on
plot(x1,y1,'o','Color','yellow','MarkerSize',10,...
    'MarkerFaceColor','yellow')
plot(x2,y2)
xlabel("X")
ylabel("Y")
%daspect([1 1])


%% Equatorial

data = xlsread("Sun-Earth-2023-2024-equotoarl(1).csv");
data_sun = xlsread("Sun-Sun-2023-2025-equatorial.csv");
data_mars = xlsread("Sun-Mars-2023-2025-ecquatorial.csv");

x = data(:,3);
y = data(:,4);
z = data(:,5);

x1 = data_sun(:,3);
y1 = data_sun(:,4);
z1 = data_sun(:,5);

x2 = data_mars(:,3);
y2 = data_mars(:,4);
z2 = data_mars(:,5);

figure(2)
plot3(x,y,z)
hold on
grid on
plot3(x1,y1,z1,'o','Color','yellow','MarkerSize',10,...
    'MarkerFaceColor','yellow')
plot3(x2,y2,z2)
xlabel("X")
ylabel("Y")
zlabel("Z")
title("Equatorial")
daspect([1 1 1])

t = "2024-12-01";
t = datetime(t);

jd = juliandate(t);

for i = 1:length(data)
    if data(i,1) == jd
        plot3(data(i,3),data(i,4),data(i,5),'o','Color','blue')
        %plot3(data_sun(i,3),data_sun(i,4),data_sun(i,5),'*','Color','yellow')
        plot3(data_mars(i,3),data_mars(i,4),data_mars(i,5),'o','Color','red')
    end
end

%% Offset angle

data = xlsread("Sun-Earth-2023-2025(1).csv");
data1 = xlsread("Sun-Earth-2023-2024-equotoarl(1).csv");

pos = [data(1,6);data(1,7);data(1,8)];
pos1 = [data1(1,6);data1(1,7);data1(1,8)];

direction = norm(pos);
direction1 = norm(pos1);

angle = acos(direction/direction1);