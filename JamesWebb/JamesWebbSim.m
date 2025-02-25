clear
close all

%JPL Horizons Data 2022/06/01

posS = [-1.342265108825705E+06 3.044515635495678E+05 2.884704529257542E+04]*1000;
posE = [-5.270751030639560E+07 -1.424211399657750E+08 3.602153098984063E+04]*1000;
posJ = [-5.256168304341657E+07 -1.440483280181565E+08 -1.734837072574794E+05]*1000;

velS = [-3.050634513857201E-03 -1.557381710015418E-02 1.925613372905045E-04]*1000;
velE = [2.755293933513340E+01 -1.021705234088606E+01 -3.346080722117506E-04]*1000;
velJ = [2.771985965841215E+01 -1.024983981129689E+01 -1.464231284627822E-01]*1000;

G = 6.674e-11; %kg^-1m^3s^-2
mE = 5.97219e24;
mS = 1988410e24;

out = sim("jamesWebb");

OS = out.OS;
OE = out.OE;
OJ = out.OJ;

figure(1)
grid on
hold on
e = plot3(OE.Data(:,1),OE.Data(:,2),OE.Data(:,3));
s = plot3(OS.Data(:,1),OS.Data(:,2),OS.Data(:,3),"Marker","*");
j = plot3(OJ.Data(:,1),OJ.Data(:,2),OJ.Data(:,3));

view(3)

hJ = animatedline("Color","red","LineStyle","--","Marker","o");
h = animatedline("Color","red","LineStyle","--");
hE = animatedline("Color","green","LineStyle","--","Marker","o");

%daspect([1 1 1])

%{
for i = 1:2000:length(OE.Data)
    addpoints(hE,OE.Data(i,1),OE.Data(i,2),OE.Data(i,3))
    addpoints(h,OJ.Data(i,1),OJ.Data(i,2),OJ.Data(i,3))
    addpoints(hJ,OJ.Data(i,1),OJ.Data(i,2),OJ.Data(i,3))

    pause(0.0001)

    clearpoints(hE)
    clearpoints(hJ)
end
%}

dataE = xlsread("Earth-Sun-2022-2024.csv");
dataJ = xlsread("JamesWebb-Sun-2022-2024.csv");

figure(2)
hold on
grid on

datE = plot3(dataE(:,3),dataE(:,4),dataE(:,5));
datJ = plot3(dataJ(:,3),dataJ(:,4),dataJ(:,5));


view(3)
%daspect([1 1 1])

hJ1 = animatedline("Color","red","LineStyle","--","Marker","o");
h1 = animatedline("Color","red","LineStyle","--");
hE1 = animatedline("Color","green","LineStyle","--","Marker","o");

%{
for b = 1:length(dataE)
    addpoints(hE1,dataE(b,3),dataE(b,4),dataE(b,5))
    addpoints(h1,dataJ(b,3),dataJ(b,4),dataJ(b,5))
    addpoints(hJ1,dataJ(b,3),dataJ(b,4),dataJ(b,5))

    pause(0.01)

    clearpoints(hE1)
    clearpoints(hJ1)
end
%}

figure(3)
hold on
grid on
j = plot3(OJ.data(:,1)-OE.data(:,1), OJ.data(:,2)-OE.data(:,2),OJ.data(:,3)-OE.data(:,3));
view(3)
daspect([1 1 1])

figure(4)
hold on
grid on
j = plot3(dataJ(:,3)-dataE(:,3), dataJ(:,4)-dataE(:,4),dataJ(:,5)-dataE(:,5));
view(3)
daspect([1 1 1])