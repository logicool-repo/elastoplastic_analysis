% This program computes elastoplastic deformation of folded bare electro-dynamic tape tether

close all; clear all

global t E Ep fy P Qc I yc

b = 25E-3; % width = 25 mm
t = 45E-6; % thickness = 45 µm
I = b*t^3/12/b; % MOI of area
E = 5.1515E+9; % Young's modulus
Ep = 7.3838E+7; % Young's modulus
fy = 88.9E+6; % yield stress

y0 = [0 0]; % theta & dtheta at origin
s0 = 0; % length-wise distance at origin
sc = 25E-3; % length-wise distance at C = 25 mm
Qc = 0; % shearing stress at C is assumed to be zero
yc = pi/2+10*pi/180; % temporary

%% 100 g
P = 100E-3*9.80665/b; % pressure applied to the tape tether for crease formation

[S,Y] = ode45(@EPA,[s0 sc],y0);

imax = numel(Y(:,1));
x1 = zeros(imax,2);
dS = zeros(imax,1);
x0 = [0 0];

dS(1) = S(1);
x1(1,1) = x0(1)+dS(1)*cos(Y(1,1)); % x
x1(1,2) = x0(2)+dS(1)*sin(Y(1,1)); % z
x0(1) = x1(1,1); x0(2) = x1(1,2); % reset initial values

for i = 2:imax
    dS(i) = S(i)-S(i-1);
    x1(i,1) = x0(1)+dS(i)*cos(Y(i,1)); % x
    x1(i,2) = x0(2)+dS(i)*sin(Y(i,1)); % z
    x0(1) = x1(i,1); x0(2) = x1(i,2); % reset initial values
end
yc = Y(imax,1)*180/pi % update
lt = max(x1(:,1))*10^3 % layer thickness
ind_lt = find(x1(:,1)==max(x1(:,1)));

%% 200 g
P = 200E-3*9.80665/b; % pressure applied to the tape tether for crease formation

[S,Y] = ode45(@EPA,[s0 sc],y0);

imax = numel(Y(:,1));
x2 = zeros(imax,2);
dS = zeros(imax,1);
x0 = [0 0];

dS(1) = S(1);
x2(1,1) = x0(1)+dS(1)*cos(Y(1,1)); % x
x2(1,2) = x0(2)+dS(1)*sin(Y(1,1)); % z
x0(1) = x2(1,1); x0(2) = x2(1,2); % reset initial values

for i = 2:imax
    dS(i) = S(i)-S(i-1);
    x2(i,1) = x0(1)+dS(i)*cos(Y(i,1)); % x
    x2(i,2) = x0(2)+dS(i)*sin(Y(i,1)); % z
    x0(1) = x2(i,1); x0(2) = x2(i,2); % reset initial values
end
yc = Y(imax,1)*180/pi % update
lt = max(x2(:,1))*10^3 % layer thickness
ind_lt = find(x2(:,1)==max(x2(:,1)));

%% 300 g
P = 300E-3*9.80665/b; % pressure applied to the tape tether for crease formation

[S,Y] = ode45(@EPA,[s0 sc],y0);

imax = numel(Y(:,1));
x3 = zeros(imax,2);
dS = zeros(imax,1);
x0 = [0 0];

dS(1) = S(1);
x3(1,1) = x0(1)+dS(1)*cos(Y(1,1)); % x
x3(1,2) = x0(2)+dS(1)*sin(Y(1,1)); % z
x0(1) = x3(1,1); x0(2) = x3(1,2); % reset initial values

for i = 2:imax
    dS(i) = S(i)-S(i-1);
    x3(i,1) = x0(1)+dS(i)*cos(Y(i,1)); % x
    x3(i,2) = x0(2)+dS(i)*sin(Y(i,1)); % z
    x0(1) = x3(i,1); x0(2) = x3(i,2); % reset initial values
end
yc = Y(imax,1)*180/pi % update
lt = max(x3(:,1))*10^3 % layer thickness
ind_lt = find(x3(:,1)==max(x3(:,1)));

%% 400 g
P = 400E-3*9.80665/b; % pressure applied to the tape tether for crease formation

[S,Y] = ode45(@EPA,[s0 sc],y0);

imax = numel(Y(:,1));
x4 = zeros(imax,2);
dS = zeros(imax,1);
x0 = [0 0];

dS(1) = S(1);
x4(1,1) = x0(1)+dS(1)*cos(Y(1,1)); % x
x4(1,2) = x0(2)+dS(1)*sin(Y(1,1)); % z
x0(1) = x4(1,1); x0(2) = x4(1,2); % reset initial values

for i = 2:imax
    dS(i) = S(i)-S(i-1);
    x4(i,1) = x0(1)+dS(i)*cos(Y(i,1)); % x
    x4(i,2) = x0(2)+dS(i)*sin(Y(i,1)); % z
    x0(1) = x4(i,1); x0(2) = x4(i,2); % reset initial values
end
yc = Y(imax,1)*180/pi % update
lt = max(x4(:,1))*10^3 % layer thickness
ind_lt = find(x4(:,1)==max(x4(:,1)));

%% 500 g
P = 500E-3*9.80665/b; % pressure applied to the tape tether for crease formation

[S,Y] = ode45(@EPA,[s0 sc],y0);

imax = numel(Y(:,1));
x5 = zeros(imax,2);
dS = zeros(imax,1);
x0 = [0 0];

dS(1) = S(1);
x5(1,1) = x0(1)+dS(1)*cos(Y(1,1)); % x
x5(1,2) = x0(2)+dS(1)*sin(Y(1,1)); % z
x0(1) = x5(1,1); x0(2) = x5(1,2); % reset initial values

for i = 2:imax
    dS(i) = S(i)-S(i-1);
    x5(i,1) = x0(1)+dS(i)*cos(Y(i,1)); % x
    x5(i,2) = x0(2)+dS(i)*sin(Y(i,1)); % z
    x0(1) = x5(i,1); x0(2) = x5(i,2); % reset initial values
end
yc = Y(imax,1)*180/pi % update
lt = max(x5(:,1))*10^3 % layer thickness
ind_lt = find(x5(:,1)==max(x5(:,1)));

%% PLOT FIG
figure(1)
set(1,'Position',[40 40 900 600],'Color',[1 1 1]);
p = uipanel('Parent',figure(1),'BorderType','none'); 
p.Title = 'Elastoplastic deformation of folded tape tether';
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

subplot(1,5,1,'Parent',p);
h1 = plot(x1(:,1)*10^3,x1(:,2)*10^3); hold on
% h3 = plot(-x(:,1)*10^3,x(:,2)*10^3);
h2 = plot(x1(ind_lt,1)*10^3,x1(ind_lt,2)*10^3); hold off
set(h1,'LineStyle','-','LineWidth',1,'Color','b')
% set(h3,'LineStyle','-','LineWidth',1,'Color','b')
% set(h2,'LineStyle','o','LineWidth',1,'Color','b')
axis equal
axis([0 sc/5 0 sc]*10^3)
set(gca,'FontName','Arial','FontSize',18);
xlabel('x (mm)','FontName','Arial','FontSize',18);
ylabel('z (mm)','FontName','Arial','FontSize',18);
title('100gw')

subplot(1,5,2);
h1 = plot(x2(:,1)*10^3,x2(:,2)*10^3); hold on
% h3 = plot(-x(:,1)*10^3,x(:,2)*10^3);
h2 = plot(x2(ind_lt,1)*10^3,x2(ind_lt,2)*10^3); hold off
set(h1,'LineStyle','-','LineWidth',1,'Color','b')
% set(h3,'LineStyle','-','LineWidth',1,'Color','b')
% set(h2,'LineStyle','o','LineWidth',1,'Color','b')
axis equal
axis([0 sc/5 0 sc]*10^3)
set(gca,'FontName','Arial','FontSize',18);
xlabel('x (mm)','FontName','Arial','FontSize',18);
ylabel('z (mm)','FontName','Arial','FontSize',18);
title('200gw')

subplot(1,5,3);
h1 = plot(x3(:,1)*10^3,x3(:,2)*10^3); hold on
% h3 = plot(-x(:,1)*10^3,x(:,2)*10^3);
h2 = plot(x3(ind_lt,1)*10^3,x3(ind_lt,2)*10^3); hold off
set(h1,'LineStyle','-','LineWidth',1,'Color','b')
% set(h3,'LineStyle','-','LineWidth',1,'Color','b')
% set(h2,'LineStyle','o','LineWidth',1,'Color','b')
axis equal
axis([0 sc/5 0 sc]*10^3)
set(gca,'FontName','Arial','FontSize',18);
xlabel('x (mm)','FontName','Arial','FontSize',18);
ylabel('z (mm)','FontName','Arial','FontSize',18);
title('300gw')

subplot(1,5,4);
h1 = plot(x4(:,1)*10^3,x4(:,2)*10^3); hold on
% h3 = plot(-x(:,1)*10^3,x(:,2)*10^3);
h2 = plot(x4(ind_lt,1)*10^3,x4(ind_lt,2)*10^3); hold off
set(h1,'LineStyle','-','LineWidth',1,'Color','b')
% set(h3,'LineStyle','-','LineWidth',1,'Color','b')
% set(h2,'LineStyle','o','LineWidth',1,'Color','b')
axis equal
axis([0 sc/5 0 sc]*10^3)
set(gca,'FontName','Arial','FontSize',18);
xlabel('x (mm)','FontName','Arial','FontSize',18);
ylabel('z (mm)','FontName','Arial','FontSize',18);
title('400gw')

subplot(1,5,5);
h1 = plot(x5(:,1)*10^3,x5(:,2)*10^3); hold on
% h3 = plot(-x(:,1)*10^3,x(:,2)*10^3);
h2 = plot(x5(ind_lt,1)*10^3,x5(ind_lt,2)*10^3); hold off
set(h1,'LineStyle','-','LineWidth',1,'Color','b')
% set(h3,'LineStyle','-','LineWidth',1,'Color','b')
% set(h2,'LineStyle','o','LineWidth',1,'Color','b')
axis equal
axis([0 sc/5 0 sc]*10^3)
set(gca,'FontName','Arial','FontSize',18);
xlabel('x (mm)','FontName','Arial','FontSize',18);
ylabel('z (mm)','FontName','Arial','FontSize',18);
title('500gw')