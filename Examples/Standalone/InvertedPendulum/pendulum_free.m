clear all
close all
clc

%% Parameters
M = 0.05;
m = 0.010;
l = 0.025;
g = 9.81;
J = 0;
c = 1e-6;

timespan = [0,20];

%% Stable position - down
y0s = [pi,0,0,0];
u = @(t) 1*sin(2*pi*t);
[ts,ys] = ode45(@(t,x) pendulum_model(t,x,u,m,l,g,J,c), timespan, y0s);
us = u(ts);
dys = zeros(size(ys));
for k = 1:length(ts)
    dys = pendulum_model(ts(k,1),ys(k,:),u,m,l,g,J,c);
end

% Linearized model
Gs = ss([0 1;-m*g*l/(J+m*l^2) -c/(J+m*l^2)],[0;-m*l/(J+m*l^2)],[1 0],[0]);
ts_lin = linspace(timespan(1), timespan(2), 2000);
us_lin = u(ts_lin);
[ys_lin] = lsim(Gs,us_lin,ts_lin);
ys_lin = ys_lin+pi;

% Plotting
figure()
subplot(211)
plot(ts,rad2deg(ys(:,1)),'b-',ts_lin,rad2deg(ys_lin),'r-');
legend('Nonlinear','Linear')
subplot(212)
plot(ts,rad2deg(ys(:,3)),'k-');

%% Unstable position - up
y0u = [0.001,0,0,0];
u = @(t) 0;
[tu,yu] = ode45(@(t,x) pendulum_model(t,x,u,m,l,g,J,c), timespan, y0u);
dyu = zeros(size(yu));
for k = 1:length(tu)
    dyu = pendulum_model(tu(k,1),yu(k,:),u,m,l,g,J,c);
end

% Linearized model
Gu = ss([0 1;m*g*l/(J+m*l^2) -c/(J+m*l^2)],[0;m*l/(J+m*l^2)],[1 0],[0]);
tu_lin = linspace(timespan(1), timespan(2), 2000);
uu_lin = tu_lin*0;
[yu_lin] = lsim(Gu,uu_lin,tu_lin,[0.001;0]);
yu_lin = yu_lin+pi;

figure()
subplot(211)
plot(tu,rad2deg(yu(:,1)),'b-',tu_lin,rad2deg(yu_lin),'r-');
legend('Nonlinear','Linear')
subplot(212)
plot(tu,yu(:,3),'k-');