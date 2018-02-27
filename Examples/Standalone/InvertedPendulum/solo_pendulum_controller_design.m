clear all
clear global
close all
clc

% Parameters
M = 0.05;
m = 0.010;
l = 0.025;
g = 9.81;
J = 0;
c = 1e-6;

% Model
Gm = DSSmod([0 0 1 0;0 0 0 1;0 0 0 0;0 m*g*l 0 -c],[0 0;0 0;1 0;0 1],...
           [1 0 0 0;0 1 0 0],[0 0;0 0],...
           [1 0 0 0;0 1 0 0;0 0 M+m -m*l;0 0 -m*l J+m*l*l]);
G = IOSystem(Gm);
K = IOSystem(3,1);

% Weights
WS = Weight.LF(0.5,2,-30);
MS = Weight.DC(8);
WT = Weight.HF(10,2);
MS_theta = Weight.DC(20);

%% Controller design
% define systems
r = Signal()

% Connections
F = G.in(1);
T = G.in(2);
x = G.out(1);
theta = G.out(2);
e = Signal();

% Controller in/outputs
c1 = (K.in == [r;x;theta]);
c2 = (K.out == F);
c3 = (e == r-x);
P = IOSystem(G,K,[c1;c2;c3]);
    
% Control problem formulation
Sx = Channel(e/r,'Sensitivity X');
St = Channel(theta/r,'Sensitivity Theta');
CS = Channel(x/r,'Complementary Sensitivity');
    
obj1 = WS*Sx;
constr1 = [MS*Sx <= 1, WT*CS <= 1];
opts.gamma.solver = 'lmilab';
opts.controller_name = 'Controller1';
[P,C1,info1] = P.solve(obj1,constr1,K,opts);

obj2 = WS*Sx;
constr2 = [MS*Sx <= 1, WT*CS <= 1, MS_theta*St <= 1];
opts.controller_name = 'Controller2';
[P,C2,info2] = P.solve(obj2,constr2,K,opts);

figure, bodemag(info1,info2);
figure, bode(K);

%% Simulation part
Ts = 0.005;
t = 0:Ts:20;

u = @(t)(((t>5)&&(t<15))*[0.01;0] + ((t>10)&&(t<20))*[0;m*g*l*sin(1*pi/180)]);
[y,t,x] = sim(P([x;theta]/[r;T]),u,t);

figure()
subplot(211)
plot(t{1},rad2deg(y{1}(:,2)),t{2},rad2deg(y{2}(:,2)));
xlabel('t [s]')
ylabel('theta [\circ]')
legend('Controller1','Controller2')
subplot(212)
plot(t{1},y{1}(:,1),t{2},y{2}(:,1));
xlabel('t [s]')
ylabel('X [m]')
legend('Controller1','Controller2')

%% Save simulation to animation
show_gif = false;
if show_gif
    param = struct('l',M/4,'wc',M/2,'hc',M/3.5,'wp',l/4,'hp',l*2);
    f = figure()
    ax = axes('Parent',f);
    filename = 'solo_pendulum.gif';
    
    dt = 0.05;
    k = 1;
    for tau = 4:dt:20
        cla(ax,'reset');
        draw_pendulum(ax,interp1(t{2},y{2}(:,1),tau),interp1(t{2},2*y{2}(:,2),tau),param,[-0.025 0.045 0 0.07]);
        drawnow;
        frame = getframe(f);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256); 
            if k == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
            end
%         save2pdf(sprintf('inverted_pendulum_solo/frame_%d',k),f,300); k = k+1;
    end
end
