clear global
clear all
clc
close all

addpath(genpath('../lti_toolbox'));

%% 1. Connecting an augmented plant 
% Simple continuous time controller design for the X-stage of the XYtable
% at PMA. Special thanks to Dora for supplying the model and testing the
% toolbox.

%% 1.1. Let's define the systems we will be using
load identifiedX.mat

SYSc = d2c(SYS,'matched'); %Choose either matched or foh discretization in order to avoid nmp zero
G = fromstd(SYSc);

% [m1,p1] = bode(SYS,2*pi*frequency);
% [mc,pc] = bode(SYSc,2*pi*frequency);
% 
% figure,subplot(211),semilogx(frequency,db(abs(squeeze([m1, mc]))))
% subplot(212),semilogx(frequency,(squeeze([p1 (pc-360)])))

MS = Weight.DC(5); % weightDC constructs a dc weight equivalent with a peak of 5db
WS = Weight.LF(400/(2*pi),1,-100); % Construct a first order low frequency roll-off weight with a cross-over frequency of 50Hz
WU = Weight.HF(600/(2*pi),1,-30); % Weight on the input sensitivity to enforce roll-off in the controller

% Make a second WS to have additional lag
Ts = (2*pi)/300;    % Lag time
alpha = 100;        % amount of lag
lag = tf(alpha*[Ts 1],[alpha*Ts 1]);  % Lag filter
WS2 = WS*lag*0.9; % Final WS2: WS+lag+(slight gain of 0.9 to cope with conservatism)

%% 1.2. Design the optimal controller using the lti_toolbox
% Let's also set some options to get the plots we want
options.output = struct('controller',1,'performance',1,'closedloop',0,'displaystyle','normal');

lti_begin(options)
    % Define variables to do the eventual controller design
    r = Signal();
    u = G.in;
    y = G.out;
    e = r - y;

    connect to
    K.in = e;               % Assign the error to the controller input
    K.out = u;              % Assign the plant input to the controller output

    % Do the controller design
    S = Channel(e/r,'Sensitivity');
    U = Channel(u/r,'Input Sensitivity');
    T = Channel(y/r,'Complementary Sensitivity');
    show(T);
    
    ctrl_begin('minimal input sensitivity')
        minimize(WU*U)
            MS*S <= 1;
            WS*S <= 1;
    ctrl_end
    
    ctrl_begin('additional lag')
        minimize(WU*U)
            MS*S <= 1;
            WS2*S <= 1;
    ctrl_end
lti_end

%% Processing & plotting
C = extract({'minimal input sensitivity';'additional lag'}); % CT controller
C = cellfun(@(x)shiftlow(std(x),1e-2),C,'UniformOutput',false);
C = cellfun(@(x)removehigh(x,1e4),C,'UniformOutput',false);

Cd = cellfun(@(x)c2d(x,1e-3,'tustin'),C,'UniformOutput',false); % DT controller
Gol = cellfun(@(x)series(x,SYS),Cd,'UniformOutput',false);
Gcl = cellfun(@(x)feedback(x,1),Gol,'UniformOutput',false);
figure(),
subplot(211)
bode(Cd{:});
subplot(212)
bode(Gol{:});

figure(), title('Closed loop')
bode(Gcl{:});
figure(), title('Compare')
bode(C{:})
figure, title('Open loop')
bode(C{1}*std(G),C{2}*std(G))