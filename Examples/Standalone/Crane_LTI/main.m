%
% script for identification and control of the overhead crane with fixed
% cable length
%
clear all; close all; clc;

%% excitation signals
fs = 50;                                % sampling frequency
Ts = 1/fs;                              % sampling period
df = 0.01;                              
N  = 2^ceil(log10(1/(Ts*df))/log10(2)); % number of samples
df = fs/N;                              % frequency resolution

% generate multisines 
number_of_experiments = 5;
for k = 1:number_of_experiments
    musin{k} = Multisine('label',['musin' num2str(k)],'fs',fs,'fwindow',[df,5],'freqres',df);
end
musin{1}.plotSignal;
musin{1}.plotSpectrum('dft');

% IDEA: Intuitive syntax for generating multiple multisines with consistent labels at once???
% Store multiple multisines in one object...

%% experimental measurements
load measurements/meas1.mat
load measurements/meas2.mat
load measurements/meas3.mat
load measurements/meas4.mat
load measurements/meas5.mat

% store measurements in a cell array
meas{1} = meas1;
meas{2} = meas2;
meas{3} = meas3;
meas{4} = meas4;
meas{5} = meas5;

% unit conversions
L = 0.4557; % cable length
for k = 1:number_of_experiments
    t{k}         = meas{k}.X.Data;                      % time vector
    u{k}         = -meas{k}.Y(:,4).Data;                % applied input signal	
    pos_cart{k}  = meas{k}.Y(:,2).Data*(0.035*pi/9000); % horizontal cart position [m]
    angle_rad{k} = meas{k}.Y(:,1).Data*(2*pi/10000);    % cable angle [rad]
    y{k}         = pos_cart{k} + L*angle_rad{k};        % estimated load position [m]
end

% create MeasurementData objects
for k = 1:number_of_experiments
    experiment{k} = TDMeasurementData('label',['experiment' num2str(k)],'excitation',musin{k},'data',[u{k}(:),y{k}(:)],'datalabels',{'cart ref velocity','load position'},'periodic',true);
end
%experiment{1}.plotSignal('load position')
%experiment{1}.plotSignal('load position','periodic')

% de-trend measurements
for k = 1:number_of_experiments
    experiment{k} = experiment{k}.detrend({'load position'});
end
%experiment{1}.plotSignal('load position')
%experiment{1}.plotSignal('load position','periodic')

% IDEA: Automated 'clipping' of the time axis of the experimental data for 
% improved identification results

% REMARK: generated excitation signal and applied input signal are
% inconsistent, since random multisines have been used. However, this is 
% irrelevant for demonstration purposes.

%% identification nonparametric model
 
% define inputs and outputs
labels.input  = experiment{1}.datalabels_{1};
labels.output = experiment{1}.datalabels_{2};

% compute nonparametric FRF
[Gnonp,diag] = nonpar_ident(experiment{1},experiment{2},experiment{3},experiment{4},experiment{5},labels,'time2frf');
Gnonp = fromstd(Gnonp);
Gnonp.name = ('NonParametric');

bode(Gnonp); hold on;

% IDEA: wouldn't it be more intuitive to provide the MeasurementData of
% multiple experiments in a cell array (1 input argument), instead of
% varargin? I.e.:
% [nonpar_model_avg, diag] = nonpar_ident(experiments,labels,'time2frf');
% instead of
% [nonpar_model_avg, diag] = nonpar_ident(experiment{1},experiment{2},experiment{3},experiment{4},experiment{5},labels,'time2frf');

% TO DO: extension to MIMO models

%% fit a nominal parametric model
opt1 = struct('nh',3,'nl',0,'Bh',3,'Bl',0,'GN',0,'cORd','d','FRFW',1./sqrt(diag.sampleTotalVariance));
opt2 = struct('maxIter',500,'relVar', 1e-20,'GN',0,'FRFW',1,'gradTol',1e-8,'sysInit',[]); % default
opt = mergestruct(opt1,opt2);
Gparam = param_ident('data',Gnonp,'method','nllsfdi','options',opt); 
bode(Gparam);

%% derive uncertainty model

% REMARK: by studying the deviation of the FRFs corresponding to the different
% experiments, and assuming a multiplicative uncertainty model, one can
% derive a robustness bound to be used in the controller design step

% REMARK: We should think about how to generalize this. One can consider 
% various uncertainty models: parametric/unstructured, ... both parametric 
% and unstructured uncertainty can be taken into account in the solver 
% algorithms that we consider.

% Design robustness bound based on maximum relative error
for k = 1:number_of_experiments
    nonpar_model{k} = nonpar_ident(experiment{k},labels,'time2frf');
    FRF{k} = squeeze(nonpar_model{k}.ResponseData);
end
FRF_fit = squeeze(Gnonp.ResponseData);

relERROR = abs([(FRF{1} - FRF_fit)./FRF_fit,... 
                (FRF{2} - FRF_fit)./FRF_fit,...
                (FRF{3} - FRF_fit)./FRF_fit,... 
                (FRF{4} - FRF_fit)./FRF_fit,...   
                (FRF{5} - FRF_fit)./FRF_fit]);

maxrelERROR = max(relERROR,[],2);

% Parray = stack(1,nonpar_model{1},nonpar_model{2},nonpar_model{3},nonpar_model{4},nonpar_model{5});
% [u_model,Info] = ucover(Parray,nonpar_model_avg,'OutputMult');
% relerr = (nonpar_model_avg-Parray)/nonpar_model_avg;
% bodemag(relerr,'b--',Info.W1,'r',{0.01,50});

% [usys_new,info_new] = ucover(nonpar_model_avg,Info,2,0)
% relerr = (nonpar_model_avg-Parray)/nonpar_model_avg;
% bodemag(relerr,'b--',info_new.W1,'r',{0.01,50});

% IDEA: ucover does not seem to work well with this spiky error data... 
% maybe we should first smoothen the data before calling ucover.

% robustness bound
wd    = 5;      % crossover frequency [rad/s]
T0    = 10^(12/20); Tinf  = 10^(-40/20);    
numWT = 1.5*[1 3*wd/T0 3*(wd/T0)^2 (wd/T0)^3]; 
denWT = [Tinf^3 3*Tinf^2*wd 3*Tinf*wd^2 wd^3];
WT    = tf([numWT],[denWT]);
WT    = c2d(WT,Ts,'tustin');

% Plot max rel error and WT
figure; semilogx(Gnonp.Frequency,db(maxrelERROR),'black-'); hold on;
xlim([Gnonp.Frequency(1) Gnonp.Frequency(end)]);
wmin = 2*pi*Gnonp.Frequency(1); 
wmax = 2*pi*Gnonp.Frequency(end);
[mag,~,wout] = bode(WT,{wmin,wmax}); mag = squeeze(mag)'; magDB = db(mag);
plot(wout/2/pi,magDB,'black--');
box on; xlabel('frequency [Hz]'); ylabel('magnitude [dB]');
legend('maximum relative uncertainty','robustness bound','location','northwest');

%% controller design
% for now just reference tracking by bounding complementary sensitivity
% (using uncertainty model) and optimizing bandwidth

G = IOSystem(1,1);      % Make a system for the plant
Gparam.name  = 'parametric';
G.add(Gparam);          % Add the derived model to the plant
G.add(Gnonp);           % Add the non parametric model to the plant
K = IOSystem(1,1);      % Make a system for the controller

% Weight to maximize the bandwidth
WS = c2d(Weight.LF(0.5,1,-40),Ts,'Tustin'); % sensitivity weight

% Make the control configuration
r = Signal(); % reference input
u = G.in;     % control input
y = G.out;    % measured output
e = r - y;    % error on the absolute position

connections = [K.in == e;K.out == u];
CL = IOSystem(G,K,connections);

% Controller design
S = Channel(e/r,'Sensitivity');
T = Channel(y/r,'Compl. Sensitivity');

% Controller design 1
obj = [WS*S];
constr = [WT*T <= 1];
CL.solve(obj,constr,K);

% Controller design 2
obj2 = [WS*S];
constr2 = [WT*T <= 3];
CL.solve(obj2,constr2,K);

% Compute closed loop response
mod = CL.model();
S = mod(e,r);
U = mod(u,r);
T = mod(y,r);

% sensitivity and complementary sensitivity
bode(mod(y,r)); xlim([1e-1,5]); 
bode(mod(e,r)); xlim([1e-1,5]);

%step(mod(y,r));