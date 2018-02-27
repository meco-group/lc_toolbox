%% Example : XY Table based on the previous measurements 

clear all
close all
clc

%% Creating the excitation signal
N = 5;
fmin = 0.1;
fmax = 50;
for i = 1:N
InputCurrent{i} = Multisine('label','Current_in',...
    'fwindow',[0.1,50],'fs',1000, 'freqres', 0.12207);
end

%% Load already measured data
load longMeasurement_01508_1.mat;
myInput{1} = myInput1;
Data_measured{1} = [mea_cur1 mea_pos1];
load longMeasurement_01508_3.mat;
myInput{2} = myInput1;
Data_measured{2} = [mea_cur1 mea_pos1];
load longMeasurement_01508_5.mat;
myInput{3} = myInput1;
Data_measured{3} = [mea_cur1 mea_pos1];
load longMeasurement_01508_7.mat;
myInput{4} = myInput1;
Data_measured{4} = [mea_cur1 mea_pos1];
load longMeasurement_01508_9.mat;
myInput{5} = myInput1;
Data_measured{5} = [mea_cur1 mea_pos1];
clearvars -except myInput InputCurrent Data_measured N fmin fmax
%% Changing the reference Inputs
% This is done since we already have the measured data from previous
% experiments. Its nothing but overwriting the multisine inputs, to match
% the reference input used, because well it was random multisine.
for j= 1:N
    InputCurrent{j}.signal_.Time = myInput{j}(1,1:8192);
    InputCurrent{j}.signal_.Value = myInput{j}(2,1:8192);
end
InputCurrent{1}.plotSpectrum('dft');
%% Creating the MeasurementData class
for k=1:N
Experiment{k} = MeasurementData('label','ExperimentXY_io',...
    'excitation', InputCurrent{k}, 'data' , Data_measured{k} ,...
    'datalabels',{'current','position'}, 'periodic', true);
end
%% Preparing for the nonParametric Identification

label.input = {'current'};
label.output = {'position'};
[model,diag] = nonpar_ident(Experiment{1},Experiment{2},Experiment{3},Experiment{4},Experiment{5},label,'time2frf');

G = fromstd(model);
G.name = ('NonParametric');
figure;
bode(G);
% chosen between the excited frequency window 

%% fit a nominal parametric model
FRF_W  = ones(size(squeeze(G.ResponseData)));
options = IdentOpt('n',2,'Bh',0,'Bl',0,'GN',0,'cORd','c','FRFW',FRF_W,'maxIter',200);
measdata = idfrd(squeeze(G.ResponseData),2*pi*G.Frequency,0.001);
Gparam = param_ident(measdata, 'nllsfdi', options);  

bode(G); % FRF parametric model
hold on;
bode(Gparam)
%% Creating Interconnections for Systems
H = IOSystem(1,1);
%H = H.add(G);
H = H.add(Gparam);
K = IOSystem(1,1);

%% Weights description
MS = Weight.DC(5); % weightDC constructs a dc weight equivalent with a peak of 5db
WS = Weight.LF(400/(2*pi),1,-100); % Construct a first order low frequency roll-off weight with a cross-over frequency of 50Hz
WU = Weight.HF(600/(2*pi),1,-30); % Weight on the input sensitivity to enforce roll-off in the controller

%% control configuration
r = Signal();u = H.in;y = H.out;e = r - y;

connections = [u==K.out;e==K.in];
I = IOSystem(H,K,connections);

CL = IOSystem(H,K,connections);

%% Controller design
S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Actuator Effort');

%%
CL.solve(WS*S,WU*U <= 1,K);


%% Compute closed loop response
mod = CL.model();
S = mod(e,r);
U = mod(u,r);
T = mod(y,r);

bode(mod([e;y],r));