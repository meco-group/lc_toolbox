%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for the ECC2016 conference controller synthesis using the LTI
% control toolbox
% Matlab's robust control (hinfstruct and systune) toolbox is required to
% run the examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

ploton = true;

%% System definition - load all systems as LTIsys so they are automatically loaded by the toolbox
load('MSM_model.mat','model')
G = model;

% Define the optimal gamma which has been retrieved via mixedHinfsynMIMO
% This way we get the same, optimal controller for mixedHinfsynMIMO and
% systune/hinfstruct
gamma = [0.9068;...
        0.8277;...
        0.9954];

% Optimal bandwidth weight both stable and unstable
WSs = Weight.LF(5,3,-90)/gamma(1); %-80dB flatgain, cross over frequency at 5 Hz, 3th order
WSu = Weight.LFu(5,3)/gamma(1);    %-infdB flatgain, cross over frequency at 5 Hz, 3th order  
% Maximum sensitivity weight
MS = Weight.DC(8);        %4dB maximum sensitivity
% Robustness weight
WTp = Weight.HF(13,2,-40);
WTi = Weight.HFu(13,2);

%% Compute controllers
% All controller specifications are listed in the following structure
list(1) = struct('WS',WSu,'MS',MS,'WT',WTp,'new',false,'func','mixedhinfsynmimo','K',[],'CL',[]);
list(2) = struct('WS',WSs,'MS',MS,'WT',WTp,'new',false,'func','hinfstruct','K',[],'CL',[]);
list(3) = struct('WS',WSs,'MS',MS,'WT',WTp,'new',false,'func','systune','K',[],'CL',[]);
list(4) = struct('WS',WSu,'MS',MS,'WT',WTi,'new',false,'func','mixedhinfsynmimo','K',[],'CL',[]);
list(5) = struct('WS',WSs,'MS',MS,'WT',WTi,'new',false,'func','hinfstruct','K',[],'CL',[]);
list(6) = struct('WS',WSs,'MS',MS,'WT',WTi,'new',false,'func','systune','K',[],'CL',[]);
list(7) = struct('WS',WSu,'MS',MS,'WT',WTp,'new',true,'func','mixedhinfsynmimo','K',[],'CL',[]);
list(8) = struct('WS',WSs,'MS',MS,'WT',WTp,'new',true,'func','hinfstruct','K',[],'CL',[]);
list(9) = struct('WS',WSs,'MS',MS,'WT',WTp,'new',true,'func','systune','K',[],'CL',[]);
list(10) = struct('WS',WSu,'MS',MS,'WT',WTi,'new',true,'func','mixedhinfsynmimo','K',[],'CL',[]);
list(11) = struct('WS',WSs,'MS',MS,'WT',WTi,'new',true,'func','hinfstruct','K',[],'CL',[]);
list(12) = struct('WS',WSs,'MS',MS,'WT',WTi,'new',true,'func','systune','K',[],'CL',[]);

% Solve for all controllers in the list
for k=[4:6,10:12]
    solver = str2func(['ecc_2016_' list(k).func]);
    tic;
    [list(k).K,list(k).CL,list(k).solver_time] = solver(G,list(k).WS,list(k).MS,list(k).WT,list(k).new);
    list(k).total_time = toc;
end

%% show results
save('ecc_2016_list_new','list');

if(ploton)
    %% Plot G
    freq = logspace(0,2,100);
    Gresp = transpose(squeeze(freqresp(G,freq,'Hz')));

    figure('Position',[100, 100, 600 380]), 
    subplot(411)
    semilogx(freq,db(Gresp(:,1)),'LineWidth',2)
    % xlabel('f [Hz]')
    ylabel('|Y_1/F| [dB]')
    subplot(412)
    semilogx(freq,(180/pi)*angle(Gresp(:,1)),'LineWidth',2)
    % xlabel('f [Hz]')
    ylabel('<Y_1/F [\circ]')
    subplot(413)
    semilogx(freq,db(Gresp(:,2)),'LineWidth',2)
    % xlabel('f [Hz]')
    ylabel('|Y_2/F| [dB]')
    subplot(414)
    semilogx(freq,(180/pi)*angle(Gresp(:,2)),'LineWidth',2)
    xlabel('f [Hz]')
    ylabel('<Y_2/F [\circ]')

    %% Plot weights
    freq = logspace(-2,3,200);
    MSresp = db(squeeze(freqresp(MS,freq,'Hz')));
    WSsresp = db(squeeze(freqresp(WSs,freq,'Hz')));
    WSuresp = db(squeeze(freqresp(WSu,freq,'Hz')));
    WTpresp = db(squeeze(freqresp(WTp,freq,'Hz')));
    WTiresp = db(squeeze(freqresp(WTi,freq,'Hz')));

    figure('Position',[100, 100, 600 380]),
    semilogx(freq,-WSsresp,'b--',...
             freq,-WSuresp,'b',...
             freq,-WTpresp,'r--',...
             freq,-WTiresp,'r',...
             freq,-MSresp,'k','LineWidth',2)
    legend('1/WS_{stable}','1/WS_{unstable}','1/WT_{proper}','1/WT_{improper}','1/MS');
    xlabel('f [Hz]')
    ylabel('magnitude [dB]')
    ylim([-120,30])
    
    %% Show controller performances
    show_list(list);    
end