%% LAB LINEAR GUIDES: DEMONSTRATION OF THE PROGRESS OF LCTOOLBOX
%  April 2017

%  Note: nor the identification experiments and control performance are 
%  too nice in this example, meant for illustrational purposes

clear; 
clc;
close all;

% 1. generate excitation signals
% ------------------------------

%     exc1 = Multisine('label', 'excitation1', ...
%                      'fwindow', [0.1 20], ...
%                      'fs', 1000, ...
%                      'freqres', 0.2);
%     exc2, exc3 and ex4 exactly the same

load('excitations.mat');    % we load the original signals: otherwise the
                            % phases would be different (since they're
                            % random)
                            
% 2. import, have a look at measurements and preprocess
% -----------------------------------------------------

    exp1 = importDSPACE('exp1.mat', exc1, true);
    exp2 = importDSPACE('exp2.mat', exc2, true);
    exp3 = importDSPACE('exp3.mat', exc3, true);
    exp4 = importDSPACE('exp4.mat', exc4, true);
    
    % for example: check the spectrum of the measured output of encoder 1
    plotSpectrum(exp1, {'enc1 (mm)'});
    
    % clip the measurements (same number of periods each time, but in a different ways)
    exp1c = clip(exp1, 'lastnper', 5);
    exp2c = clip(exp2, 'firstnper', 5);
    exp3c = clip(exp3, 'sample', [501 5500]);
    exp4c = clip(exp4, 'lastnper', 5);
    
% 3. non-parametric identification (motor mass)
% ---------------------------------------------
    
    IO.input = {'input (V)'};
    IO.output = {'enc1 (mm)'};
    Gnp = nonpar_ident(exp1c, exp2c, exp3c, exp4c, IO,'time2frf');
    
% 4. parametric identification
% ----------------------------

    FRF_W  = ones(size(squeeze(Gnp.ResponseData)));
    FRF_W(1:10) = 0;
    
    settings = struct('denh',4,'denl',0,'numh',2,'numl',0,'cORd','c','FRFW',FRF_W);
    Gp = param_ident('data', Gnp, 'method', 'nllsfdi', 'settings', settings); 
    
    figure;
    bode(Gnp); hold on; bode(Gp);
    
% 5. cast into an uncertain model
% -------------------------------

    % fit a multiplicative uncertainty weigth (bit artificial, only excited
    % up to 20 Hz here)
    
    nom = FRDmod(Gp,Gnp.Frequency,'Hz');
    
    resp1 = nonpar_ident(exp1c, IO, 'time2frf');
    resp2 = nonpar_ident(exp2c, IO, 'time2frf');
    resp3 = nonpar_ident(exp3c, IO, 'time2frf');
    resp4 = nonpar_ident(exp4c, IO, 'time2frf');
    resp = {resp1; resp2; resp3; resp4};
    relerror = cellfun(@(x) (x-nom)*inv(nom), resp, 'un', 0);
    
    s = tf('s');
    W1 = 5*((s+40)/(s+62))^8;
    W2 = tf(1,1);
    
    figure; hold on;
    bodemag(relerror{:});
    bodemag(W1);

% 6. design a controller
% ----------------------

    % the system's block with multiple models

        G = IOSystem(1,1);              % define the system's block

        par_model.name  = 'nominal';    % add the identified
        G.add(Gp);                      % (parametric) model

        resp1.name = 'measurement 1';   %
        G.add(resp1);                   %
        resp2.name = 'measurement 2';   %%
        G.add(resp2);                   %%%  add the non-parametric
        resp3.name = 'measurement 3';   %%%  models
        G.add(resp3);                   %% 
        resp4.name = 'measurement 4';   %
        G.add(resp4);                   %

    % the controller block
    
        K = IOSystem(1,1);
        
    % the control configuration
    
        r = Signal();
        u = G.in;
        y = G.out;
        e = r - y;
    
        connections = [K.in == e; K.out == u];
        CL = IOSystem(G,K,connections);

    % define the channels and the specificatins
    
        S = Channel(e/r,'Sensitivity');
        T = Channel(y/r,'Complementary sensitivity');
        
        MS = Weight.DC(2);
        WS = Weight.LF(20,4,-100);
        WT = W1;
    
        obj = [WS*S];
        constr = [WT*T <= 1, MS*S <= 1];
        
    % solve

        CL.solve(obj,constr,K);
        mod = CL.model();
        
    % check the results

        figure; bode(K); title('Controller');
        figure; bodemag(CL(S)); title('Sensitivity');
        figure; bodemag(CL(T)); title('Complementary sensitivity');
        figure; step(CL(T)); title('Step response');