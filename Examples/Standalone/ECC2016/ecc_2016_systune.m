function [K,CL,time] = ecc_2016_systune(G,WS,MS,WT,new)
%ECC_2016_RC_SYSTUNE2 Summary of this function goes here
%   Do Hinf controller design with systune

%% Construct controller block
Ny = 1;
if(new)
    Nx = 5;    
    Nu = 3;
    K = ltiblock.ss('K',Nx,Ny,Nu);
    K.u = {'e';'y';'ya'}; K.y = 'u';
else
    Nx = 7 + length(pole(WT));    
    Nu = 1;
    K = ltiblock.ss('K',Nx,Ny,Nu);
    K.u = 'e'; K.y = 'u';
end

%% Make the interconnection
WS.u = 'e';  WS.y = 'z1';
MS.u = 'e';  MS.y = 'z2';
WT.u = 'y'; WT.y = 'z3';
G.u = 'u'; G.y = {'y';'ya'};

% Specify summing junctions
Sum1 = sumblk('e = r - y');

% Connect the blocks together
P = connect(G,K,WS,MS,WT,Sum1,{'r'},{'z1','z2','z3'})

%% synthetize controller
ReqWS = TuningGoal.Gain('r','z1',1);
ReqMS = TuningGoal.Gain('r','z2',1);
ReqWT = TuningGoal.Gain('r','z3',1);

systuneOptions('RandomStart',5);
tic;
[CL,fSoft,fHard,info] = systune(P,[ReqWS],[ReqMS;ReqWT]);
time = toc;
K = getBlockValue(CL,'K');
if(new)
    K.u = {'e';'y';'ya'}; K.y = 'u';
else
    K.u = 'e'; K.y = 'u';
end

CL = connect(G,K,Sum1,{'r'},{'e','y'});
end

