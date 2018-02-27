function [K,CL,time] = ecc_2016_hinfstruct(G,WS,MS,WT,new)
%ECC_2016_RC_HINFSTRUCT synthetize controller with hinfstruct
%   Detailed explanation goes here

%% Make controller block 
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
opt = hinfstructOptions('Display','final','RandomStart',5);
tic;
[Ps,gam,info] = hinfstruct(P,opt);
time = toc;
K = getBlockValue(Ps,'K');
if(new)
    K.u = {'e';'y';'ya'}; K.y = 'u';
else
    K.u = 'e'; K.y = 'u';
end

CL = connect(G,K,Sum1,{'r'},{'e','y'});
end

