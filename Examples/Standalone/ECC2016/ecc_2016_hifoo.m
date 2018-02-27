function [K,CL,time] = ecc_2016_hifoo(G,WS,MS,WT,new)
%ECC_2016_LTI_MIXEDHINFSYNMIMO Summary of this function goes here
%   Detailed explanation goes here

options.output = struct('controller',0,'performance',1,'closedloop',0,'displaystyle','normal')

G = LTIsys(G);
WS = LTIsys(WS);
MS = LTIsys(MS);
WT = LTIsys(WT);

lti_begin(options)
    % Define the exogeneous input r
    signal r 

    % Define the connections
    u = G.in;       % Control input
    y = G.out;      % Control output
    e = r - y(1);   % tracking error
    
    % Define the control problems
    S = Channel(e/r, 'Sensitivity');
    U = Channel(u/r, 'Input sensitivity');
    T = Channel(y(1)/r, 'Complementary sensitivity');
    
    show(S,U,T)
    
%     ctrl_begin('classic')
%         K.in = e;       % Controller inputs: reference
%         K.out = u;      % Controller output
%         
%         minimize(WS*S)
%         subject to
%             MS*S <= 1
%             WT*T <= 1
%     ctrl_end
    
    options.FullOrderSolver = 'HIFOO';
    options.ReducedOrderSolver = 'HIFOO';
    ctrl_begin('new',options)
        if(new)
            K.in = [e;y];   % Controller inputs: reference, y(1), y(2)
        else
            K.in = [e];
        end
        K.out = u;      % Controller output
        
        minimize(WS*S)
        subject to
            MS*S <= 1
            WT*T <= 1
            if(new)
                order = 5;
            end
    ctrl_end
lti_end

K = extract('new');
CL = extract([e;y(1)]/r);
time = lti_problem.ctrl_prob(1).solver.info.time;
end

