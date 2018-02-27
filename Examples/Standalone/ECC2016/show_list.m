function show_list(list)
%SHOW_LIST Displays the controller performance
%   Displays the performance of different controller designs
%   Figure 1: Controller for 3 single input controllers with improper WT
%   Figure 2: Closed-loop performance for single input controllers with
%   improper WT
%   Figure 3: Closed-loop performance for multiple input controllers with
%   improper WT

    if(isstr(list))
        load(list);
    end

    freq = logspace(-2,4,200);
    K1 = [list(4:6).K];
    CL1S = [list(4).CL(1,1),list(5).CL(1,1),list(6).CL(1,1)];
    CL1T = [list(4).CL(2,1),list(5).CL(2,1),list(6).CL(2,1)];
    CL2S = [list(10).CL(1,1),list(11).CL(1,1),list(12).CL(1,1)];
    CL2T = [list(10).CL(2,1),list(11).CL(2,1),list(12).CL(2,1)];
    K1resp = transpose(squeeze(freqresp(K1,freq,'Hz')));
    CL1Sresp = transpose(squeeze(freqresp(CL1S,freq,'Hz')));
    CL1Tresp = transpose(squeeze(freqresp(CL1T,freq,'Hz')));
    CL2Sresp = transpose(squeeze(freqresp(CL2S,freq,'Hz')));
    CL2Tresp = transpose(squeeze(freqresp(CL2T,freq,'Hz')));
    
    %% First controller response
    h1 = figure('Name','Controller 1 - bode');
    set(h1, 'Position', [100, 100, 600 380]);
    subplot(211),
    semilogx(freq,db(K1resp),'LineWidth',2);
    xlabel('f [Hz]');
    ylabel('magnitude [dB]');
    legend(list(4).func,list(5).func,list(6).func);
    subplot(212),
    K1resp_angle = (180/pi)*unwrap(angle(K1resp));
    K1resp_angle(:,1) = K1resp_angle(:,1)-360;
    semilogx(freq,K1resp_angle,'LineWidth',2);
    xlabel('f [Hz]');
    ylabel('phase [\circ]');
    
    %% First controller response
    h2 = figure('Name','Controller 1 - closed loop');
    set(h2, 'Position', [100, 100, 600 380]);
    subplot(121),
    semilogx(freq,db(CL1Sresp),'LineWidth',2);
    xlabel('f [Hz]');
    ylabel('Sensitivity [dB]');
    legend(list(4).func,list(5).func,list(6).func,'Location','southeast');
    subplot(122),
    semilogx(freq,db(CL1Tresp),'LineWidth',2);
    xlabel('f [Hz]');
    ylabel('Complementary Sensitivity [dB]');
    
    %% First controller response
    h3 = figure('Name','Controller 2 - closed loop');
    set(h3, 'Position', [100, 100, 600 380]);
    subplot(121),
    semilogx(freq,db(CL2Sresp),'LineWidth',2);
    xlabel('f [Hz]');
    ylabel('Sensitivity [dB]');
    legend(list(10).func,list(11).func,list(12).func,'Location','southeast');
    subplot(122),
    semilogx(freq,db(CL2Tresp),'LineWidth',2);
    xlabel('f [Hz]');
    ylabel('Complementary Sensitivity [dB]');

end

