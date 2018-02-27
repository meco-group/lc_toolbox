% This file is part of LCToolbox.
% (c) Copyright 2018 - MECO Research Team, KU Leuven. 
%
% LCToolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% LCToolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with LCToolbox. If not, see <http://www.gnu.org/licenses/>.

function h = Hinfsyn_design_plot(G,K,WS,WU,WT)
%HINFSYN_DESIGN_PLOT Make performance plots of the H infinity design
%procedure

%% Extract tfs
G = balreal(G);
P = augw(G,1,1,1);
H = lft(P,K,1,1); % get size right
S = H(1,1); U = H(2,1); T = H(3,1);

%% Compute frfs
p = pole(K); pmax = max(abs(p)); pmin = min(abs(p));
z = zero(K); zmax = max(abs(z)); zmin = min(abs(z));

logf_range = log10([min(pmin,zmin)/50;max(pmax,zmax)*50]/(2*pi)); N = (logf_range(2)-logf_range(1))*20; % 20 points per decade
f = logspace(logf_range(1),logf_range(2),N);

if(~isempty(WS)), FRF_WS = squeeze(freqresp(WS,f,'Hz')); else FRF_WS = []; end
if(~isempty(WU)), FRF_WU = squeeze(freqresp(WU,f,'Hz')); else FRF_WU = []; end
if(~isempty(WT)), FRF_WT = squeeze(freqresp(WT,f,'Hz')); else FRF_WT = []; end

FRF_S = squeeze(freqresp(S,f,'Hz'));
FRF_U = squeeze(freqresp(U,f,'Hz'));
FRF_T = squeeze(freqresp(T,f,'Hz'));


%% Display plots
% set plotting style
color_G = [0,0,0];
color_W = [0,0,1];
linestyle_G = '-';
linestyle_W = '--';

% calculate ranges
x_range = [f(1) f(end)]; logx0 = [ceil(logf_range(1)) floor(logf_range(2))]; n = 4;
df = floor((logx0(2)-logx0(1))/(n-1)); logx0(2) = logx0(1)+(n-1)*df;
x_ticks = logspace(logx0(1),logx0(2),n);
y_range_S = [min(db(abs([FRF_S(:);1./FRF_WS(:)])))-10, max(db(abs([FRF_S(:);1./FRF_WS(:)])))+10];
y_range_U = [min(db(abs([FRF_U(:);1./FRF_WU(:)])))-10, max(db(abs([FRF_U(:);1./FRF_WU(:)])))+10];
y_range_T = [min(db(abs([FRF_T(:);1./FRF_WT(:)])))-10, max(db(abs([FRF_T(:);1./FRF_WT(:)])))+10];

h(1) = figure('Name','Closed loop FRFs');
subplot(131)
semilogx(f, db(abs(FRF_S)), 'Color', color_G, 'LineWidth',1,'LineStyle',linestyle_G)
if(~isempty(FRF_WS)), hold on, semilogx(f, db(abs(1./FRF_WS)), 'Color', color_W, 'LineWidth',1,'LineStyle',linestyle_W), hold off, end
axis([x_range y_range_S]);
set(gca,'XTick',x_ticks);
%axis([1,1000,min(db(abs(FRF_S(f>1)))),20]), grid on
%set(gca, 'YTick', [-40;-20;0], 'YTickMode', 'manual')
xlabel('f [Hz]')
ylabel('|S| [dB]')
subplot(132)
semilogx(f, db(abs(FRF_U)), 'Color', color_G, 'LineWidth',1,'LineStyle',linestyle_G)
if(~isempty(FRF_WU)), hold on, semilogx(f, db(abs(1./FRF_WU)), 'Color', color_W, 'LineWidth',1,'LineStyle',linestyle_W), hold off, end
axis([x_range y_range_U]);
set(gca,'XTick',x_ticks);
% axis([1,1000,-40,45]), grid on
% set(gca, 'YTick', (-20:20:20)', 'YTickMode', 'manual')
xlabel('f [Hz]')
ylabel('|U| [dB]')
subplot(133)
semilogx(f, db(abs(FRF_T)), 'Color', color_G, 'LineWidth',1,'LineStyle',linestyle_G)
if(~isempty(FRF_WT)), hold on, semilogx(f, db(abs(1./FRF_WT)), 'Color', color_W, 'LineWidth',1,'LineStyle',linestyle_W), hold off, end
axis([x_range y_range_T]);
set(gca,'XTick',x_ticks);
% axis([1,1000,-80,25]), grid on
% set(gca, 'YTick', (-40:40:20)', 'YTickMode', 'manual')
xlabel('f [Hz]')
ylabel('|T| [dB]')

end

