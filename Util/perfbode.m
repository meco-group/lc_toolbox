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

function [h] = perfbode(G,W,varargin)
%PERFBODE 
%   Make bode diagram of the function G, and weighting functions W.
%   inputs:
%   G = transferfunction
%   W = Weighting function on G. (W = [] also supported)
%   options = structure specifying plotting style options
%       color_G,color_W: plotting color of G and W (default
%       [0,0,0],[0,0,1])
%       linestyle_G,linestyle_W: linestyle of G and W (default '-','--')
%       marker_G,marker_W: marker style (default 'n','n' <= no marker)
%       ppdec: points per decade (default 20)
%       legendlabel: label for the legends (default ignored)
%       tfname: label on y-axis (default ignored)
%   h = axes handles on which to plot the new graphs
%   output:
%   h = axes handles (can be a new object or the h passed to the function)

options = struct('ppdec',20,...
                 'legendlabel','',...
                 'tfname','',...
                 'plotmode',0,...
                 'draw_style', struct('color_G',[0,0,0],...
                                      'color_W',[0,0,1],...
                                      'linestyle_G','-',...
                                      'linestyle_W','--',...
                                      'marker_G','none',...
                                      'marker_W','none'));

if(nargin>2)    
    if(~isempty(varargin{1}))
        options = mergestruct(varargin{1},options);
    end
end
if(nargin>3)
    h = varargin{2};
    if(~isempty(get(h,'Children'))) % Hold plot on if there are already plots on the axis
        hold(h,'on');
    end
else
    fig = figure();
    if(options.plotmode == 1)
        p = get(fig, 'Position');
        p(3:4) = [300,600];
        set(fig, 'Position', p);
    end
    
    h = axes('Parent',fig);
    axis(h,[1,2,1,2]);
end

%% Check range
% G = minrealr(G);
% G = std(G);
ap = pole(G);
az = zero(G);
if(isct(G))
    ap = abs(horzcat(ap{:}))/(2*pi);
    az = abs(horzcat(az{:}))/(2*pi);
    pmax = max(ap(ap<1e10)); pmin = min(ap(ap>1e-6));
    zmax = max(az(az<1e10)); zmin = min(az(az>1e-6));   
    logf_range = [-6;10];
else
    ap = abs(horzcat(ap{:}))/(G.Ts*2*pi);
    pmax = max(ap(ap<1e10)); pmin = min(ap(ap>1e-6));
    az = abs(horzcat(az{:}))/(G.Ts*2*pi);
    zmax = max(az(az<1e10)); zmin = min(az(az>1e-6));
    logf_range = [-6;log10(1/(G.Ts*2))];
end

% logf_range = log10([min([pmin zmin])/50;max([pmax zmax])*50]/(2*pi)); N = (logf_range(2)-logf_range(1))*(options.ppdec); % 20 points per decade
% f = logspace(logf_range(1),logf_range(2),N);
N = round((logf_range(2)-logf_range(1))*(options.ppdec));
f = logspace(logf_range(1),logf_range(2),N);
sigma = zeros(1,length(f));

%% Compute frequency responses
FRF_G = freqresp(G,f,'Hz');
size_G = size(FRF_G);
if(~(all(size_G(1:2)==1)))
    for(k = 1:size(FRF_G,3))
        [~,sigmas,~] = svd(FRF_G(:,:,k));
        n = min(size(sigmas));
        sigmas = sigmas(1:n,1:n);
        sigma(k) = max(diag(sigmas));
    end
    FRF_G = sigma;
else
    FRF_G = squeeze(FRF_G);
end

if(~isempty(W))
    FRF_W = zeros(length(W),length(f));
    for(j=1:length(W))
        v = freqresp(W{j},f,'Hz');
        size_W = size(v);
        if(~(all(size_W(1:2)==1)))
            for(k = 1:size(v,3))
                [~,sigmas,~] = svd(v(:,:,k)); 
                sigma(k) = max(sigmas(:));
            end
            FRF_W(j,:) = sigma;
        else
            FRF_W(j,:) = squeeze(v);
        end
    end
end

%% Display plots
% calculate ranges
x_range = [min([pmin zmin]);max([pmax zmax])]; 

a = axis(h);
a(1) = 10^(floor(log10(min(x_range(1)*0.1,a(1))))); a(2) = 10^(ceil(log10(max(x_range(2)*10,a(2))))); %(10-9*isdt(G))

logx0 = [ceil(log10(a(1))) floor(log10(a(2)))]; n = 4;
df = floor((logx0(2)-logx0(1))/(n-1)); logx0(2) = logx0(1)+(n-1)*df;

if(logx0(1)~=logx0(2))
    x_ticks = logspace(logx0(1),logx0(2),n);
else
    m_logx0 = mean(logx0);
    x_ticks = logspace(m_logx0-1,m_logx0+n-2,n);
end

y_range_G = [min(db(abs(FRF_G(:))))-15, max(db(abs(FRF_G(:))))+15];
a(3) = min(y_range_G(1),a(3)); a(4) = max(y_range_G(2),a(4));

% Plot G
p = semilogx(h,f, db(abs(FRF_G)), 'Color', options.draw_style.color_G, 'LineWidth',1,'LineStyle',options.draw_style.linestyle_G); hold(h,'off');
d = ceil(length(f)/10); fm = f(ceil(d/2):d:end); 
if(options.draw_style.marker_G ~= 'n')
    hold(h,'on');
    FRF_Gm = FRF_G(ceil(d/2):d:end,:); 
    p = semilogx(h,fm, db(abs(FRF_Gm)), 'Color', options.draw_style.color_G, 'LineWidth',1,'LineStyle',options.draw_style.marker_G); hold(h,'off');
end

if(isdt(G))
    hold(h,'on');
    semilogx([1/(2*G.Ts) 1/(2*G.Ts)],[-400;400],'k');
end

% Plot W
if(~isempty(W))
    hold(h,'on');
    semilogx(h,f, db(abs(1./FRF_W)), 'Color', options.draw_style.color_W, 'LineWidth',1,'LineStyle',options.draw_style.linestyle_W), hold(h,'off')
    if(options.draw_style.marker_W ~= 'n')
        hold(h,'on');
        FRF_Wm = FRF_W(:,ceil(d/2):d:end);
        semilogx(h,fm, db(abs(1./FRF_Wm)), 'Color', options.draw_style.color_W, 'LineWidth',1,'LineStyle',options.draw_style.marker_W), hold(h,'off')
    end
end

% set(h,'XTick',x_ticks);

xlabel(h,'f [Hz]')
ylabel(h,sprintf('%s [dB]',options.tfname))
if(~strcmp(options.legendlabel,''))
    l = findobj(get(h,'Parent'),'Type','axes','Tag','legend');
    if(~isempty(l))
        [~,~,plots,strings] = legend(l);
        plots(length(plots)+1) = p;
        strings{length(strings)+1} = options.legendlabel;
        delete(l);
        legend(h,plots,strings,'Location','NorthOutside','Orientation','horizontal');
    else
        legend(h,p,options.legendlabel,'Location','NorthOutside','Orientation','horizontal');
    end
end

axis(h,a);
end

