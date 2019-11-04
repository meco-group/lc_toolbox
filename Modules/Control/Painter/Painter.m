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

classdef Painter
    %PAINTER Drawing of TFs
    %   Class to handle drawing of performance plots for controller design
    
    properties
        channels = struct('io',{},'name',{},'data',struct('frf',{},'weights',{},'ctrl_prob_ID',{}))
        ctrl_prob_data = struct('name',{}, ...
                                'frf',{},...
                                'cl',{},...
                                'draw_style',{});

        style = PaintStyleConference();
        canvas_size = [600,450];
        controller_canvas_size = [600,450];
    end
    
    methods
        function obj = Painter(varargin)
            if nargin > 0
                obj.setStyle(varargin{1});
            end
        end
                
        function obj = addCtrlProb(obj,ctrl_prob_list)
        % Adds a complete control problem or a list of problems to the channels array
            for ctrl_prob = ctrl_prob_list' % row vector needed to iterate over 
                if(isa(ctrl_prob,'ControllerAnalysis') || ctrl_prob.solver.solved) % Check if the control problem has been solved correctly
                    ctrl_prob_ID = length(obj.ctrl_prob_data)+1;
                    obj.ctrl_prob_data(ctrl_prob_ID,1)...
                        = struct('name',ctrl_prob.name,...
                                 'frf',ctrl_prob.sys_K,...
                                 'cl',[],... % ctrl_prob.sys_H
                                 'draw_style',struct('color_G',[0 0 1],...
                                                     'color_W',[0 0 0],...
                                                     'linestyle_G','-',...
                                                     'linestyle_W','--',...
                                                     'marker_G','none',...
                                                     'marker_W','none'));
                    [new_channels] = getChannels(ctrl_prob);
                    for k = 1:length(new_channels)
                        obj = addChannel(obj,new_channels{k},ctrl_prob_ID);
                    end
                end
            end
        end
        
        function obj = addChannel(obj,new_channel,ctrl_prob_ID)
        % Adds a single instance of performance measures to one of the
        % channels
            channelID = findChannel(obj,new_channel.io);
            if(channelID>0) %this means there is a channel with the same in and outputs
                obj.channels(channelID,1).data(length(obj.channels(channelID,1).data)+1,1) = ...
                    struct('frf',new_channel.frf,'weights',{new_channel.weights},'ctrl_prob_ID',ctrl_prob_ID);
                if(isempty(obj.channels(channelID,1).name))
                    obj.channels(channelID,1).name = new_channel.name;
                end
            else
                channelID = length(obj.channels)+1;
                obj.channels(channelID,1).io = new_channel.io;
                obj.channels(channelID,1).name = new_channel.name;
                obj.channels(channelID,1).data(1) = struct('frf',new_channel.frf,'weights',{new_channel.weights},'ctrl_prob_ID',ctrl_prob_ID);
            end
        end  
        
        function channelID = findChannel(obj,io)
        % Searches the channel in the channels array based on provided io
            channelID = 0;
            if(~isempty(obj.channels))
                for k = 1:length(obj.channels)
                    if(isequal(obj.channels(k,1).io{1,1},io{1,1})&&isequal(obj.channels(k,1).io{1,2},io{1,2}))
                        channelID = k;
                        break;
                    end
                end
            end
        end
                
        function fig = renderControllers(obj,varargin)
            % make figure
            fig = figure('Name','Controller');
            p = get(fig, 'Position');
            set(fig, 'Position', [p(1:2) obj.controller_canvas_size]);
            
            % plot bode diagram of all controllers
            N = length(obj.ctrl_prob_data);
            labels = cell(1,N);
            args = cell(1,2*N);
            for k = 1:N
                labels{k} = obj.ctrl_prob_data(k,1).name;
                args{2*k-1} = obj.ctrl_prob_data(k,1).frf;
                args{2*k} = obj.style.channelstyle(k);
            end
            bode(args{:});
            legend(labels,'Orientation','horizontal','position',[0.45 0.95 0.1 0.05]);
            
            % set style of all graphs - nasty findobj because bode does not
            % allow easy property access
%             av = findobj(fig,'Type','axes');
%             for a = av'
%                 l = findobj(a,'Type','line');
%                 for il = 0:(k-1)
%                     set(l(end-il),'Color',obj.ctrl_prob_data(1+il,1).draw_style.color_G);
%                     set(l(end-il),'LineStyle',obj.ctrl_prob_data(1+il,1).draw_style.linestyle_G);
%                     set(l(end-il),'Marker',obj.ctrl_prob_data(1+il,1).draw_style.marker_G);
%                 end
%             end
        end
        
        function fig = renderChannels(obj,varargin)
            % Render all performance plots
            if nargin>1
                channelIDs = varargin{1};
            else
                channelIDs = 1:length(obj.channels);
            end
            
            for channelID = channelIDs
                % Set channel name
                if(~isempty(obj.channels(channelID,1).name))
                    channelname = obj.channels(channelID,1).name;
                else
                    channelname = sprintf('Channel %d',channelID);
                end
                    
                % Make figure
                fig = figure('Name',channelname);
                p = get(fig, 'Position');
                set(fig, 'Position', [p(1:2) obj.style.size()]);
                
                % Get data
                data = obj.channels(channelID,1).data;
                N = length(data);
                % Make argument lists
                syss = cell(1,2*N); weights = {}; npfrfs = {}; labels = cell(1,N);
                for k = 1:N
                    % NPFRF
                    if data(k,1).frf.hasnpfrf()
                        npfrfs{end+1} = data(k,1).frf.npfrf;
                        data(k,1).frf.npfrf = [];
                        npfrfs{end+1} = obj.style.npfrfstyle(k);
                    end
                    % SYSS
                    syss{2*k-1} = data(k,1).frf;
                    syss{2*k} = obj.style.channelstyle(k);
                    % WEIGHTS
                    for l = 1:length(data(k,1).weights)
                        if all(size(data(k,1).weights{l}) == [1,1])
                            % Try to invert the model only for siso
                            weights{end+1} = inv(data(k,1).weights{l});
                        else
                            % Compute the singular values on frequency grid
                            w = 2*pi*logspace(-8,8,160);
                            sv = sigma(data(k,1).weights{l},w);
                            weights{end+1} = frd(1./sv(1,:),w);
                        end
                        weights{end+1} = obj.style.weightstyle(l);
                    end
                    labels{k} = obj.ctrl_prob_data(data(k,1).ctrl_prob_ID,1).name;
                end
                bodemag(syss{:},weights{:},npfrfs{:});
                title(obj.channels(channelID,1).name);
                legend(labels{:});
            end
        end
        
        function fig = renderClosedloop(obj)
            fig = figure('Name','Closed Loop');
            p = get(fig, 'Position');
            set(fig, 'Position', [p(1:2) obj.controller_canvas_size]);

            Nu = size(obj.ctrl_prob_data(1,1).cl.b,2);
            Ny = size(obj.ctrl_prob_data(1,1).cl.c,1);
            
            h = zeros(Ny,Nu); p = 1;
            for in = 1:Nu
                for out = 1:Ny
                    h(out,in) = subplot(Ny,Nu,p);
                    axis(h(out,in),[1,2,1,2]); % Set some value so that it is easy to change it later
                    p = p+1;
                end
            end
            
            legend_data = cell(length(obj.ctrl_prob_data),1);
            for k = 1:length(obj.ctrl_prob_data)
                legend_data{k,1} = obj.ctrl_prob_data(k,1).name;
                options = struct('draw_style',obj.ctrl_prob_data(k,1).draw_style);
                for in = 1:Nu
                    for out = 1:Ny
                        h(out,in) = perfbode(obj.ctrl_prob_data(k,1).cl(out,in),[],options,h(out,in));
                    end
                end
            end
            legend(legend_data,'Orientation','horizontal','position',[0.45 0.95 0.1 0.05]);
        end
        
        function obj = setStyle(obj,style)
            switch style
                case 'conference'
                    obj.style = PaintStyleConference();
                case 'normal'
                    obj.style = PaintStyleNormal();
                case 'paper_plain'
                    obj.style = PaintStylePaperPlain();
                case 'paper_fancy'
                    obj.style = PaintStylePaperFancy();
                otherwise
                    error(['PaintStyle ''' style ''' not known']);
            end
        end        
    end
    
end

