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

classdef (InferiorClasses = {?LPVLFTmod,?LTILFTmod,?LPVDSSmod,?LTIDSSmod,?FRDmod}) Umod < Model
    %UMOD Uncertain model class
    %   Uncertain model. Some of the underlying model's channels are
    %   uncertain (the dimensions correspond to delta_).
    %   Input arguments:
    %   Model,(type,W1,W2)
    
    properties
        M_               % plant with additional uncertainty channel(s) (ss)
        delta_           % normalized uncertainty (uss or umat)
    end
    
    methods
        function self = Umod(varargin)
            if isa(varargin{1},'uss')
                [M, delta] = lftdata(varargin{1});
                self.M_ = fromstd(M);
                self.delta_ = delta;
            elseif isa(varargin{1},'Model')
                if isa(varargin{1},'Umod')
                    error('Umod of Umod not implemented for now');                    
                else
                    model = varargin{1};
                    switch(nargin)
                        case 1
                            self.M_ = model;
                            self.delta_ = umat([]);
                        case 2
                            self.M_ = model;
                            self.delta_ = varargin{2};
                        case 4
                            type = varargin{2}; assert(any(validatestring(type,{'InputMult','OutputMult','Additive'})),'The uncertainty type should be of type ''InputMult'', ''OutputMult'' or ''Additive''.');
                            W1 = varargin{3}; assert(isa(W1,'numlti') || isa(W1,'AbstractLTImod') || isnumeric(W1),'Your weight is not an LTI system.');
                            W2 = varargin{4}; assert(isa(W2,'numlti') || isa(W2,'AbstractLTImod') || isnumeric(W2),'Your weight is not an LTI system.');
                            assert(size(W2,1) == size(W1,2),'Uncertainty should be square');
                            usz = size(W2,1);
                            
                            switch type
                                case 'InputMult'
                                    self.M_ = [[zeros(usz) W2]; model*[W1,eye(size(W1,1))]];
                                    norm = ucomplex('delta',0,'Radius',1);
                                    self.delta_ = norm*eye(usz);

                                case 'OutputMult'
                                    self.M_ = [[zeros(usz);W1],[W2;eye(size(W2,2))]*model];
                                    norm = ucomplex('delta',0,'Radius',1);
                                    self.delta_ = norm*eye(usz);  

                                case 'Additive'
                                    self.M_ = [zeros(usz) W2 ; W1 model];
                                    norm = ucomplex('delta',0,'Radius',1);
                                    self.delta_ = norm*eye(usz);  

                            % TO DO: add inverse I/O multiplicative
                            % uncertainty, other types, ... 

                            end
                        otherwise
                            error('The syntax for Umod is: Umod(uss) or Umod(Model, type, W1, W2).');
                    end
                end
            end
        end
        
        % dimensions
        function out = nout(self) 
            out = size(self.M_,1) - size(self.delta_,2);
        end
        
        function out = nin(self) 
            out = size(self.M_,2) - size(self.delta_,1);
        end
        
        function [dout,din] = usize(self)
            [dout,din] = size(self.delta_);
        end
        
        % standard MATLAB format
        function out = std(self)
            out = lft(std(self.M_),self.delta_);
        end
        
        function mod = evalme(self,delta)
            mod = lft(fromstd(delta),self.M_);
        end
        
        % calculate LFT
        function self = lft(self,other,nu,ny)
            if ~isempty(other) % return other when empty
                if nargin == 2
                    if all(size(other) <= size(self))
                        [nu,ny] = size(other);
                    else
                        [nu,ny] = size(self);
                    end
                end
                if ~isa(self,'Umod'), self = Umod(self); end
                if ~isa(other,'Umod'), other = Umod(other); end
                % Rearrange other uncertainty to bottom
                [uouts,uins] = usize(self);
                [uouto,uino] = usize(other);
                idxout = [(uouto+1):size(other.M_,1),1:uouto];
                idxin = [(uino+1):size(other.M_,2),1:uino];
                otherM = other.M_(idxout,idxin);
                % Rearrange other uncertainty back to top
                M = lft(self.M_,otherM,nu,ny);
                
                uout1 = 1:uouts;
                cout1 = (uouts+1):(size(self.M_,1)-ny);
                cout2 = size(self.M_,1)-ny+(1:(size(other.M_,1)-nu-uouto));
                uout2 = (size(self.M_,1)+size(other.M_,1)-nu-ny-uouto)+(1:uouto);
                uin1 = 1:uins;
                cin1 = (uins+1):(size(self.M_,2)-nu);
                cin2 = size(self.M_,2)-nu+(1:(size(other.M_,2)-ny-uino));
                uin2 = (size(self.M_,2)+size(other.M_,2)-nu-ny-uino)+(1:uino);
                
                idxout = [uout1,uout2,cout1,cout2];
                idxin = [uin1,uin2,cin1,cin2];
                
                M = M(idxout,idxin);
                delta = blkdiag(self.delta_,other.delta_);

                self = Umod(M,delta);
            end
        end
        
        % Worst-case sigma plot
        
        function wcsigma(self,varargin)
            assert(exist('wcsigma','builtin'),'wcsigma is only supported in MATLAB R2016b or higher.');
            wcsigma(self.std,varargin{:});
        end
        
        function blk = blkdiag(varargin)
            if all(cellfun(@(x)isa(x,'Umod'),varargin))
                Ms = cellfun(@(x)x.M_,varargin,'UniformOutput', false);
                M = blkdiag(Ms{:});
                
                US = cell2mat(cellfun(@usize,varargin,'un',0));
                Sout = cell2mat(cellfun(@(x)nout(x.M_),varargin,'un',0));
                Sin = cell2mat(cellfun(@(x)nin(x.M_),varargin,'un',0));
                % pick up from here
                USo = US;
                USi = US;
                for i = 2:nargin
                    USo(i:end) = USo(i:end) + Sout(i-1:end-1);
                    USi(i:end) = USi(i:end) + Sin(i-1:end-1);
                end
                sum_o = 1:sum(Sout);
                sum_i = 1:sum(Sin);
                outlist = [USo setdiff(sum_o,USo)];
                inlist = [USi setdiff(sum_i,USi)];
                outlist(outlist==0)=[];inlist(inlist==0)=[];
                M = M(outlist,inlist);
                del = cellfun(@(x)x.delta_,varargin,'UniformOutput', false);
                Delta = blkdiag(del{:});
                blk = Umod(M,Delta);

            else
                blk = varargin{1};
                for i = 2:nargin
                    blk = blkdiag_model(blk,varargin{i});
                end
            end
        end
        
        % Bode diagram
        % especially for MIMO systems, this implementation is terribly slow
        % maybe we should reimplement wcgain ourselves?
        function ubode(self,varargin)
            
            warning('The current implementation for fancy Bode diagrams is very slow. Faster alternatives are bode(std(...)), sigma(...) and wcsigma(...).');
            
            % get frequency response for each I/O channel
            
                model = self.std;
                
                if nargin==2 
                    wout = varargin{1};
                else
                    [~,~,wout] = bode(model);
                end

                % elementwise inversion for calculating the lower bound
                for i = 1:size(model,1)
                    for j = 1:size(model,2)
                        try
                        ew_inv(i,j) = minreal(inv(model(i,j)));
                        catch  % minreal only for proper systems
                        ew_inv(i,j) = inv(model(i,j));
                        end
                    end
                end

                % actual calculation of the worst-case scenario for each
                % I/O channel
                for k = 1:size(model,1)
                    for l = 1:size(model,2)
                          this_mod = model(k,l);
                          this_inv = ew_inv(k,l);
                          for i = 1:length(wout)
                              [~,~,info] = wcgain(this_mod,wout(i));
                              wcuparam = info.BadUncertainValues; 
                              try
                                respu(k,l,i) = evalfr(minreal(usubs(this_mod,wcuparam),[],false),1i*wout(i));
                              catch % minreal only for proper systems
                                respu(k,l,i) = evalfr(usubs(this_mod,wcuparam),1i*wout(i));
                              end
                              [~,~,info] = wcgain(this_inv,wout(i));
                              wclparam = info.BadUncertainValues; 
                              try
                                 respl(k,l,i) = evalfr(minreal(usubs(this_mod,wclparam),[],false),1i*wout(i));
                              catch % minreal only for proper systems
                                 respl(k,l,i) = evalfr(usubs(this_mod,wclparam),1i*wout(i));
                              end
                              respn(k,l,i) = evalfr(this_mod.NominalValue,1i*wout(i));
                          end
                    end
                end

                % plot all stuff
                figure;
                width = size(model,2);
                height = 2*size(model,1);
                wout = wout/2/pi; % convert to Hz
                
                respl(abs(respl) < 1e-8) = min(0.1*abs(respn)); % the grey area cannot be drawn to zero on a logarithmic scale
                
                for k = 1:size(model,1)
                    for l = 1:size(model,2)
                        
                        subplot(height,width,2*width*(k-1)+l);
                        set(gca,'XScale','log');
                        hold on; fill([wout ; flipud(wout)],[db(squeeze(respl(k,l,:))) ; flipud(db(squeeze(respu(k,l,:))))],[0.8 0.8 0.8],'EdgeColor','None');
                        semilogx(wout, db(squeeze(respn(k,l,:))),'b','LineWidth',1); 
                        xlabel('frequency (Hz)');
                        ylabel(['from input ' num2str(k) ' \newline magnitude (dB)']);
                        title(['to output ' num2str(l)],'FontWeight','Normal');
                        xlim([min(wout),max(wout)]);
                        ax = gca;
                        limits = ax.YLim;
                        ylim([min(db(squeeze(respl(k,l,:)))) limits(2)]);
                        
                        subplot(height,width,2*width*(k-1)+l+width);
                        set(gca,'XScale','log');
                        % ignore the phase for the time being
                        %hold on; fill([wout ; flipud(wout)],180/pi*[angle(squeeze(respl(k,l,:))) ; flipud(angle(squeeze(respu(k,l,:))))],[0.8 0.8 0.8],'EdgeColor','None');
                        semilogx(wout, 180/pi*angle(squeeze(respn(k,l,:))),'b','LineWidth',1); 
                        xlabel('frequency (Hz)');
                        ylabel('phase (�)');
                        xlim([min(wout),max(wout)]);
                        
                    end
                end

        end
    end
    
    methods (Access = protected)
        function submod = submodel(self,idxout,idxin)
            usz = usize(self);
            M = self.M_([1:usz,usz+idxout],[1:usz,usz+idxin]);
            submod = Umod(M,self.delta_);
        end
    end
        
end

