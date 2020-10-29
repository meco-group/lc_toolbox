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

classdef TDMeasurementData < MeasurementData
    %MEASUREMENTDATA Class -- Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = private)
        data_            % matrix N x labels.length
        datalabels_ = {} % {strings} ?!? or should these be signals? 
        excitation_      % TimeSignal !!! currently only scalar signal !!!
        periodicity_     % if the measurement is periodic or not
    end
    
    methods
        
        function self = TDMeasurementData(varargin)
            
            required = {'label','excitation','data','datalabels','periodic'};
            optional = {'date','info'};
            
            calls = varargin(1:2:end);
            
            assert(mod(nargin,2) == 0,'Odd number of arguments.');
            assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
            assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
            
            for i = 1:nargin/2
                if strcmp(varargin{2*i-1},'label'); assert(~exist('label','var'),'You cannot define a label twice.'); label = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'excitation'); assert(~exist('excitation','var'),'You cannot define an excitation signal twice.'); excitation = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'data'); assert(~exist('data','var'),'You cannot define the data twice.'); data = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'datalabels'); assert(~exist('datalabels','var'),'You cannot define data labels twice.'); datalabels = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'periodic'); assert(~exist('periodic','var'),'You cannot define a periodicity flag twice.'); periodic = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'date'); assert(~exist('date','var'),'You cannot define a date twice.'); date = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'info'); assert(~exist('info','var'),'You cannot define info twice.'); info = varargin{2*i}; end;
                if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
            end
            
            if ~exist('date','var'); date = datestr(today); end;
            if ~exist('info','var'); info = ''; end;
            
            p = inputParser;
            addRequired(p,'label',@(x) assert(ischar(x), 'Your label should be a string.'));
            addRequired(p,'excitation',@(x) assert(isa(x,'TimeSignal'),'Your excitation signal should be a TimeSignal object.'));
            addRequired(p,'data',@(x) assert(isnumeric(x),'Your data should be numeric.'));
            addRequired(p,'datalabels',@(x) assert(iscellstr(x),'Your data labels should be provided in a cell of strings.'));
            addRequired(p,'periodic',@(x) assert(islogical(x) && isscalar(x),'The periodicity flag should be either true or false.'));
            addRequired(p,'date',@(x) assert(ischar(x) || isa(x,'datetime'),'The date should be a string or datetime object.'));
            addRequired(p,'info',@(x) assert(ischar(x) || isempty(x),'Additional information should be provided in string format.'));

            parse(p,label, excitation, data, datalabels, periodic, date, info);
            
            % assign properties
            
            msdata = {p.Results.label, p.Results.info, p.Results.date};
            self@MeasurementData(msdata{:});
            
            self.excitation_ = p.Results.excitation;
            self.data_ = p.Results.data;
            self.datalabels_ = p.Results.datalabels;
            self.periodicity_ = p.Results.periodic;
            
            
            % additional 'inter-argument' checks
            
            [~,numcols] = size(self.data_);
            assert(numcols == length(self.datalabels_), 'Number of labels does not match number of measured parameters.'); 
            
        end
        
        function disp(self)
            if ~ischar(self.date_); ds = datestr(self.date_); else ds = self.date_; end;
          	disp(['Time domain measurement data ''' self.label_ ''' recorded on ' ds ' with these properties:']);
            disp(['   excitation_ signal: ''' self.excitation_.label '''']);             
            disp(['   sample frequency:   ' num2str(self.fs) ' Hz']);  
            disp(['   number of samples:  ' num2str(self.length)]);
            if ~isempty(self.datalabels_)
                labels = '   datalabels:         ';
                for i = 1:length(self.datalabels_)
                    labels = [labels '''' self.datalabels_{i} ''' '];
                end
                disp(labels)
            end
        end
        
        function out = type(self)     % type of excitation_ 
            out = self.excitation_.type;
        end
        
        function out = length(self)   % number of time samples
            out = size(self.data_,1);
        end
        
        function out = Ts(self)       % sampling period
            out = 1/self.excitation_.fs;
        end
        
        function out = fs(self)       % sampling frequency
            out = self.excitation_.fs;
        end
                    
        function out = fmin(self)     % min excitation_ frequency
            out = self.excitation_.fmin;
        end
           
        function out = fmax(self)     % max excitation_ frequency
            out = self.excitation_.fmax;
        end
        
        function out = nop(self)     % number of periods
            T = self.excitation_.length;
            out = floor(self.length / T);
        end
        
       function [val,t,aux] = signal(self,varargin)

            if nargin == 1 || isempty(varargin{1})
                indices = (1:size(self.data_,2))';
            elseif isnumeric(cell2mat(varargin{1}))
                indices = cell2mat(varargin{1});
                assert(all(indices <= size(self.datalabels_,2)),'At least one of your column numbers was out of range.');
            else
                indices = self.getIndex(varargin{1}{:});
            end

            if nargin < 3
                mode = 'full';
            elseif strcmp(varargin{2}, 'full') || strcmp(varargin{2}, 'periodic')
                mode = varargin{2};
            else
                error('illegal mode for ''signal''; the options are ''full'' and ''periodic''')
            end

            switch mode
                case 'full'
                    t = self.Ts * (0 : self.length-1)';
                    val = self.data_(:, indices);
                    aux = zeros(0,0,0);
                case 'periodic'
                    ppp = self.excitation_.length;
                    nop  = self.nop;
                    aux = permute(reshape(self.data_(end-ppp*nop+1:end,indices)', length(indices), ppp, nop),[2,1,3]);
                    val = mean(aux,3);
                    t = self.Ts * (0 : ppp-1)';
            end
            
        end

        function plotSignal(self, varargin)
            
            mode = 'full';
            indices = (1:size(self.data_,2))';
            if ~isempty(varargin)
                if iscell(varargin{1}) && isnumeric(cell2mat(varargin{1}))
                    indices = cell2mat(varargin{1});
                    assert(all(indices <= size(self.datalabels_,2)),'At least one of your column numbers was out of range.');
                elseif ischar(varargin{1})
                    mode = varargin{1};
                    if nargin >= 3 && iscell(varargin{2})
                        if isnumeric(cell2mat(varargin{2}))    
                            indices = cell2mat(varargin{2});
                        else
                            indices = getIndex(self,varargin{2}{:});
                        end
                        assert(all(indices <= size(self.datalabels_,2)),'At least one of your column numbers was out of range.');
                    end
                else
                    indices = self.getIndex(varargin{1}{:});
                end
            end
            
            if isempty(self.datalabels_)
                legends = cellfun(@num2str, mat2cell(indices,ones(size(indices))), 'UniformOutput', false);
            else
                legends = self.datalabels_(indices);
            end
			
            assert(strcmp(mode,'full') || strcmp(mode,'periodic'), 'Illegal mode for ''plotSignal''; the options are ''full'' and ''periodic''.');
            
            [val,t,aux] = self.signal(num2cell(indices), mode);
            alls = [reshape(aux, size(aux,1),size(aux,2)*size(aux,3)), val];
            
            map = colormap(lines); map = map(1:length(indices), :);
            map = kron([0.5*ones(size(aux,3),1);1], map) + kron([0.5*ones(size(aux,3),1);0], ones(size(map)));
			figure
			h = plot(t, alls);
            for i = 1:length(h); set(h(i),'Color',map(i,:)); end
			xlabel('t'), ylabel('data'), title(['Measured data in ''' self.label_ ''''])
            legend(h(end-length(indices)+1:end),legends), axis tight
        end

        function [F,f,aux] = spectrum(self, varargin)
            if nargin == 1 || isempty(varargin{1})
                indices = (1:size(self.data_,2))';
            elseif iscell(varargin{1}) && isnumeric(cell2mat(varargin{1}))
                indices = cell2mat(varargin{1});
                assert(all(indices <= size(self.datalabels_,2)),'At least one of your column numbers was out of range.');
            else
                indices = self.getIndex(varargin{1}{:});
            end
            
            if nargin < 3 
                mode = 'full';
            elseif strcmp(varargin{2}, 'full') || strcmp(varargin{2}, 'periodic')
                mode = varargin{2};
            else
                error('illegal mode for ''spectrum''; the options are ''full'' and ''periodic''')
            end
            
            fs = self.fs;
            nop = self.nop; 
            ppp = self.excitation_.length;
            N = nop * ppp;
            
            switch mode
                case 'full'
                    if nop<=1 ; warning('Noise cannot be distinguished if you use only one period.'); end
                    f = fs/N * (0 : N-1)';
                    F = fft(self.data_(end-N+1:end,indices))/self.nop;
                    F = F(f<=fs/2,:); f = f(f<=fs/2);
                    i_tot = (1:length(f))';
                    i_per = (1:self.nop:length(f))';
                    i_nonper = setdiff(i_tot, i_per);
                    i_exc = i_tot(ismembertols(f, self.excitation_.excf,sqrt(eps)));
                    i_per_exc = intersect(i_per, i_exc);
                    i_per_nonexc = setdiff(i_per, i_per_exc);
                    if ~any(ismember(i_per_nonexc, min(i_exc):max(i_exc)));  warning('Nonlinearity in the frequency range of interest cannot be distinguished if you excite all frequencies in this range.'); end
                    aux = struct('excitation',i_per_exc,'nonlinearity',i_per_nonexc,'noise',i_nonper);
                case 'periodic'
                    f = fs/ppp * (0 : ppp-1)';
                    i_per = (1:ppp); 
                    i_exc = i_per(ismembertols(f, self.excitation_.excf, sqrt(eps)));
                    dmat = permute(reshape(self.data_(end-N+1:end,indices)', length(indices), ppp, nop),[2,1,3]);
                    Dmat = fft(dmat,[],1);
                    f = f(i_exc);
                    aux = Dmat(i_exc,:,:);
                    F = mean(aux,3);
            end
        end

        function plotSpectrum(self,varargin)
            if nargin == 1 || isempty(varargin{1})
                indices = (1:size(self.data_,2))';
            elseif isnumeric(cell2mat(varargin{1}))
                indices = cell2mat(varargin{1});
                assert(all(indices <= size(self.datalabels_,2)),'At least one of your column numbers was out of range.');
            else
                indices = self.getIndex(varargin{1}{:});
            end

            if isempty(self.datalabels_)
                legends = cellfun(@num2str, mat2cell(indices,ones(size(indices))), 'UniformOutput', false);
            else
                legends = self.datalabels_(indices);
            end
			
            if nargin < 3
                mode = 'full';
            elseif strcmp(varargin{2}, 'full') || strcmp(varargin{2}, 'periodic')
                mode = varargin{2};
            else
                error('illegal mode for ''plotSpectrum''; the options are ''full'' and ''periodic''')
            end

            switch mode
                case 'full' % full spectrum, with indices indicating excitation_, nonlinearity, noise
					[F, f, indices] = self.spectrum(num2cell(indices), mode);
					F(db(F)<-180) = NaN;  % dirty trick to prevent -Inf
					figure
                  if ~isempty(indices.nonlinearity) && ~isempty(indices.noise)
                        legends = [cellfun(@(x) ['noise ' x], legends, 'un', 0) cellfun(@(x) ['nonlinearity ' x], legends, 'un', 0) legends];
                  elseif ~isempty(indices.nonlinearity) && isempty(indices.noise)
                        legends = [cellfun(@(x) ['nonlinearity ' x], legends, 'un', 0) legends];
                  elseif isempty(indices.nonlinearity) && ~isempty(indices.noise) 
                        legends = [cellfun(@(x) ['noise ' x], legends, 'un', 0) legends];
                  end
				  try set(gca, 'ColorOrderIndex', 1); end
                  semilogx(f(indices.noise), db(F(indices.noise,:)), 'g.'), hold on, axis tight
                  semilogx(f(indices.nonlinearity), db(F(indices.nonlinearity,:)), 'ro'), try set(gca, 'ColorOrderIndex', 1); end
                  semilogx(f(indices.excitation), db(F(indices.excitation,:)), 'bx')
                  ylabel('|data| [dB]'), xlabel('f  [Hz]'), title(['Measured data in ''' self.label_ '''']);
                  legend(legends);

				case 'periodic' % only spectrum at excitation_ frequencies, possibly with spectra of all periods
					[F, f, F_all] = self.spectrum(num2cell(indices), mode);
					F(db(F)<-180) = NaN;  % dirty trick to prevent -Inf
					F_all(db(F_all)<-180) = NaN;
                    
                  alls = [reshape(F_all, size(F_all,1),size(F_all,2)*size(F_all,3)), F];
            
                  map = colormap(lines); map = map(1:length(indices), :);
                  map = kron([0.5*ones(size(F_all,3),1);1], map) + kron([0.5*ones(size(F_all,3),1);0], ones(size(map)));
                  figure
                  h = semilogx(f, db(alls));
                  for i = 1:length(h); h(i).Color = map(i,:); end
                  ylabel('|data| [dB]'), xlabel('f  [Hz]'), title(['Measured data in ''' self.label_ ''''])
                  legend(h(end-length(indices)+1:end),legends), axis tight
            end
       end

        function chkinpconsist(self,label)
            
            index = find(ismember(self.datalabels_,label));
            assert(length(index)==1,'Your label was not recognized.');
            data1 = self.data_(:,index);
            data1 = data1/abs(max(data1)); % normalize
            [~,data2] = self.excitation_.signal;
            
            if self.periodicity_ == 1  % repeat excitation signal
                r = ceil(length(data1)/length(data2));
                data2 = repmat(data2,r);
                data2 = data2(1:length(data1));
                data2 = data2'/abs(max(data2));
                shift = -finddelay(data1,data2);
            else if self.periodicity_ == 0 % fill with zeros
                     tmp = zeros(length(data1),1);
                     tmp(1:length(data2)) = data2';
                     data2 = tmp;
                     shift = -finddelay(data1,data2);
                 else
                     error('I ended up in an undefined state. Shouldn''t happen.');
                 end
            end
            
            data2 = circshift(data2,shift);
            
            figure;
            subplot(211)
            plot(data1); hold on; plot(data2);
            title(['Visual check for similarity of ''' label ''' and ''' self.excitation_.label ''' (normalized)']);
            set(gca,'xtick',[]);
            xlim([1 length(data1)]);
            ylim([-1 1]);
            legend(['recorded data ''' label ''''],['signal data ''' self.excitation_.label '''']);
            subplot(212)
            title(['Difference between ''' label ''' and ''' self.excitation_.label ''' (normalized)']);
            set(gca,'xtick',[]);
            if self.periodicity_ == 1
                plot(data1-data2);
                xlim([1 length(data1)]);
                ylim([-2 2]);
            else if self.periodicity_ == 0
                    tmp = data1-data2;
                    if shift < 0
                        shift = length(data1)-shift;
                    end
                    plot(shift+1:self.excitation_.length+shift,tmp(shift:self.excitation_.length+shift));
                    xlim([1 length(data1)]);
                    ylim('auto');
                 else
                    error('I ended up in an undefined state. Shouldn''t happen.');
                 end
            end
            legend(['difference between ''' label ''' and ''' self.excitation_.label '''']);
            
        end
        
        function self = detrend(self, labels)
            
            if nargin == 1
                indices = (1:size(self.data_,2))'; %if labels not provided, we detrend all channels
            else
                assert(sum(ismember(self.datalabels_,labels))==length(labels),['At least one of the labels your provided was not a valid data label. Valid labels are: ' strjoin(self.datalabels_)]);
                indices = find(ismember(self.datalabels_,labels));
            end
            self.data_(:,indices) = detrend(self.data_(:,indices));     
            
        end
        
        function labels = getDataLabel(self, varargin)
           % return the labels corresponding to data column numbers
           assert(nargin > 1,'You should provide at least one column number.');
           assert(all(cellfun(@(x) (isnumeric(x) && isscalar(x) && mod(x,1) == 0 && x >= 1),varargin)), 'Column numbers can only be positive integers.');
           assert(all(cellfun(@(x) x <= length(self.datalabels_),varargin)),'At least one of your column numbers was out of range.');
           
           labels = cell(nargin-1,1);
           for i = 1:(nargin-1)
                labels{i} = self.datalabels_{varargin{i}};
           end
        end
        
        function idx = getIndex(self,varargin)
           % return the data column numbers corresponding to labels 
            assert(all(cellfun(@(x) ischar(x),varargin)), 'Labels can only be of type string.');
            assert(all(cellfun(@(x) ismember(x,self.datalabels_),varargin)), ['I could not recognize some of your labels. The possibilities are : ' strjoin(self.datalabels_)]);
            idx = zeros(nargin-1,1);
            for i = 1:(nargin-1)
                idx(i) = find(ismember(self.datalabels_,varargin{i}));
            end
        end
        
        function self = clip(self, mode, varargin) % create limited view on data
           % mode: sample range -- time range -- number of periods (first so many / last so many)
           
           assert(ischar(mode),'Your mode was not recognized. The possible modes are: ''sample'', ''value'', ''firstnper'' and ''lastnper''.');
           
           if strcmp(mode,'sample')
               assert(nargin >= 3,'Not enough input arguments for mode ''sample''.');
               val = varargin{1};
               assert(all(val > 0) && all(mod(val,1) == 0) && isnumeric(val) ...
                   && (all(size(val)-[1 2]==0) || all(size(val)-[2 1]==0)) && val(2) >= val(1) ...
                   && size(self.data_,1) >= val(2),'The sampling range you provided is not valid.');
               
               self.data_ = self.data_(val(1):val(2),:);
               self.info_ = [self.info_ ' // Clipped from ''' self.label_ ''': samples ' num2str(val(1)) ' to ' num2str(val(2)) '.'];
               self.label_ = [self.label_ '_clipped'];
               
           elseif strcmp(mode,'value')
               assert(nargin >= 4,'Not enough input arguments for mode ''value''.');
               if ischar(varargin{1});
                    col = self.data_(:,self.getIndex(varargin{1}));
               elseif isnumeric(varargin{1})
                    assert(varargin{1} <= size(self.datalabels_,2),'Your column number was out of range.');
                    col = self.data_(:,varargin{1});
               else
                   error('Invalid input for mode ''value''.');
               end
               val = varargin{2};
               assert(isnumeric(val) && (all(size(val)-[1 2]==0) || all(size(val)-[2 1]==0)) ...
                   && val(1) >= min(col) && max(col) >= val(2),'The range you provided is not valid.');
               
               [~,imin] = min(abs(col-val(1)));
               [~,imax] = min(abs(col-val(2)));
               imax = max(imin, imax);

               self.data_ = self.data_(imin:imax,:);
               self.info_ = [self.info_ ' // Clipped from ''' self.label_ ''': ''' varargin{1} ''' values from ' num2str(col(imin)) ' to ' num2str(col(imax)) '.'];
               self.label_ = [self.label_ '_clipped'];

               if ~all(diff(col) >= 0)
                   warning('The data which you used for clipping is not monotonically increasing. This might lead to unexpected results.');
               end
               
               if(col(imin)~=val(1) || col(imax~=val(2)))
                   warning(['Could not find exact matches of your interval with the data. The closest match was: [' num2str(col(imin)) ' , ' num2str(col(imax)) '].']);
               end

           elseif strcmp(mode,'firstnper');
               assert(nargin >= 3,'Not enough input arguments for mode ''firstnper''.');
               assert(self.periodicity_ == 1,'You cannot ask for periods if your measurement was not periodic.');
               val = varargin{1};
               assert(isscalar(val) && isnumeric(val) && mod(val,1) == 0,'You can only ask for an integer number of periods.');
               p = self.excitation_.length;
               nofp = floor(self.length/self.excitation_.length);
               
               if val > nofp
                   warning(['You asked more periods than your measurement actually contains. I returned the first ' num2str(nofp) ' periods instead.']);
                   val = nofp;
               end
               
               self.data_ = self.data_(1:p*val,:);
               self.info_ = [self.info_ ' // Clipped from ''' self.label_ ''': first ' num2str(val) ' periods.'];
               self.label_ = [self.label_ '_clipped'];
               
           elseif strcmp(mode,'lastnper');
               assert(nargin >= 3,'Not enough input arguments for mode ''lastnper''.');
               assert(self.periodicity_ == 1,'You cannot ask for periods if your measurement was not periodic.');
               val = varargin{1};
               assert(isscalar(val) && isnumeric(val) && mod(val,1) == 0,'You can only ask for an integer number of periods.');
               p = self.excitation_.length;
               nofp = floor(self.length/self.excitation_.length);
               
               if val > nofp
                   warning(['You asked more periods than your measurement actually contains. I returned the last ' num2str(nofp) ' periods instead.']);
                   val = nofp;
               end
               
               self.data_ = self.data_(end-p*val+1:end,:);
               self.info_ = [self.info_ ' // Clipped from ''' self.label_ ''': last ' num2str(val) ' periods.'];
               self.label_ = [self.label_ '_clipped'];
               
           else 
               error('Your mode was not recognized. The possible modes are: ''sample'', ''value'', ''firstnper'' and ''lastnper''.');
           end
        end
        
        function out = split(self)
                % splits the measurementdata in its periods
                if mod(self.length(), self.excitation_.length())~=0
                    warning(['Your measurement data does not contain an integer number of periods. I will ignore the first incomplete period.']);
                end
                for i = 1:nop(self)
                    out{nop(self)+1-i} = clip(self, 'sample', [self.length()-i*self.excitation_.length()+1, self.length()-(i-1)*self.excitation_.length()]);
                end
        end
           
    end
end

