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

function [ varargout ] = param_ident( varargin )

required = {'data','settings'};
optional = {'method'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);

for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'data'); assert(~exist('n','var'),'You cannot define the system order twice.'); data = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'settings'); assert(~exist('settings','var'),'You cannot define the highest numerator order twice.'); settings = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'method'); assert(~exist('method','var'),'You cannot define the method twice.'); method = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end

if ~exist('method','var'); method = []; end;

optd = struct('numl',0,'denl',0,'maxIter',500,'relVar', 1e-20,'GN',0,'FRFW',[],'gradTol',1e-8,'sysInit',[],'Ts',0);
settings = mergestruct(settings,optd);

if isa(data,'idfrd') || isa(data,'frd') || isa(data,'FRDmod')
    if all(size(settings.numl) == 1), settings.numl = repmat(settings.numl,nout(data),nin(data)); end
    if all(size(settings.numh) == 1), settings.numh = repmat(settings.numh,nout(data),nin(data)); end
    
    freq = data.Frequency;
    if ~strcmp(data.FrequencyUnit,'Hz')
        freq = freq/(2*pi);
    end
            
    switch method
    
        case 'nllsfdi'
            FRF = permute(reshape(data.ResponseData,[],length(data.Frequency)),[2,1])

            if settings.Ts == 0
                [Bn,An,~,~,~,~] = nllsfdi(FRF, freq, settings.FRFW, settings.denh, settings.denl, settings.numh(:), settings.numl(:), settings.maxIter, settings.relVar, settings.GN, 'c');
            else
                [Bn,An,~,~,~,~] = nllsfdi(FRF, freq, settings.FRFW, settings.denh, settings.denl, settings.numh(:), settings.numl(:), settings.maxIter, settings.relVar, settings.GN, 'd', 1/settings.Ts);           
            end
            B = reshape(num2cell(Bn,2),[nout(data),nin(data)]);
            model = TFmod(B,An,settings.Ts);
            
        case 'MIMO_ML'
            assert(isa(data,'IdentFRDmod'),'MIMO_ML expects the data to be an IdentFRDmod. Try nllsfdi instead');
            
            % initialize estimation process
            if settings.Ts == 0, plantplane = 's';
            else plantplane = 'z'; end
            [Sel,~,ModelVar,IterVar] = MIMO_ML_DefaultValues(settings.denh, settings.numh, nin(data), nout(data), plantplane);
            Sel.A(1:settings.denl) = 0;
            for i=1:nout(data)
                for j=1:nin(data)
                    Sel.B(i,j,1:settings.numl(i,j)) = 0;
                end
            end
            ModelVar.Struct = 'OE'; % noise free input
            
            % data for parametric estimation
            ident_data.freq = freq;
            ident_data.Ts = settings.Ts;
            ident_data.CU = data.CYU_.n(2,2,:);
            ident_data.CYU = data.CYU_.n(1,2,:);
            ident_data.CY = data.CYU_.n(1,1,:);
            ident_data.Y = data.Y_.mean;
            ident_data.U = data.U_.mean;

            [ThetaWGTLS,smax,smin,wscale] = MIMO_WGTLS(ident_data,Sel,ModelVar);
            [ThetaBTLS,CostBTLS,smax,smin,wscale] = MIMO_BTLS(ident_data,Sel,ThetaWGTLS,ModelVar,IterVar);

            [ThetaML, CostML, smax, smin, wscale] = MIMO_ML(ident_data, Sel, ThetaBTLS, ModelVar, IterVar);
            B = reshape(ThetaML.B,[],size(ThetaML.B,3))
            B = num2cell(fliplr(B),2);
            B = reshape(B,[nout(data),nin(data)]);
            A = fliplr(ThetaML.A);
            model = TFmod(B,A,settings.Ts);
            
        case 'MIMO_NLS'
            assert(all(size(data)==1),'Only implemented for siso due to problems with Fourier coefficient on Y');
            data = IdentFRDmod(data.ResponseData,data.Frequency,'FrequencyUnit',data.FrequencyUnit);
            F = length(freq); 
            data.CYU_.n = repmat(eye(nout(data)+nin(data)),[1,1,F]);
            data.Y_.mean = data.ResponseData;
            data.U_.mean = ones(nin(data),F);
            
            if ~isempty(settings.FRFW)
                assert(all(size(settings.FRFW) == [F,1]),'MIMO_NLS expects a weight for each frequency line only');
                for k = 1:F
                    data.CYU_.n(:,:,k) = data.CYU_.n(:,:,k)/settings.FRFW(k);
                end
            end
            
            model = param_ident('data',data,'method','MIMO_ML','settings',settings);
            
        otherwise
            try
                fh = str2func(method);
                if strcmp(data.FrequencyUnit,'Hz')
                    data.Frequency = data.Frequency*(2*pi);
                    data.FrequencyUnit = 'rad/s';
                end
                model = fh(data,settings.nh);
                
            catch
                warning('No appropriate method for the given data type.')
            end
        
    end
    varargout{1} = model;

    
elseif isa(data,'iddata') || isa(data,'MeasurementData')
    
    if isa(data,'MeasurementData')   
        labels = varargin{1};
        data = iddata(getSignal(data,labels.output),getSignal(data,labels.input),Ts(data));
    end
    
    switch method
        
    otherwise
        try
            fh = str2func(method);
            model = fh(data,settings.n);
        catch
            warning('No appropriate method for the given data type.')
        end
    end
    varargout{1} = fromstd(model);
    
end

end


