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

function [model, diag, allFRFmods] = nonpar_ident(varargin)
  
% data data data ... labels method (options)


% method
% labels: input output scheduling  
% model = frd model with Input Name and  

% diag = all FRFs, (crbounds), sample total variance, sample noise variance 

if iscell(varargin{1})
    if nargin < 3 || (~isa(varargin{1}{1}, 'TDMeasurementData') && ~isa(varargin{1}{1}, 'FDMeasurementData'))
        error('Incorrect/unsupported call to nonpar_ident: should be nonpar_ident(TDMeasurementData,..., labels, method, (options))')
    end
end
    
if isstruct(varargin{end})
    options = varargin{end};
    method = varargin{end-1};
    labels = varargin{end-2};
    ndata = nargin-3;
elseif ischar(varargin{end})
    options = {};
    method = varargin{end};
    labels = varargin{end-1};
    ndata = nargin-2;
else
    error('Incorrect/unsupported call to nonpar_ident: should be nonpar_ident(TDMeasurementData,..., labels, method, (options))')
end

if ~isfield(labels, 'input') || ~(iscellstr(labels.input)||ischar(labels.input)) || ~isfield(labels, 'output') || ~(iscellstr(labels.output)||ischar(labels.output)) 
    error('Incorrect parsing of labels struct: should contain fields ''input'' and ''output'' and can contain ''scheduling'' ')
end

if isfield(labels, 'scheduling') && ~isempty(labels.scheduling)
    error('Too bad, nonpar_ident not yet implemented for LPV')
end

if ~iscell(labels.input);  labels.input = {labels.input}; end
if ~iscell(labels.output); labels.output = {labels.output}; end

data = cell(ndata,1);
for i = 1:ndata
    if iscell(varargin{i})
        for j = 1:numel(varargin{i}) 
            if ~isa(varargin{i}{j}, 'TDMeasurementData') && ~isa(varargin{i}{j}, 'FDMeasurementData')
                error('Incorrect data type!')
            else
                data{i}{j} = varargin{i}{j};
                assert(all(ismember([labels.input(:);labels.output(:)], data{i}{j}.datalabels_)), 'input or output labels not present in datalabels of measurement data');
            end
        end
    else
        data{i}{1} = varargin{i};   
    end
end

switch method
    case 'time2frf'
        [~,f] = data{1}{1}.spectrum(labels.input,'periodic');
        FRFs = zeros(length(f), ndata); % mean for each realization 
        noiseVar = zeros(length(f), ndata); 
        np = zeros(ndata,1); % total number of periods for each realization
     
        for i = 1:ndata
            FRFmp = cell(1,numel(data{i}));
            for j = 1:numel(data{i})
                [U, f1, U_all_j] = data{i}{j}.spectrum(labels.input , 'periodic');
                [Y, f1, Y_all_j] = data{i}{j}.spectrum(labels.output, 'periodic');
                assert( all(f == f1), 'measurements with different ranges of excitation frequencies') 
                np(i) = np(i) + size(U_all_j,3);
                FRFs(:,i) = FRFs(:,i) + Y./U;
                FRFmp{j} = zeros(length(f),size(U_all_j,3));  
                for p = 1:size(U_all_j,3)
                    FRFmp{j}(:,p) = Y_all_j(:,:,p)./U_all_j(:,:,p);
                end     
            end
            FRFs(:,i) = FRFs(:,i)/numel(data{i});
            FRFmp_all = [FRFmp{:}];
        end 
        
        for i = 1:ndata
           if np(i) > 1
               for p = 1:np(i)
                   noiseVar(:,i) = noiseVar(:,i) + (abs(FRFmp_all(:,p) - FRFs(:,i))).^2;
               end
           end
           noiseVar(:,i) = noiseVar(:,i)/(np(i)*(np(i) - 1));  
        end 
        noiseVar(:,~any(noiseVar,1)) = []; 
        
        if size(noiseVar,2)>= 1 % if there was at least one realization with more than one period
            blaNoiseVar = zeros(length(f),1);
            for i = 1:size(noiseVar,2)
                blaNoiseVar = blaNoiseVar + noiseVar(:,i);
            end
            blaNoiseVar = blaNoiseVar/(size(noiseVar,2)^2); 
            blaNoiseSTDMod = FRDmod(sqrt(blaNoiseVar), f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', [], 'OutputName', []);
        else
            blaNoiseVar = [];
            blaNoiseSTDMod = [];  
        end
        FRFbla = mean(FRFs, 2);
        if ndata > 1
            blaVar = zeros(length(f),1);
            for i = 1:ndata
                blaVar = blaVar + (abs(FRFs(:,i) - FRFbla)).^2;
            end
            blaVar = blaVar/(ndata*(ndata-1));
            blaSTDMod = FRDmod(sqrt(blaVar), f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', [], 'OutputName', []);
        else 
            blaVar = [];
            blaSTDMod = [];
        end   
        diag = struct('allFRFestimate', FRFs, 'sampleTotalVariance', blaVar, 'sampleTotalVarianceModel', blaSTDMod, 'sampleNoiseVariance', blaNoiseVar, 'sampleNoiseVarianceModel', blaNoiseSTDMod);
        model = FRDmod(FRFbla, f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', labels.input{1}, 'OutputName', labels.output{1});  
    
    case 'nonlinDetect'   
        f = data{1}{1}.Frequency; % the same for each period of each realization
        allFRFmods = cell(ndata,1);
        FRFs = zeros(length(f), ndata); % ndata is the number of realizations (mean in each realization)
        noiseVar = zeros(length(f), ndata); 
        np = zeros(ndata,1);
        
        for i = 1:ndata  
            np(i) = numel(data{i});
            allFRFmods{i} = cell(1,np(i));
            FRFmp = zeros(length(f),np(i));
            for p = 1:np(i)
                FRFmp(:,p) = squeeze(data{i}{p}.ResponseData);
                allFRFmods{i}{p} = FRDmod(FRFmp(:,p), f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', labels.input{1}, 'OutputName', labels.output{1});
            end
            FRFs(:,i) = mean(FRFmp,2);
            if np(i) > 1
                for p = 1:np(i)
                    noiseVar(:,i) = noiseVar(:,i) + (abs(FRFmp(:,p) - FRFs(:,i))).^2;
                end
                noiseVar(:,i) = noiseVar(:,i)/(np(i)*(np(i) - 1));
            end
        end 
        noiseVar(:,~any(noiseVar,1)) = []; 
        
        if size(noiseVar,2) >= 1 % if there is at least one realization with more than one period
            blaNoiseVar = zeros(length(f),1);
            for i = 1:size(noiseVar,2)
                blaNoiseVar = blaNoiseVar + noiseVar(:,i);
            end
            blaNoiseVar = blaNoiseVar/(size(noiseVar,2)^2); 
            blaNoiseSTDMod = FRDmod(sqrt(blaNoiseVar), f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', [], 'OutputName', []);
        else
            blaNoiseVar = [];
            blaNoiseSTDMod = [];  
        end
        
        FRFbla = mean(FRFs, 2);
        if ndata > 1
            blaVar = zeros(length(f),1);
            for i = 1:ndata
                blaVar = blaVar + (abs(FRFs(:,i) - FRFbla)).^2;
            end
            blaVar = blaVar/(ndata*(ndata-1));
            blaSTDMod = FRDmod(sqrt(blaVar), f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', [], 'OutputName', []);
        else 
            blaVar = [];
            blaSTDMod = [];
        end   
        diag = struct('sampleTotalVariance', blaVar, 'sampleTotalSTDModel', blaSTDMod, 'sampleNoiseVariance', blaNoiseVar, 'sampleNoiseSTDModel', blaNoiseSTDMod);
        model = FRDmod(FRFbla, f, data{1}{1}.Ts, 'FrequencyUnit', 'Hz', 'InputName', labels.input{1}, 'OutputName', labels.output{1});  
        
    case 'Robust_NL'
        
        assert(length(labels.input) == 1 && length(labels.output) == 1, 'Robust_NL only supports SISO systems.'); 
        
        % split up all measurement data in periods
        periods = cellfun(@(x) cellfun(@split, x, 'un', 0), data, 'un', 0);
        nops = cellfun(@(x) cellfun(@length, x, 'un', 0), periods, 'un', 0);
        nops = cell2mat([nops{:}]);
        
        if ~all(diff(nops)==0)
            warning('Some measurements have more periods than others. I will only consider the maximal number of periods that every measurement has in common.');
            assert(min(nops)>0, 'At least one period is required.');
            periods = cellfun(@(x) cellfun(@(y) y(end-min(nops)+1:end), x, 'un', 0), periods, 'un', 0);
        end
            
        % calculate the I/O Fourier coefficients
        spectra = cellfun(@(x) cellfun(@(y) cellfun(@(z) spectrum(z, [labels.input labels.output], 'periodic'), y, 'un', 0), x, 'un', 0), periods, 'un', 0);
        [~,f] = spectrum(periods{1}{1}{1},{1},'periodic');
        
        % cast into the right input structure
        Yall = zeros(length(spectra), min(nops), length(f));
        Uall = zeros(length(spectra), min(nops), length(f));
        
        for i=1:length(spectra)
            for j=1:min(nops)
                s = spectra{i}{1}{j};
                Yall(i,j,:) = s(:,2);
                Uall(i,j,:) = s(:,1);
            end
        end
     
        Rall = Uall; % we assume noise-free input
        
        % Robust_NL_Anal
        [G, Y, U, CYU] = Robust_NL_Anal(Yall, Uall, Rall);
        
        % parse output into the right toolbox structure
        model = IdentFRDmod(G.mean, f, 'FrequencyUnit', 'Hz'); % the VUB toolbox runs in Hz
        model.Y_ = Y;
        model.U_ = U;
        model.CYU_ = CYU;
        
        diag = []; allFRFmods = [];
        
    case 'RobustLocalPolyAnal'
        
        % keep the same number of periods for every measurement
        nops = cellfun(@(x) cellfun(@nop, x, 'un', 0), data, 'un', 0);
        nops = cell2mat([nops{:}]);
        
        if ~all(diff(nops)==0)
            warning('Some measurements have more periods than others. I will only consider the maximal number of periods that every measurement has in common.');
            assert(min(nops)>0, 'At least one period is required.');
        end
        data = cellfun(@(x) cellfun(@(y) clip(y,'lastnper',min(nops)), x, 'un', 0), data, 'un', 0);

        % cast into the right input structure
        inputdata.u = zeros(length(labels.input), length(labels.input), ndata, size(signal(data{1}{1},{1}),1));
        inputdata.y = zeros(length(labels.output), length(labels.output), ndata, size(signal(data{1}{1},{1}),1));
        
        for i=1:length(labels.input)
            for j=1:ndata
                inputdata.u(:,i,j,:) = data{j}{i}.signal(labels.input)';
            end
        end
        
         for i=1:length(labels.output)
            for j=1:ndata
                inputdata.y(:,i,j,:) = data{j}{i}.signal(labels.output)';
            end
        end
     
        inputdata.r = inputdata.u; % we assume noise-free input
        inputdata.N = length(data{1}{1});
        inputdata.Ts = data{1}{1}.Ts;
        inputdata.ExcitedHarm = round((data{1}{1}.excitation_.excf*data{1}{1}.Ts*length(data{1}{1})))';
        
        % RobustLocalPolyAnal
        method = struct(); % for now, no options for the user to specify
        [CZ, Z, freq, G, CvecG, dof, CL] = RobustLocalPolyAnal(inputdata, method);
        
        % parse output into the right toolbox structure
        model = IdentFRDmod(G, freq, 'FrequencyUnit', 'Hz'); % the VUB toolbox runs in Hz
        % other parameters ?? -> don't know what these things mean, to be discussed with Jan?
        
        diag = []; allFRFmods = [];
        
    otherwise
        error('The non parametric identification method you are looking for is not supported (yet?)')
end
        
end
    






