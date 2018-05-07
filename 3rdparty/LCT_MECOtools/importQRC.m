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

function [ out ] = importQRC( in )

    % Parse the file
    % ==============
    if ischar(in)
        rec = readlog(in);
    else
        rec = in;
    end
        
    e = in.excitation;
    if strcmp(e.type,'multisine')
        inputSignal = Multisine('label','multisine','fs',e.fs,'fwindow',[e.fmin,e.fmax],'ppp',e.period*e.fs);
    else
        inputSignal = Multisine('label','multisine','fs',e.fs,'fwindow',[e.fmin,e.fmax],'ppp',e.period*e.fs,'type','odd');
    end
    
    date = sprintf('%i-%i-%i_%i-%i',rec.time.year,rec.time.month,rec.time.day,rec.time.hour,rec.time.minute);
    
    out = TDMeasurementData('label',rec.file,'excitation',inputSignal,'data',rec.data,'datalabels',rec.labels,'date',date,'info',rec.comment,'periodic',true);
end
