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

clear all
close all
clc

%% Declare systems and derive state space representation
F{1,1} = zpk([2],[1 -10 2 -6 9],1);
F{2,1} = zpk([3 2],[4 -4],1);
F{3,1} = zpk([1 -4 5 60],[],1);
F{4,1} = zpk([40 5],[3],1);

for k = 1:length(F)
    Fss{k,1} = ss(F{k,1});
    npz(k,1) = length(pole(F{k,1}));
    npz(k,2) = length(zero(F{k,1}));
end

%% Derive all combinations of the defined systems and compose combined systems description
mult = {};
N = length(F);
for(k = 1:N)
    C = combnk(1:N,k);
    for(j = 1:size(C,1))
        mult{length(mult)+1,1} = C(j,:);
    end
end

for(k = 1:length(mult))
    H{k,1} = F{mult{k,1}(1,1),1};
    Hss{k,1} = Fss{mult{k,1}(1,1),1};
    for(j = 2:length(mult{k,1}))
        H{k,1} = H{k,1}*F{mult{k,1}(1,j),1};
        Hss{k,1} = Hss{k,1}*Fss{mult{k,1}(1,j),1};
    end
    H{k,1} = ss(H{k,1});
end

%% Use dss2ss to derive a state space description if possible
for(k = 1:length(Hss))
    tic;
    [Gss{k,1},is_ss(k,1)] = dss2ss(Hss{k,1});
    time(k,1) = toc;
    is_ss_pz(k,1) = (sum(npz(mult{k,1},1)-npz(mult{k,1},2))>=0);
    if(is_ss_pz(k,1)==1)
        [ss_comp(k,1),~] = sscomp(Gss{k,1},H{k,1});
    else
        ss_comp(k,1) = 1; % don't check when the model is improper
    end
end

%% Check if the models were correctly converted from dss to ss
if(all(is_ss==is_ss_pz))
    disp('All systems are correctly identified as ss or dss: check passed');
else
    disp('Error when identifying systems: check not passed');
end

if(all(ss_comp==1))
    disp('All models describe the original systems correctly: check passed');
else
    disp('Error when comparing original and converted models: check not passed');
end


