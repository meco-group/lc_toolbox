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

function smileMdl = BSplineSMILE(varargin)

meco_binaries('cpp_splines','develop')
import splines.*

required = {'griddedMdl','splineBasis','schParam'};
calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'griddedMdl'); assert(~exist('griddedMdl','var'),'You cannot define griddedMdl twice.'); griddedMdl = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'splineBasis'); assert(~exist('splineBasis','var'),'You cannot define griddedFRFs twice.'); splineBasis = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'schParam'); assert(~exist('schParam','var'),'You cannot define schParam twice.'); schParam = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},required); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end

%% Fetching local LTI models

sysd = cell(1,numel(griddedMdl.grid_));
for i = 1:numel(griddedMdl.grid_)
    %sysd{i} = balreal(std(griddedMdl.grid_{i}));
    sysd{i} = std(griddedMdl.grid_{i});
end

%% Computing observability/controllability matrices and selecting the reference model

O = zeros(size(sysd{1}.A,1)*size(sysd{1}.C,1),size(sysd{1}.A,1),numel(griddedMdl.grid_));
C = zeros(size(sysd{1}.A,1),size(sysd{1}.A,1)*size(sysd{1}.B,2),numel(griddedMdl.grid_));
for i = 1:numel(griddedMdl.grid_)  
    O(:,:,i) = obsv(sysd{i}.A,sysd{i}.C); 
    C(:,:,i) = ctrb(sysd{i}.A,sysd{i}.B);
end
[kappaOmax, kappaCmax] = deal(zeros(1,numel(griddedMdl.grid_)));
[To, Tc, T_o, T_c] = deal(zeros(size(sysd{1}.A,1),size(sysd{1}.A,2),numel(griddedMdl.grid_)));
for i = 1:numel(griddedMdl.grid_)
   [kappaO, kappaC] = deal(zeros(1,numel(griddedMdl.grid_)));
   for k = 1:numel(griddedMdl.grid_)   
      To(:,:,k) = (((O(:,:,i)')*O(:,:,i))^(-1))*(O(:,:,i)')* O(:,:,k);
      Tc(:,:,k) =  C(:,:,k)*(((C(:,:,i)')*C(:,:,i))^(-1))*(C(:,:,i)');
      kappaO(i) = cond(To(:,:,k));
      kappaC(i) = cond(Tc(:,:,k));
   end
   kappaOmax(i) = max(kappaO);
   kappaCmax(i) = max(kappaC); 
end
[minKappaO, iO] = min(kappaOmax);
[minKappaC, iC] = min(kappaCmax);

%% Similarity matrices for observability and controllability

[sys_co(numel(griddedMdl.grid_)), sys_cc(numel(griddedMdl.grid_))] = deal(struct('A',[],'B',[],'C',[],'D',[]));
 for i = 1:numel(griddedMdl.grid_)
    T_o(:,:,i) = ((((O(:,:,iO)')*O(:,:,iO))^(-1))* (O(:,:,iO)'))* O(:,:,i);
    sys_co(i).A = T_o(:,:,i)*sysd{i}.A*(T_o(:,:,i)^(-1));
    sys_co(i).B = T_o(:,:,i)*sysd{i}.B;
    sys_co(i).C = sysd{i}.C*(T_o(:,:,i)^(-1));
    sys_co(i).D = sysd{i}.D;   
    T_c(:,:,i) = C(:,:,iC)*((((C(:,:,i)')*C(:,:,i))^(-1))*(C(:,:,i)'));
    sys_cc(i).A = T_c(:,:,i)*sysd{i}.A*(T_c(:,:,i)^(-1));
    sys_cc(i).B = T_c(:,:,i)*sysd{i}.B;
    sys_cc(i).C = sysd{i}.C*(T_c(:,:,i)^(-1));
    sys_cc(i).D = sysd{i}.D;    
 end

%% Choosing the set easier to interpolate

if minKappaC < minKappaO
    sys_c = sys_cc;
    fprintf('Controllability based transformation \n')
else
    sys_c = sys_co;
    fprintf('Observability based transformation \n')
end

%% Fitting

Am = zeros(numel(griddedMdl.params_{2}),size(sys_c(1).A,1),size(sys_c(1).A,2));
Bm = zeros(numel(griddedMdl.params_{2}),size(sys_c(1).B,1),size(sys_c(1).B,2));
Cm = zeros(numel(griddedMdl.params_{2}),size(sys_c(1).C,1),size(sys_c(1).C,2));
Dm = zeros(numel(griddedMdl.params_{2}),size(sys_c(1).D,1),size(sys_c(1).D,3));
for i = 1:numel(griddedMdl.params_{2})
   Am(i,1:size(sys_c(i).A,1),1:size(sys_c(i).A,2)) = sys_c(i).A;
   Bm(i,1:size(sys_c(i).B,1),1:size(sys_c(i).B,2)) = sys_c(i).B;
   Cm(i,1:size(sys_c(i).C,1),1:size(sys_c(i).C,2)) = sys_c(i).C;
   Dm(i,1:size(sys_c(i).D,1),1:size(sys_c(i).D,2)) = sys_c(i).D;
end
AmBm = cat(3,Am,Bm);
CmDm = cat(3,Cm,Dm);
Tm = cat(2,AmBm,CmDm);

opti = OptiSpline();
if ~isa(splineBasis,'splines.TensorBasis')
    splineBasis = TensorBasis(splineBasis,argument(schParam));
end
A = opti.Function(splineBasis,[size(sys_c(1).A,1),size(sys_c(1).A,2)]);
B = opti.Function(splineBasis,[size(sys_c(1).B,1),size(sys_c(1).B,2)]);
C = opti.Function(splineBasis,[size(sys_c(1).C,1),size(sys_c(1).C,2)]);
D = opti.Function(splineBasis,[size(sys_c(1).D,1),size(sys_c(1).D,2)]);

T = [A B; C D];
T_eval = T.list_eval(griddedMdl.params_{2});
T_eval.dims()
e = T.list_eval(griddedMdl.params_{2}) - Tm;
en = matrix(e(:));
obj = (en)'*en;
opti.minimize(obj)
opti.solver('ipopt');
sol = opti.solve();
T = sol.value(T);

%% LPV model

As = slice(T,1:A.size(1),1:A.size(2));
Bs = slice(T,1:B.size(1),A.size(2)+1:A.size(2)+B.size(2));
Cs = slice(T,A.size(1)+1:A.size(1)+C.size(1),1:C.size(2));
Ds = slice(T,A.size(1)+1:A.size(1)+D.size(1),C.size(2)+1:C.size(2)+D.size(2));
smileMdl = SSmod(As,Bs,Cs,Ds,schParam,sysd{1}.Ts);

