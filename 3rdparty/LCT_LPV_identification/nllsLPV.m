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

function [nllsMdl] = nllsLPV( griddedFRFs, basisFcns, model, settings )

meco_binaries('cpp_splines','develop')
import splines.*
schParam = SchedulingParameter(griddedFRFs.params_{1}, [min(griddedFRFs.params_{2}), max(griddedFRFs.params_{2})]);
p_local = basisFcns.list_eval(griddedFRFs.params_{2})';

freqLines = 1;
OmegaConc = [];
FRFmConc  = [];  
As = zeros(numel(griddedFRFs.params_{2}),size(model.A.coeff_tensor,2),size(model.A.coeff_tensor,3));
Bs = zeros(numel(griddedFRFs.params_{2}),size(model.B.coeff_tensor,2),size(model.B.coeff_tensor,3));
Cs = zeros(numel(griddedFRFs.params_{2}),size(model.C.coeff_tensor,2),size(model.C.coeff_tensor,3));
Ds = zeros(numel(griddedFRFs.params_{2}),size(model.D.coeff_tensor,2),size(model.D.coeff_tensor,3));

for k = 1:numel(griddedFRFs.params_{2})
    modEval = model.evalme(griddedFRFs.params_{2}(k),{schParam});
    As(k,1:size(model.A.coeff_tensor,2),1:size(model.A.coeff_tensor,3)) = modEval.A;
    Bs(k,1:size(model.B.coeff_tensor,2),1:size(model.B.coeff_tensor,3)) = modEval.B;
    Cs(k,1:size(model.C.coeff_tensor,2),1:size(model.C.coeff_tensor,3)) = modEval.C;
    Ds(k,1:size(model.D.coeff_tensor,2),1:size(model.D.coeff_tensor,3)) = modEval.D;
    OmegaConc = [OmegaConc; 2*pi*griddedFRFs.grid_{k}.Frequency];
    FRFmConc = [FRFmConc; vec(griddedFRFs.grid_{k}.frdata)];
    freqLines = [freqLines; (freqLines(end) + numel(griddedFRFs.grid_{k}.Frequency))];
end

W_freq = diag((1/sqrt(sumsqr(abs(FRFmConc))))*ones(1,2*numel(FRFmConc)));

%% Generating mex files 

cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
if (~isempty(griddedFRFs.grid_))
    codegen -config cfg getFRFjacobian_LCT -args { zeros(size(p_local,2),size(model.A.coeff_tensor,2),size(model.A.coeff_tensor,3)), zeros(size(p_local,2),size(model.B.coeff_tensor,2),size(model.B.coeff_tensor,3)), zeros(size(p_local,2),size(model.C.coeff_tensor,2),size(model.C.coeff_tensor,3)), zeros(size(p_local,2),size(model.D.coeff_tensor,2),size(model.D.coeff_tensor,3)),0, zeros(numel(griddedFRFs.grid_)+1,1), zeros(numel(griddedFRFs.grid_)*numel(griddedFRFs.grid_{1}.frdata),1), zeros(size(p_local)) }                                                                                           
end

%% Levenberg-Marquardt

lambda = -1;
incFac = sqrt(10);
err_old = W_freq*getFRFerror_LCT(freqLines,OmegaConc,FRFmConc,griddedFRFs.params_{2},schParam,model);
J = W_freq*getFRFjacobian_LCT_mex(As,Bs,Cs,Ds,model.Ts,freqLines,OmegaConc,p_local);
K_old = err_old'*err_old;
warning('off','MATLAB:pack:InvalidInvocationLocation');
pack;

iteration = 0;
gradient = 2*J'*err_old;
gradNorm = norm(gradient);
disp([num2str(iteration) ' - cost: ' num2str(K_old) ' - gradNorm: ' num2str(gradNorm) ' - stable: ' num2str(isstable(model))]);
disp('Starting L.M. Optimization...')

while ((iteration < settings.MaxCount) && (gradNorm > settings.gradTol) && (lambda < 1e9))
    
    if lambda == -1  
        [Uj,Sj,Vj] = svd(J,'econ');
        lambda = sqrt(Sj(1,1)); 
        clear Uj Sj Vj
        pack   
    end
    if (isnan(K_old) || isinf(lambda)), break; end  
    [J0,scalingJ0] = normalise(J);
    nowexit = false;
    
    while not(nowexit) 
         
        q = - (J0'*J0 + (lambda^2)*eye(size(J0,2))) \ (J0'*err_old);
        dtheta = q./scalingJ0';
        
        dcoeff_tensor_A = permute(reshape(dtheta(1:numel(model.A.coeff_tensor)),size(model.A.coeff_tensor,2),size(model.A.coeff_tensor,3),size(model.A.coeff_tensor,1)),[3 1 2]);
        dcoeff_tensor_B = permute(reshape(dtheta((numel(model.A.coeff_tensor)+1):(numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor))),size(model.B.coeff_tensor,2),size(model.B.coeff_tensor,3),size(model.B.coeff_tensor,1)),[3 1 2]);
        dcoeff_tensor_C = permute(reshape(dtheta((numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+1):(numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor))),size(model.C.coeff_tensor,2),size(model.C.coeff_tensor,3),size(model.C.coeff_tensor,1)),[3 1 2]);
        dcoeff_tensor_D = permute(reshape(dtheta((numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor)+1):(numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor)+numel(model.D.coeff_tensor))),size(model.D.coeff_tensor,2),size(model.D.coeff_tensor,3),size(model.D.coeff_tensor,1)),[3 1 2]);
        
        A_new = splines.Function(model.A.tensor_basis,model.A.coeff_tensor + dcoeff_tensor_A);
        B_new = splines.Function(model.B.tensor_basis,model.B.coeff_tensor + dcoeff_tensor_B);
        C_new = splines.Function(model.C.tensor_basis,model.C.coeff_tensor + dcoeff_tensor_C);
        D_new = splines.Function(model.B.tensor_basis,model.D.coeff_tensor + dcoeff_tensor_D);
        
        model_new = SSmod(A_new,B_new,C_new,D_new,schParam,model.Ts);
        err = W_freq*getFRFerror_LCT(freqLines,OmegaConc,FRFmConc,griddedFRFs.params_{2},schParam,model_new);
        K = err'*err;
      
        if (K >= K_old) 
            lambda = lambda*incFac;
        else
            lambda = lambda/2;
            nowexit = true;
        end
        
    end
    
    iteration = iteration + 1;
    K_old = K;         
    err_old = err;
    model = model_new; 
    
    for k = 1:numel(griddedFRFs.params_{2})
        modEval = model.evalme(griddedFRFs.params_{2}(k),{schParam});
        As(k,1:size(model.A.coeff_tensor,2),1:size(model.A.coeff_tensor,3)) = modEval.A;
        Bs(k,1:size(model.B.coeff_tensor,2),1:size(model.B.coeff_tensor,3)) = modEval.B;
        Cs(k,1:size(model.C.coeff_tensor,2),1:size(model.C.coeff_tensor,3)) = modEval.C;
        Ds(k,1:size(model.D.coeff_tensor,2),1:size(model.D.coeff_tensor,3)) = modEval.D;
    end
    J = W_freq*getFRFjacobian_LCT_mex(As,Bs,Cs,Ds,model.Ts,freqLines,OmegaConc,p_local);
    gradient = 2*J'*err_old;
    gradNorm = norm(gradient);
    disp([num2str(iteration) ' - K: ' num2str(K_old) ' - grad: ' num2str(gradNorm) ' - lambda: ' num2str(lambda) ' - stable: ' num2str(isstable(model))]);     
    
end

disp('End of L.M. Optimization.')
nllsMdl = model;

