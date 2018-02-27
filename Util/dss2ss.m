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

function [sysout,is_pss,r] = dss2ss(sysin)
%DSS2SS Converts a descriptor state space model into a proper state space
%model.
%   Converts a descriptor state space model into a proper state space
%   model. States that are unnecessary are removed. This does not include a
%   possible cancellation of poles and zeros.
%   In case the underlying model is improper, the the output state
%   space model is still in the descriptor form.
%   input:
%   ipss = (im)proper state space model to be converted
%   output:
%   pss = (im)proper output state space model
%   is_pss = 0 if the input system is improper, 1 if the input is proper and
%   the output therefore is in a regular state space form
%   r = amount of removed states

flag = false;
if isa(sysin,'AbstractDSSmod')
    ipss = std(sysin);
    flag = true;
else
    ipss = sysin;
end

aug_state = [];

Nstate = size(ipss.E,2);
Nstate_ex = 0;
Nin = size(ipss.B,2);
Nin_ex = Nin;
Nout = size(ipss.C,1);
r = 0;

% Diagonalize E
[Ad,Bd,Cd,Ed,~,~] = diagE_ip(ipss.A,ipss.B,ipss.C,ipss.E);
Dd = ipss.D;

is_pss = all(diag(Ed)==1);
if(~is_pss)
    M = [eye(Nout+Nstate,Nout) [zeros(Nout,Nstate);Ed] (-[Cd Dd;Ad Bd])];

    % Do staircase algorithm
    done = 0;
    while(~done)
        I = find(all(M(1:(Nout+Nstate+Nstate_ex),1:(Nout+Nstate))==0,2)); % equation numbers
        if(~isempty(I))
            % Find states that might require augmented dynamics
            aug_state_prop = zeros(size(I));
            for(k=1:length(I))
                t = find(M(I(k),Nout+Nstate+(1:Nstate))==1,1);
                if(~isempty(t))
                    aug_state_prop(k) = t;
                end
            end
            [~,i] = setdiff(aug_state_prop,aug_state); aug_state = aug_state_prop; % Check which states are new in the list

            % Augment the dynamics or stop augmenting
            if(~isempty(i))
                Nin_ex_flag = 0;
                for(k=1:length(i))
                    Nstate_ex = Nstate_ex+1;
                    M(Nout+Nstate+Nstate_ex,Nout+(1:Nstate)) = M(I(i(k)),Nout+Nstate+(1:Nstate));
                    if(~all(M(I(i(k)),Nout+2*Nstate+(1:Nin_ex))==0))
                        Nin_ex_flag = 1;
                        M(Nout+Nstate+Nstate_ex,Nout+2*Nstate+Nin+(1:Nin_ex)) = M(I(i(k)),Nout+2*Nstate+(1:Nin_ex));
                    end
                end
                if(Nin_ex_flag)
                    Nin_ex = Nin_ex + Nin;
                end
                M = rref(M);
                k = size(M,1);
                while(all(M(k,:)==0))
                    M(k,:) = [];
                    Nstate_ex = Nstate_ex-1;
                    k = k-1;
                end     
            else
                done = 1;
            end
        else
            done = 1;
        end
    end

    %% Clean up the system
    M(all(M(1:(Nstate+Nout),1:(Nstate+Nout))==0,2),:) = 0; % Check if there are unresolved state equations and ignore them  
    state_rm = find(all(M(1:(Nstate+Nout),Nout+Nstate+(1:Nstate))==0,1)); 
    r = length(state_rm);
    M([(Nout+state_rm) (Nout+Nstate+(1:Nstate_ex))],:) = [];
    M(:,[(Nout+state_rm) (Nout+Nstate+state_rm)]) = [];

    Nstate = Nstate - r;

    A = -M(Nout+(1:Nstate),Nout+Nstate+(1:Nstate)); B = -M(Nout+(1:Nstate),Nout+2*Nstate+(1:Nin_ex));
    C = -M(1:Nout,Nout+Nstate+(1:Nstate)); D = -M(1:Nout,Nout+2*Nstate+(1:Nin_ex));
    E = M(Nout+(1:Nstate),Nout+(1:Nstate));

    Bt = zeros(size(B)); Bt(:,(Nin_ex-Nin+1):Nin_ex) = B(:,(Nin_ex-Nin+1):Nin_ex);
%     Dt = zeros(size(D)); Dt(:,1:Nin) = D;
    Dt = D;
    for(k=(Nin_ex/Nin-2):-1:0)
        Bt(:,((k)*Nin+1):((k+1)*Nin)) = B(:,((k)*Nin+1):((k+1)*Nin)) + A*Bt(:,((k+1)*Nin+1):((k+2)*Nin));
        Dt(:,((k)*Nin+1):((k+1)*Nin)) = D(:,((k)*Nin+1):((k+1)*Nin)) + C*Bt(:,((k+1)*Nin+1):((k+2)*Nin));
    end

    Dt_zero = Dt(:,(Nin+1):end);
    is_pss = all(Dt_zero(:)==0) & ~isempty(Dt_zero);

%     if(is_pss)
        Ad = A; Bd = Bt(:,1:Nin);
        Cd = C; Dd = Dt(:,1:Nin);
        Ed = eye(size(Ad));
%         pss = ss(At,Bt,Ct,Dt,ipss.Ts);
%     else
%         pss = ipss;
%     end
% else
%     pss = ss(Ad,Bd,Cd,Dd,ipss.Ts);
end

% Assign output
if flag
    sysout = sysin;
    sysout = sysout.setdssdata(Ad,Bd,Cd,Dd,Ed);
else
    sysout = dss(Ad,Bd,Cd,Dd,Ed,sysin);
end
end

