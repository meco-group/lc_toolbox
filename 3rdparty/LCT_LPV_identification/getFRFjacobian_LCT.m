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

function [ jacobianFRF ] = getFRFjacobian_LCT( As, Bs, Cs, Ds, Ts, freqStartLine, OmegaConc, p_local )

na = size(As,2);
nu = size(Bs,3);
ny = size(Cs,2);
nb = size(p_local,1);
nrofL = size(p_local,2);
inouts = ny*nu;
Jf = complex(zeros(ny*nu*numel(OmegaConc),nb*(na*na + na*nu + ny*na + ny*nu)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A
for j = 1:na*nb  
    for i = 1:na
        for k = 1:nrofL        
            for freqLine = freqStartLine(k):(freqStartLine(k+1)-1)    
                o = OmegaConc(freqLine); 
                if Ts > 0
                    z = exp(1i*o*Ts);  
                else
                    z = 1i*o; % continuous-time, for the local approach only
                end
                X = permute(Cs(k,:,:),[2 3 1])*((z*eye(na) - permute(As(k,:,:),[2 3 1]))^(-1));
                Y = ((z*eye(na) - permute(As(k,:,:),[2 3 1]))^(-1))*permute(Bs(k,:,:),[2 3 1]);
                Jf(((freqLine-1)*inouts+1):freqLine*inouts, (j-1)*na + i) = vec(X*Iij(na,na*nb,i,j)*kron(p_local(:,k),eye(na))*Y);   
            end   
        end   
    end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B
for j = 1:nu*nb  
    for i = 1:na
        for k = 1:nrofL      
            for freqLine = freqStartLine(k):(freqStartLine(k+1)-1)   
                o = OmegaConc(freqLine); 
                if Ts > 0
                    z = exp(1i*o*Ts);  
                else
                    z = 1i*o; % continuous-time, for the local approach only
                end   
                X = permute(Cs(k,:,:),[2 3 1])*((z*eye(na) - permute(As(k,:,:),[2 3 1]))^(-1));
                Jf(((freqLine-1)*inouts+1):freqLine*inouts, na*na*nb + ((j-1)*na + i)) = vec(X*Iij(na,nu*nb,i,j)*kron(p_local(:,k),eye(nu)));       
            end    
        end   
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C
for j = 1:na*nb  
    for i = 1:ny
        for k = 1:nrofL        
            for freqLine = freqStartLine(k):(freqStartLine(k+1)-1)     
                o = OmegaConc(freqLine); 
                if Ts > 0
                    z = exp(1i*o*Ts);  
                else
                    z = 1i*o; % continuous-time, for the local approach only
                end 
                Y = ((z*eye(na) - permute(As(k,:,:),[2 3 1]))^(-1))*permute(Bs(k,:,:),[2 3 1]);
                Jf(((freqLine-1)*inouts+1):freqLine*inouts,  na*(na+nu)*nb + ((j-1)*ny + i)) = vec(Iij(ny,na*nb,i,j)*(kron(p_local(:,k), eye(na)))*Y);      
            end  
        end     
    end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D
for j = 1:nu*nb  
    for i = 1:ny
        for k = 1:nrofL            
            for freqLine = freqStartLine(k):(freqStartLine(k+1)-1)
                Jf(((freqLine-1)*inouts+1):freqLine*inouts, na*(na+nu)*nb + ny*na*nb + ((j-1)*ny + i)) = vec(Iij(ny,nu*nb,i,j)*(kron(p_local(:,k), eye(nu))));
            end    
        end    
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jacobianFRF = [real(Jf); imag(Jf)];

end
