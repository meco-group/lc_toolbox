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

function [Orec, Crec] = ss_get_OC(sys, varargin)
%SS_GET_OC Get observability and controlability matrix of LTI SS system
%   Detailed explanation goes here

if(~isempty(varargin))
    ip = varargin{1};
else
    if(isempty(sys.e)||(rank(sys.e,1e-6)==size(sys.a,1)))
        ip = 'p';
    else
        ip = 'i';
    end
end

switch ip
    case 'p'
        Orec = obsv(sys);
        Crec = 
        
    case 'i'
        P = sys;
        Orec = zeros(size(P.A));
        Oeig = zeros(size(P.A));
        [Vae,Dae] = eig(P.A\P.E);
        Dae = diag(Dae);
        Vaeinv = pinv(Vae)%\eye(size(Vae));
        for k=0:(size(P.A,1)-1)
            Orec(k+1,:) = P.C*((P.A\P.E)^k);
            Oeig(k+1,:) = P.C*(Vaeinv*(diag(Dae.^k))*Vae);
        end
%         rorec = rank(Orec,1e-6)
%         roeig = rank(Oeig,1e-6)

        Crec = zeros(size(P.A));
        Ceig = zeros(size(P.A));
        Ainv = P.a\eye(size(P.a));
        [Va,Da] = eig(Ainv);
        Da = diag(Da);
        Vainv = pinv(Va);%\eye(size(Va));
        for k=0:(size(P.A,1)-1)
            Crec(:,k+1) = (Ainv^k)*P.B; 
            Ceig(:,k+1) = Vainv*diag(Da.^(k))*Va*(P.B);
        end
%         rcrec = rank(Crec,1e-6)
%         rceig = rank(Ceig,1e-6)

%         r=min(ro,rc)
        
        
end



end

