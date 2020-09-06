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

function [equal,error] = sscomp(ss1,ss2,varargin)
%SSCOMP Compares 2 realizations in state space
%   SSCOMP outputs whether or not 2 state space models describe the
%   same system. 
%   Errors may still occur for improper systems. TODO: fix this bug..
%   Inputs
%   ss1,ss2: 2 models to be compared
%   tol (optional): tolerance for the comparison and rank determination (default is 1e-9)
%   Outputs
%   equal: true if the models describe the same system
%   error: error between the 2 models, calculated as max(maxdiff(ssC*ssA*ssB),maxdiff(ssD))

if(nargin>2)
    tol = varargin{1};
else
    tol = 1e-6;
end

if ~isa(ss1,'Model')
    ss1 = fromstd(ss1);
end
if ~isa(ss2,'Model')
    ss2 = fromstd(ss2);
end

sz1 = size(ss1);
sz2 = size(ss2);

if all(sz1==sz2)
    if isa(ss1,'Gridmod') && isa(ss2,'Gridmod')
        if all(ss1.gridsize() == ss2.gridsize())
            equal = true;
            error = false;
            for k = 1:prod(ss1.gridsize())
                equal = equal && sscomp(ss1.grid_{k},ss2.grid_{k});
            end
        else
            equal = true;
            error = -1;
        end
    else
        % Make full descriptor systems
        if(isempty(ss1.E))
            ss1.E = eye(size(ss1.A));
        end
        if(isempty(ss2.E))
            ss2.E = eye(size(ss2.A));
        end

        if xor((rank(ss1.E,tol)~=size(ss1.E,1)),(rank(ss2.E,tol)~=size(ss2.E,1)))
            disp('1 of the systems is improper. Try to put it in a proper form using dss2ss and check again.')
            error = false;
            equal = false;
        else
            if((rank(ss1.E,tol)~=size(ss1.E,1))&&(rank(ss2.E,tol)~=size(ss2.E,1)))
                X1 = ss1.E; B1 = ss1.B; C1 = ss1.C; Y1 = ss1.A;
                X2 = ss2.E; B2 = ss2.B; C2 = ss2.C; Y2 = ss2.A;
            else
                X1 = ss1.A; B1 = ss1.B; C1 = ss1.C; Y1 = ss1.E;
                X2 = ss2.A; B2 = ss2.B; C2 = ss2.C; Y2 = ss2.E;
            end
            [X1,B1,C1,~,~,~] = diagE_ip(X1,B1,C1,Y1);
            [X2,B2,C2,~,~,~] = diagE_ip(X2,B2,C2,Y2);
            state_error = C1*X1*B1-C2*X2*B2;
            output_error = ss1.D-ss2.D;

            error = max(abs([state_error(:);output_error(:)]));
            equal = (error < tol);
        end
    end
else
    error = Inf;
    equal = false;
end
end