function [obj,g,T,Ti] = balreal(obj,opts,zerotol)
% BALREAL Gramian balanced realization of the slow subsystem of a descriptor system
%
% Since balancing the realization of a descriptor system eventually ends up
% in a Kronecker-Weierstrass canonical form of which the slow subsystem is
% balanced and the subspaces spanned by the pencil of the fast subsystem
% are deflated, the Kronecker-Weierstrass canonical form is calculated
% first, after which the slow subsystem is balanced by the standard MATLAB
% command. See T. Stykel, "Analysis and Numerical Solution of Generalized
% Lyapunov Equations", PhD Thesis, Faculty of Mathematics and Natural 
% Sciences, TU Berlin, 2002, for more details.
%
% [b,g,T,Ti] = balreal(a) returns a Kronecker-Weierstrass realization b of 
% a of which the stable portion of the slow subsystem has a balanced
% realization, i.e. its controllability and observability Gramians g are
% equal and diagonal. T is the state transformation matrix from the slow
% subsystem of a to the slowsubsystem of b, Ti = inv(T). 
%
% [b,g,T,Ti] = balreal(a,opts) returns the same as the previous syntax, but
% additionally takes options as a hsvdOptions object, e.g. for tolerances
% for the stable/unstable decomposition, Gramian computation, etc. 
%
% [b,g,T,Ti] = balreal(a,opts,zerotol) returns the same as the previous
% syntax, but also allows to specify the tolerance zerotol (default: 
% sqrt(eps)) for small numbers to be considered as exact zeros while 
% calculating the Kronecker-Weierstrass canonical form.
%
% See the BALREAL documentation for more details on the functionality of
% balreal. 
%
% See also BALREAL, HINFCD.KRONREAL. 

% This file is part of hinfcd.
% Copyright (c) 2019, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% hinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% hinfcd is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with hinfcd. If not, see <https://www.gnu.org/licenses/>.

    
    switch nargin
        case 1
            opts = hsvdOptions;
            zerotol = sqrt(eps); 
        case 2
            if isempty(opts) || ~isa(opts,'hsvd')
                opts = hsvdOptions;
            end
            zerotol = sqrt(eps); 
    end
    
    % 1. Kronecker-Weierstrass decomposition
    [obj,~,~,nslow] = kronreal(obj,zerotol);
    Gslow = hinfcd.dss(obj.A(1:nslow,1:nslow),obj.B(1:nslow,:),obj.C(:,1:nslow),obj.D,obj.E(1:nslow,1:nslow),obj.Ts);
    Gfast = hinfcd.dss(obj.A(nslow+1:end,nslow+1:end),obj.B(nslow+1:end,:),obj.C(:,nslow+1:end),zeros(size(obj.D)),obj.E(nslow+1:end,nslow+1:end),obj.Ts);
    
    % 2. Balanced realization
    [slow,g,T,Ti] = balreal@DynamicSystem(Gslow,opts);
    if ~isprop(slow,'E') || isempty(slow.E); slow.E = eye(size(slow.A)); end
    obj.set('A',blkdiag(slow.A,Gfast.A),'B',[slow.B ; Gfast.B],'C',[slow.C Gfast.C],'D',slow.D,'E',blkdiag(slow.E,Gfast.E));

end