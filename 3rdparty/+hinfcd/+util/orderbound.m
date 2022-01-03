function [nb,b,i,lambda] = orderbound(P,ny,nu)
% ORDERBOUND Returns a guaranteed H-infinity controller order upperbound
%
% [nb,b,i] = ORDERBOUND(P, ny, nu) returns an upperbound nb on the 
% required controller order for the descriptor plant P with ny measured 
% outputs and nu control inputs. The output i is 0 if the (potential) 
% order reductions stem from unstable or infinite zeros in Pzu, and is 1 
% if they stem from unstable or infinite zeros in Pyw. The number of 
% potential order reductions is b. lambda is the unstable or infinite
% invariant zero corresponding to nb and b. This function implements the 
% characterization of X. Xin, S. Hara and M. Kaneda, "Reduced-Order Proper 
% H-infinity Controllers for Descriptor Systems: Existence Conditions and 
% LMI-Based Design Algorithms", IEEE Transactions on Automatic Control, 
% vol. 53, no. 5, June 2008. 

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

    % get the subplants
    Pzu = P(1:(end-ny),(end-nu+1):end);
    Pyw = P((end-ny+1):end,1:(end-nu)); 
    
    % bound based on Pzu
    izzu = tzero(Pzu);
    izzu = izzu(real(izzu)>=0); 
    nzu = rank([zeros(size(Pzu.E)) Pzu.E zeros(size(Pzu.B)); 
                Pzu.E              Pzu.A Pzu.B             ; 
                zeros(size(Pzu.C)) Pzu.C Pzu.D             ]); 
    lambdazu = Inf; 
    for i=1:length(izzu)
        cand = rank(Pzu.E) + rank([Pzu.A-Pzu.E*izzu(i)  Pzu.B ; 
                                   Pzu.C                Pzu.D]);
        if nzu > cand
            nzu = cand;
            lambdazu = izzu(i);
        end
    end
    nzu = nzu - rank([Pzu.E zeros(size(Pzu.B)); 
                      Pzu.A Pzu.B             ; 
                      Pzu.C Pzu.D             ]);
                  
    % bound based on Pyw
    izyw = tzero(Pyw);
    izyw = izyw(real(izyw)>=0); 
    nyw = rank([zeros(size(Pyw.E))  Pyw.E zeros(size(Pyw.B)); 
                Pyw.E               Pyw.A Pyw.B             ; 
                zeros(size(Pyw.C))  Pyw.C Pyw.D             ]);
    lambdayw = Inf; 
    for i=1:length(izyw)
        cand = rank(Pyw.E) + rank([Pyw.A-Pyw.E*izyw(i)   Pyw.B ; 
                                   Pyw.C                 Pyw.D]);
        if nyw > cand
            nyw = cand;
            lambdayw = izyw(i);
        end
    end
    nyw = nyw - rank([Pyw.E                Pyw.A   Pyw.B;
                      zeros(size(Pyw.C))   Pyw.C   Pyw.D]);
    
    % return the smallest upperbound
    if(nzu <= nyw)
        nb = nzu;
        i = 0;
        lambda = lambdazu;
    else
        nb = nyw;
        i = 1;
        lambda = lambdayw; 
    end
    b = size(P.A,1)-nb;
end