function Z = isectspace(varargin)
% ISECTSPACE Returns an orthonormal basis for the intersection of spaces
% Returns an orthonormal basis for the intersection of vector spaces. 
%
% This approach is adopted from G.H. Golub, C.F. Van Loan, "Matrix
% Computations", 4th edition, The John Hopkin University Press, Baltimore,
% 1996.
%
% Z = ISECTSPACE(A,B,C,...) returns Z such that the columns of Z are an
% orthonormal basis for the intersection of Col(A), Col(B), ...

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

    assert(nargin>=1, 'ISECTSPACE: Not enough input arguments.'); 
    n = size(varargin{1},1);
    assert(all(cellfun(@(x) isnumeric(x) && size(x,1)==n, varargin)), 'ISECTSPACE: Inputs must be matrices with the same number of rows.');
    tol = sqrt(eps); 
    
    if nargin==1    % trivial case
        Z = orth(varargin{1});
    elseif nargin>2 % recursion
        Z = isectspace(varargin{2:end});
        Z = isectspace(varargin{1},Z);
    else            % actual algorithm
        A = varargin{1}; 
        B = varargin{2};
        if size(A,2)<size(B,2); C = A; A = B; B = C; end
        [QA,~] = qr(A,0); [QB,~] = qr(B,0); 
        [UC,SC,~] = svd(QA'*QB,0); 
        Z = QA*UC(:,1:size(B,2));
        Z = orth(Z(:,1:sum(abs(diag(SC)-1)<tol))); 
    end