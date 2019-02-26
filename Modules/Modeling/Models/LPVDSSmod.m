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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) LPVDSSmod < AbstractDSSmod & AbstractLPVmod & Model
    %LPVDSSMOD Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function self = LPVDSSmod(A,B,C,D,E,parameters,varargin)
            self@AbstractDSSmod(A,B,C,D,E,varargin{:});
            self@AbstractLPVmod(parameters);
        end

        function sys = std(self)
            sys = std(gridme(self));
        end
        
        function sys = simplify(self)
            sys = self;
        end
        
        function mod = gridme(self,varargin)
            function chk = isvalidgrid(grid)
                isc = iscell(grid); 
                hastwocol = size(grid,2) == 2;
                if ~(hastwocol && isc); chk = false; return; end
                containsparam = all(cell2mat(cellfun(@ischar, grid(:,1), 'un', 0)));
                containsnumval = all(cell2mat(cellfun(@isnumeric, grid(:,2), 'un', 0)));
                chk = isc && containsparam && containsnumval;
            end
            
            if nargin>1 && isvalidgrid(varargin{1})
                args = varargin{1}(:,1)';
                grid = varargin{1}(:,2)';
            else
                [grid,args] = makegrid(self,varargin{:});
            end
            [As,Bs,Cs,Ds,Es] = grid_eval(self,grid,args);
            
            % Make cellgrid of LTIDSSmodels
            density = cellfun(@length,grid);
            d = num2cell(density);
            if length(d) == 1, d{2} = 1; end
            C = cell(d{:});
            for k = 1:prod(density)
                C{k} = LTIDSSmod(As(:,:,k),Bs(:,:,k),Cs(:,:,k),Ds(:,:,k),Es(:,:,k),self.Ts);
            end
            
            mod = Gridmod(C,[args(:),grid(:)]);
        end
        
        function mod = evalme(self,val,args)
            if ~iscell(args); args = {args}; end;
            if ~all(cellfun(@ischar,args))
                args = cellfun(@(x) x.tensor_basis.arguments,args,'un',0); % behind the scenes args are always strings
            end
            [As,Bs,Cs,Ds,Es] = self.eval(val,args);
            mod = LTIDSSmod(squeeze(As),squeeze(Bs),squeeze(Cs),squeeze(Ds),squeeze(Es),self.Ts);
        end
        
        function out = norm(self,varargin)
            sys.A = self.A; sys.B = self.B; sys.C = self.C; sys.D = self.D; sys.Ts = self.Ts;
            assert(rank(self.E)==self.nx,'E matrix does not have full rank, cannot proceed further')
            switch nargin
                case 1
                    options = [];
                case 2
                    options = struct('objective',varargin{1});
                case 3
                    options = struct('objective',varargin{1},'spec',varargin{2});
                otherwise
                    error('Invalid number of arguments, it should be norm(LPVDSSmod,objective,spec)')
            end
            [Primal,Dual] = LPV_analysis(sys,self.parameters,options);
            out = Primal.objective;
        end
        
        function product = mtimes(self,other)
            if isnumeric(self) || isnumeric(other)
                product = mtimes@AbstractDSSmod(self,other);
            else
                product = mtimes@Model(self,other);
            end
        end
        
        function s = saveobj(self)
            function [bases,args,coeffs] = getdata(spline)
                bases = spline.tensor_basis.bases;
                bases = cellfun(@(x) {x.knots, x.degree}, bases, 'un', 0);
                args = cellfun(@(x) spline.tensor_basis.argument(x), num2cell(1:length(bases)), 'un', 0);
                coeffs = spline.coeff.data;
            end
            
            function data = tostruct(matrix)
                if ~isnumeric(matrix)
                    [bases,args,coeffs] = getdata(matrix);
                    data = struct();
                    data.bases = bases; 
                    data.args = args;
                    data.coeffs = coeffs;
                else
                    data = matrix;
                end
            end
            
            s = struct(); 
            s.A = tostruct(self.A);
            s.B = tostruct(self.B);
            s.C = tostruct(self.C);
            s.D = tostruct(self.D);
            s.E = tostruct(self.E);
            s.Ts = self.Ts;
            s.parameters = cellfun(@(x) tostruct(x), self.parameters, 'un', 0); 
            for j=1:length(s.parameters)
                l = properties(self.parameters{j});
                for k=1:length(l)
                    s.parameters{j}.(l{k}) = self.parameters{j}.(l{k});
                end
            end
        end
    end
    
    methods(Static)
        function self = loadobj(s)
            function spline = setdata(bases,args,coeffs)
                bases = cellfun(@(x) splines.BSplineBasis(x{1},x{2}), bases, 'un', 0);
                b = splines.TensorBasis(bases,args);
                c = splines.Coefficient(coeffs);
                spline = splines.Function(b,c);
            end
            
            function sp = tospline(strct)
                if isstruct(strct)
                    sp = setdata(strct.bases,strct.args,strct.coeffs);
                else
                    sp = strct;
                end
            end
            
            A = tospline(s.A); 
            B = tospline(s.B); 
            C = tospline(s.C); 
            D = tospline(s.D); 
            E = tospline(s.E); 
            parameters = cellfun(@(x) SchedulingParameter(x.args{1},splines.Function(splines.TensorBasis({splines.BSplineBasis(x.bases{1}{1},x.bases{1}{2})},x.args),splines.Coefficient(x.coeffs)),x.rate_), s.parameters, 'un', 0);
            self = LPVDSSmod(A,B,C,D,E,parameters,s.Ts);
        end
        
    end
end