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

classdef(InferiorClasses = {?zpk,?tf,?ss,?frd}) LPVLFTmod < AbstractLFTmod & AbstractLPVmod & Model
    %LFTmod Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function self = LPVLFTmod(M,Nu,Nl,E,parameters,varargin)
            self@AbstractLFTmod(M,Nu,Nl,E,varargin{:});
            self@AbstractLPVmod(parameters);
        end
        
        function sys = std(self)
            sys = std(gridme(self));
        end
        
        function sys = simplify(self)
            S = transpose(lft2ss(self));
            sys = LPVDSSmod(S{:},self.E,self.parameters(),self.Ts);
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
            [Ms,Nus,Nls,Es] = grid_eval(self,grid,args);
            
            % Make cellgrid of LTILFTmodels
            density = cellfun(@length,grid);
            d = num2cell(density);
            if length(d) == 1, d{2} = 1; end
            C = cell(d{:});
            for k = 1:prod(density)
                C{k} = LTILFTmod(Ms(:,:,k),Nus(:,:,k),Nls(:,:,k),Es(:,:,k),self.Ts);
            end
            
            mod = Gridmod(C,[args(:),grid(:)]);
        end
        
        function mod = evalme(self,val,args)
            function p = getargument(p)
                if isa(p,'SchedulingParameter')
                    p = p.tensor_basis.arguments;
                end
            end
            
            if ~iscell(args)
                args = {args};
            end
            args = cellfun(@(x) getargument(x),args,'un',0); % behind the scenes args are always strings
            [Ms,Nus,Nls,Es] = self.eval(val,args);
            mod = LTILFTmod(squeeze(Ms),squeeze(Nus),squeeze(Nls),squeeze(Es),self.Ts);
        end
        
        function product = mtimes(self,other)
            if isnumeric(self) || isnumeric(other)
                product = mtimes@AbstractLFTmod(self,other);
            else
                product = mtimes@Model(self,other);
            end
        end
        
        function dmod = c2d(self, Ts, method)
            
            % parse input
            switch nargin
                case 1
                    error('Not enough input arguments.');
                case 2
                    method = 'tustin';
                case 3
                    assert(any(strcmp(method,{'tustin','fweuler','bweuler'})), 'Unknown discretization method. Supported methods are: ''fweuler'', ''bweuler'' and ''tustin''.');
                otherwise
                    warning('c2d for LPV systems has the following syntax: c2d(mod,Ts) or c2d(mod,Ts,method). I''m ignoring your additional input arguments.');
            end
            assert(isscalar(Ts) && isreal(Ts) && isnumeric(Ts) && Ts>=0, 'Sampling time should be a positive real number.');
            if Ts==0; dmod = self; return; end
            
            if any(self.x0); warning('Initial states are lost while discretizing a parameter-dependent system.'); end
            
            % discretize the state-space matrices
            switch method
                case 'fweuler'
                    self.M(2,:) = cellfun(@(x) x*Ts, self.M(2,:), 'un', 0);
                    self.M{2,2} = self.M{2,2}+self.E;
                    fb = [];

                case 'bweuler'
                    self.M(2,:) = cellfun(@(x) x*Ts, self.M(2,:), 'un', 0);
                    ATs = self.M{2,2};
                    self.M{2,2} = self.E;
                    self.E = self.E-ATs;
                    fb = [];

                case 'tustin'
                % See Doyle, Packard and Zhou (1991): 'Review of LFTs, LMIs and mu', 
                % 30th Conference on Decision and Control (CDC), Brighton, England.
                    assert(rank(self.E)==length(self.E),'tustin''s method is currently only implemented for proper systems.'); 
                    I = eye(self.nx);
                    M_ = SSmod([I sqrt(2)*I ; sqrt(2)*I I], Ts);
                    Nl_ = SSmod(zeros(self.nx,self.nx), I, I, zeros(self.nx,self.nx), Ts);
                    fb = lft(Nl_,M_)*self.E*(Ts/2);
                    self.E = zeros(0,0); x0_old = self.x0; self.x0 = zeros(0,1); % to avoid problems in lft
                    self.M(3,:) = cellfun(@vertcat, self.M(2,:), self.M(3,:), 'un', 0);
                    self.M(2,:) = cellfun(@(x) zeros(0,size(x,2)), self.M(2,:), 'un', 0);
                    self.M(:,3) = cellfun(@horzcat, self.M(:,2), self.M(:,3), 'un', 0);
                    self.M(:,2) = cellfun(@(x) zeros(size(x,1),0), self.M(:,2), 'un', 0);
            end
            
            % return the discretized model
            self.Ts = Ts;
            dmod = lft(fb,self);
            
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
            s.M = cellfun(@(x) tostruct(x), self.M, 'un', 0);
            s.Nu = tostruct(self.Nu);
            s.Nl = tostruct(self.Nl);
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
            
            M = cellfun(@(x) tospline(x), s.M, 'un', 0); 
            Nu = tospline(s.Nu); 
            Nl = tospline(s.Nl); 
            E = tospline(s.E); 
            parameters = cellfun(@(x) SchedulingParameter(x.args{1},splines.Function(splines.TensorBasis({splines.BSplineBasis(x.bases{1}{1},x.bases{1}{2})},x.args),splines.Coefficient(x.coeffs)),x.rate_), s.parameters, 'un', 0);
            self = LPVLFTmod(M,Nu,Nl,E,parameters,s.Ts);
        end
        
    end
        
end