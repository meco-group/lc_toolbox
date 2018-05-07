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

classdef AbstractLFTmod 
% AbstractLFTmod provides a blueprint consisting of the properties and
% methods every LFT model of LCToolbox should have. 
    
    properties
        M   % main system matrix
        Nu  % upper feedback matrix
        Nl  % lower feedback matrix
        E   % needed to represent improper systems
        x0  % initial state
        Ts = 0
    end
    
    methods (Abstract)
        mergeparameters(varargin);  % implemented by AbstractLPVmod
        simplify(self);             % implemented by LTILFTmod and LPVLFTmod
    end
    
    methods
        function self = AbstractLFTmod(M,Nu,Nl,E,Ts)
        % Constructor for AbstractLFTmod objects. 
        %
        % Parameters: 
        %  M : the system's M matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  Nu : the system's Nu matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  Nl : the system's Nl matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  E : the system's E matrix, either \c double or \c Function (in
        %  case of parameter dependency) 
        %  Ts : sampling time @type double
            if ~iscell(M)
                n = size(M) - size(Nu) - size(Nl) - size(E);
                M = mat2cell(M,[size(Nu,1),size(E,1),n(1),size(Nl,1)],[size(Nu,2),size(E,2),n(2),size(Nl,2)]);
            end
            self.M = M;
            self.Nu = Nu;
            self.Nl = Nl;
            self.E = E;
            self.x0 = zeros(size(M{2,2},1),1);
            if nargin > 4
                self.Ts = Ts;
            end
            check(self);
        end
        
        function self = setdata(self,M,Nu,Nl,E)
        % Sets the model's linear fractional transform matrices.
        %
        % Parameters: 
        %  M : the system's M matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  Nu : the system's Nu matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  Nl : the system's Nl matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  E : the system's E matrix, either \c double or \c Function (in
        %  case of parameter dependency) 
        %
        % Return values:
        %  self : returns the same object with new linear fraction
        %  transform matrices @type AbstractLFTmod
            self.M = M;
            self.Nu = Nu;
            self.Nl = Nl;
            self.E = E;
            check(self);
        end
        
        function [M,Nu,Nl,E] = getdata(self)
        % Returns the model's linear fractional transform matrices.
        %
        % Return values:
        %  M : the model's M matrix, either \c double or \c Function
        %  Nu : the model's Nu matrix, either \c double or \c Function
        %  Nl : the model's Nl matrix, either \c double or \c Function
        %  E : the model's E matrix, either \c double or \c Function
            M = self.M;
            Nu = self.Nu;
            Nl = self.Nl;
            E = self.E;
        end
        
        function [e] = e(self,varargin)
        % Returns the model's E matrix.
        %
        % @note In contrast with A(), B(), C() and D() for AbstractDSSmod
        % objects, E() does not exist, as E is already a property of
        % AbstractLFTmod. 
        %
        % Parameters:
        %  varargin : allows subscripting of the E matrix 
        %
        % Return values: 
        %  E : the model's E matrix, either \c double or \c Function
            e = self.E;
            if nargin>1, e = e(varargin{:}); end
        end
        
        %% Dimensions
        
        function nx_ = nx(self)
        % Returns the number of states of the model.
        %
        % Return values: 
        %  nx_ : the model's number of states @type double
            nx_ = size(self.M{2,2},1);
        end
        
        function nl_ = nlower(self)
        % Returns the size of the lower feedback matrix Nl. 
        %
        % Return values: 
        %  nl_ : the model's Nl matrix size @type double
            nl_ = size(self.Nl,1);
        end
        
        function nu_ = nupper(self)
        % Returns the size of the upper feedback matrix Nu. 
        %
        % Return values: 
        %  nu_ : the model's Nu matrix size @type double
            nu_ = size(self.Nu,1);
        end
        
        function nin_ = nin(self)
        % Returns the number of inputs of the model. 
        %
        % Return values: 
        %  nin_ : the model's number of inputs @type double
            siz_ = subsize(self);
            nin_ = siz_(2);
        end
        
        function nout_ = nout(self)
        % Returns the number of outputs of the model. 
        %
        % Return values: 
        %  nout_ : the model's number of outputs @type double
            siz_ = subsize(self);
            nout_ = siz_(1);
        end
        
        function siz = subsize(self)
        % Returns the size of the feedthrough matrix M33 of the model.
        %
        % Return values: 
        %  siz : the model's M33 matrix size @type double
            siz = size(self.M{3,3});
        end
        
        % Tools
        function dt_ = isdt(self)
        % Checks wheter the model is a discrete-time model.
        %
        % Return values: 
        %  dt_ : boolean reflecting whether the model is a discrete-time model @type logical
            dt_ = any(self.Ts>0);
        end
        
        function ct_ = isct(self)
        % Checks wheter the model is a continuous-time model.
        %
        % Return values: 
        %  ct_ : boolean reflecting whether the model is a continuous-time model @type logical
            ct_ = ~isdt(self);
        end
        
        function S = lft2ss(self)
        % Tries to get rid of the upper and lower feedback matrices (if possible). 
        %
        % Return values: 
        %  S : same model with only M22, M23, M32 and M33 nonzero (= a
        %  descriptor state-space model) @AbstractLFTmod
            assert(isnumeric(self.Nu) && isnumeric(self.Nl),'Cannot convert the LFT because Nu and/or Nl are not constant');
            assert(isnumeric(self.M{1,1}),'Cannot convert the LFT to a SS because M11 is not constant');
            assert(isnumeric(self.M{4,4}),'Cannot convert the LFT to a SS because M44 is not constant');
            S = AbstractLFTmod.clft({self.Nu},AbstractLFTmod.clft(self.M,{self.Nl}));
        end
        
        function self = transpose(self)
        % Transposes the model.
        %
        % Return values: 
        %  self : transpose of the model @type AbstractLFTmod
            M_ = transpose(reshape(cellfun(transpose,self.M(:),'un',0),[4,4]));
            Nu_ = transpose(self.Nl);
            Nl_ = transpose(self.Nu);
            E_ = transpose(self.E);
            self = LFTmod(M_,Nu_,Nl_,E_,self.Ts);
        end
        
%% Connections

        
%         function self = c2d(self,Ts,varargin)
%             assert(nargin>1, 'c2d requires at least two input arguments.');
%             assert(self.Ts ~= 0, 'Your system is already discrete!');
%             if nargin==2
%                 assert(isnumeric(Ts) && isscalar(Ts) && Ts>0, 'The sampling time should be a positive number.');
%                 method = 'tustin';
%             else
%                 method = varargin{1};
%                 assert(any(strcmp(varargin{1},{'fweuler','tustin'})),'Only ''fweuler'' and ''tustin'' are supported for the time being.');
%             end
%             
%             switch method
%                 case 'tustin'   
%                     % See Doyle, Packard and Zhou (1991): 'Review of LFTs, LMIs and ?',
%                     % 30th Conference on Decision and Control (CDC), Brighton, England.
% 
%                     % discretizator
%                     I = eye(self.nx);
%                     M_ = SSmod([I sqrt(2)*I ; Ts/sqrt(2)*I Ts/2*I],Ts);
%                     Nl_ = SSmod(zeros(self.nx,self.nx), self.E, eye(self.nx,self.nx), zeros(self.nx,self.nx), Ts);
%                     discretizator = lft(Nl_,M_);
%                     
%                     % remove the continuous time dynamics
%                     self.E = zeros(0,0);
%                     self.Ts = Ts;
%                     self.M(3,:) = cellfun(@vertcat, self.M(2,:), self.M(3,:), 'un', 0);
%                     self.M(2,:) = cellfun(@(x) zeros(0,size(x,2)), self.M(2,:), 'un', 0);
%                     self.M(:,3) = cellfun(@horzcat, self.M(:,2), self.M(:,3), 'un', 0);
%                     self.M(:,2) = cellfun(@(x) zeros(size(x,1),0), self.M(:,2), 'un', 0);
%                     x0_old = self.x0;
%                     self.x0 = zeros(0,1);
%                     
%                     % add the discrete time dynamics
%                     self = lft(discretizator,self);
%                     
%                     % transform the initial state to discrete time
%                     % TO DO - difficulty: depends on u(0) and p(0), cf. ODEmod
%                     self = self.setx0(zeros(self.nx,1)); 
%                     if any(x0_old); warning('Initial states are lost, not yet implemented!'); end
%                     
%                 case 'fweuler'
%                    
%                     % discretize
%                     self.Ts = Ts;
%                     self.M(2,:) = cellfun(@(x) x*Ts, self.M(2,:), 'un', 0);
%                     self.M{2,2} = self.M{2,2}+eye(size(self.nx));
%                     
%                     % initial states remain the same
%                     
%                 % case 'bweuler'
%                 % TO DO - similar to Tustin 
%                 
%             end
%         end
        
        function self = lft(self,other,nu,ny)
        % Calculates the star interconnection of two AbstractLFTmod objects.
        % Overloads MATLAB's \c lft() for LCToolbox models. 
        % 
        % Parameters:
        %  self : first model @type AbstractLFTmod
        %  other : other model @type AbstractLFTmod
        %  nu : first \c nu outputs of \c other are connected to last \c nu
        %  inputs of \c self @type double
        %  ny : last \c ny outputs of \c self are connected to first \c ny
        %  inputs of \c other @type double
        %
        % Return values:
        %  self: the new linear fraction transformation model @type
        %  AbstractLFTmod
            if ~isempty(other) % return other when empty
                if nargin == 2
                    if all(size(other) <= size(self))
                        [nu,ny] = size(other);
                    else
                        [nu,ny] = size(self);
                    end
                end
                
                idxout = [1:(size(self,1)-ny),size(self,1)+((nu+1):size(other,1)),...
                          (size(self,1)-ny)+(1:ny),size(self,1)+(1:nu)];
                idxin = [1:(size(self,2)-nu),size(self,2)+((ny+1):size(other,2)),...
                          size(self,2)+(1:ny),size(self,2)-nu+(1:nu)];

                if isnumeric(self), self = SSmod(self); end
                if isnumeric(other), other = SSmod(other); end      

                self = blkdiag(self,other);
                
                % rearrange the M matrix
                self = submodel(self,idxout,idxin);
                nout = size(self.M{3,1},1); nin = size(self.M{1,3},2);
                self.M(4,:) = cellfun(@(x,y)vertcat(x(nout-(nu+ny)+1:nout,:),y),self.M(3,:),self.M(4,:),'un',0);
                self.M(3,:) = cellfun(@(x)x(1:nout-(nu+ny),:),self.M(3,:),'un',0);
                self.M(:,4) = cellfun(@(x,y)horzcat(x(:,nin-(nu+ny)+1:nin),y),self.M(:,3),self.M(:,4),'un',0);
                self.M(:,3) = cellfun(@(x)x(:,1:nin-(nu+ny)),self.M(:,3),'un',0);
                self.Nl = blkdiag(eye(nu+ny),self.Nl);
                check(self);
            end
        end

        function blk = blkdiag(varargin)
        % Calculates the block diagonal of the models that are provided (stacking inputs and outputs in the order that is provided). 
        % Overloading MATLAB's \c blkdiag() for LCToolbox models. 
        %
        % @note This implementation of \c blkdiag() is only applicable
        % AbstractLFTmod objects. In case other model types are present in
        % \c varargin, it will automatically refer to another appropriate implementation.
        % 
        % Parameters:
        %  varargin : list of models that you want to stack (\c Model)
        %
        % Return values:
        %  blk : the new model @type AbstractLFTmod

            if all(cellfun(@(x)isa(x,'AbstractLFTmod'),varargin))
                % 1. Combine the dynamical systems
                % Combine M's
                Split = cellfun(@(x)reshape(x.M,[16,1]),varargin,'un',0);
                Split = horzcat(Split{:});
                M_ = cell(16,1);
                for k = 1:16
                    row = Split(k,:);
                    M_{k} = blkdiag(row{:});
                end
                M_ = reshape(M_,[4,4]);

                % Combine Nl/Nu
                Nu_ = cellfun(@(x)x.Nu,varargin,'un',0); Nu_ = blkdiag(Nu_{:});
                Nl_ = cellfun(@(x)x.Nl,varargin,'un',0); Nl_ = blkdiag(Nl_{:});
                E_ = cellfun(@(x)x.E,varargin,'un',0); E_ = blkdiag(E_{:});

                % Check sample time
                dt = cellfun(@(x)isdt(x),varargin);
                if any(dt)
                    Ts_ = varargin{find(dt,1)}.Ts;
                else
                    Ts_ = 0;
                end

                % Check the parameters
                p = mergeparameters(varargin{:});
                if ~isempty(p)
                    blk = LFTmod(M_,Nu_,Nl_,E_,p,Ts_);
                else
                    blk = LFTmod(M_,Nu_,Nl_,E_,Ts_);
                end
                
                % Combine the initial states
                x0_ = reshape(cell2mat(cellfun(@(x) x.x0, varargin, 'un', 0)'),[sum(cell2mat(cellfun(@(x) x.nx, varargin, 'un', 0))),1]); 
                blk = blk.setx0(x0_);
                
            else
                blk = varargin{1};
                for i = 2:nargin
                    blk = blkdiag_model(blk,varargin{i});
                end
            end
        end
        
        function [Ms,Nus,Nls,Es] = grid_eval(self,grid,args)
        % Evaluates the model for a grid of parameter values.
        %
        % Parameters:
        %  grid : cell with grid{i} a vector (\c double) with the gridpoints for
        %  scheduling parameter i
        %  args : cell with args{i} the name (\c char) of scheduling parameter i @type cell
        %
        % Return values:
        %  Ms : array of M matrices with the appropriate size @type double
        %  Nus : array of Nu matrices with the appropriate size @type double
        %  Nlss : array of Nl matrices with the appropriate size @type double
        %  Es : array of E matrices with the appropriate size @type double
            nargs = length(args);
            handle = @(x)permute(safegrideval(x,grid,args),[nargs+(1:2),1:nargs]);
            Ms = reshape(cellfun(handle,self.M(:),'un',0),[4,4]);
            Ms = cell2mat(Ms);
            Nus = permute(safegrideval(self.Nu,grid,args),[nargs+(1:2),1:nargs]);
            Nls = permute(safegrideval(self.Nl,grid,args),[nargs+(1:2),1:nargs]);
            Es = permute(safegrideval(self.E,grid,args),[nargs+(1:2),1:nargs]);
        end
        
        function [Ms,Nus,Nls,Es] = eval(self,val,args)
        % Evaluates the model for one combination of scheduling parameters.
        %
        % Parameters:
        %  val : vector with val(i) the value of scheduling parameter i
        %  @type double
        %  args : cell with args{i} the name (\c char) of scheduling parameter i @type cell
        %
        % Return values:
        %  Ms : array of M matrices with the appropriate size @type double
        %  Nus : array of Nu matrices with the appropriate size @type double
        %  Nls : array of Nl matrices with the appropriate size @type double
        %  Es : array of E matrices with the appropriate size @type double
            if size(val,1) == length(val); val = val'; end
            grid = num2cell(val);
            [Ms,Nus,Nls,Es] = grid_eval(self,grid,args);
        end
        
        function product = mtimes(self,other)
        % Multiplies a AbstractLFTmod with a numeric matrix (input or output
        % scaling). 
        %
        % Either \c self or \other should be a numerical matrix of
        % appropriate size. If two LFT models are multiplied,
        % this goes through Model::mtimes. 
        %
        % Parameters:
        %  self : first factor, either \c double or \c AbstractLFTmod
        %  other : second factor, either \c double or \c AbstractLFTmod
        %
        % Return values:
        %  product : scaled model @type AbstractLFTmod
            if isnumeric(self)
                product = other;
                product.M(3,:) = cellfun(@(x) self*x, product.M(3,:),'un',0);
            elseif isnumeric(other)
                product = self;
                product.M(:,3) = cellfun(@(x) x*other, product.M(:,3),'un',0);
            end
        end
        
        function self = setx0(self,x0_)
            assert(length(x0_) == self.nx() && any(size(x0_)==1), 'The number of initial states does not match the order of your system.');
            assert(length(x0_) == self.nx() && any(size(x0_)==1), ['The number of initial states (' num2str(length(x0_)) ') does not match the order of your system (' num2str(self.nx) ').']);
            self.x0 = x0_(:);
        end

        
    end
    
    methods(Access=protected)
        function self = submodel(self,idxout,idxin)
        % Only keeps the specified inputs and outputs of the model.
        %
        % Parameters: 
        % idxout : indices specifying which outputs to keep @type double
        % idxin : indices specifying which inputs to keep @type double
        %
        % Return values:
        %  self : the new model @type AbstractLFTmod
            if isempty(idxout), idxout = 1:nout(self); end
            if isempty(idxin), idxin = 1:nin(self); end
                    
            self.M(:,3) = transpose(cellfun(@(x)x(:,idxin),self.M(:,3),'un',0));
            self.M(3,:) = cellfun(@(x)x(idxout,:),self.M(3,:),'un',0);
        end
        
        function check(self)
        % Checks whether all dimensions of the linear fractional transform
        % matrices are consistent. 
            assert(all(subsize(self)>0),'One or more dimensions of the system is <0');

            % check dimensions of M
            S = cell2mat(cellfun(@size,self.M(:),'un',0));
            S1 = reshape(S(:,1),[4,4]);
            assert(all(all(repmat(S1(:,1),[1,3])==S1(:,2:4))),'Row dimensions mismatch');
            S2 = reshape(S(:,2),[4,4]);
            assert(all(all(repmat(S2(1,:),[3,1])==S2(2:4,:))),'Row dimensions mismatch');
            assert(all(size(self.M{1,1})==size(self.Nu)),'Nu dimensions mismatch');
            assert(all(size(self.M{4,4})==size(self.Nl)),'Nl dimensions mismatch');
        end
        
    end
    
    methods (Static,Access=private)
        function L = clft(M,N)
        % Calculates the star interconnection of an AbstractLFTmod object
        % and a numeric matrix of appropriate dimensions. Special case of
        % AbstractLFTmod::lft.
        % 
        % Parameters:
        %  N : constant feedback matrix @type double
        %
        % Return values:
        %  L : the new linear fraction transformation model @type
        %  AbstractLFTmod
            function m = cellmult(a,b)
                m = cell(size(a,1),size(b,2));
                for i = 1:size(a,1)
                    for j = 1:size(b,2)                        
                        t = cellfun(@mtimes,a(i,:),b(:,j),'un',0);
                        m(i,j) = t(1);
                        for k = 2:length(t)
                            m{i,j} = m{i,j} + t{k};
                        end
                    end
                end
            end
            
            function m = cellplus(a,b)
                m = reshape(cellfun(@plus,a(:),b(:),'un',0),size(a));
            end
            
            if all(size(M)==1)
                % Upper LFT - square M
                X = {M{1,1} / (eye(size(M{1,1}))-N{1,1}*M{1,1})};
                L = cellplus(N(2:end,2:end),cellmult(N(2:end,1),cellmult(X,N(1,2:end))));
            else
                % Lower LFT - square N
                X = {N{1,1} / (eye(size(N{1,1}))-M{end,end}*N{1,1})};
                L = cellplus(M(1:(end-1),1:(end-1)),cellmult(M(1:(end-1),end),cellmult(X,M(end,1:(end-1)))));
            end
        end
    end
    
end

