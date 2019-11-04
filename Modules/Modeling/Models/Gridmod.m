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

classdef (InferiorClasses = {?LPVLFTmod,?LTILFTmod,?LPVDSSmod,?LTIDSSmod,?FRDmod}) Gridmod < Model
% Create a grid of models. 
    
    properties
        grid_       % the actual grid of models
        params_     % the parameter values corresponding to the grid entries
    end
    
    methods
        
        
        function self = Gridmod(varargin)
        % Constructor for a grid of models.
        %
        % Parameters:
        %  varargin :  
        %   - for existing grid of models, use these as input arguments:
        %     - \c grid : cell representing the grid of models
        %     - \c params : cell, of which column 1 contains the name of
        %     the scheduling parameter (char) and column 2 contains a
        %     vector with parameter values corresponding to \c grid
        %   - for creating a default grid of, use this as input argument:
        %     - \c mod : LPVDSSmod or LPVLFTmod
        % 
        % Return values:
        %  self : the grid of models @type Gridmod
        
            if (~isa(varargin{1},'AbstractLPVmod')) || isa(varargin{1},'AbstractLTImod')
                grid = varargin{1};
                params = varargin{2};
                % construct a Gridmod object
                % <MAARTEN>: this is a very strange piece of code... what is
                % intended for? It repeats a model?
                % Answer: necessary to connect grid models to standard lti
                % models
                if isa(grid,'Model') && ~isa(grid,'Gridmod')
                    dims = cellfun(@length,params(:,2));
                    if length(dims) == 1, dims = [dims,1]; end
                    grid = repmat({grid},dims(:)');
                end
                % </MAARTEN>

                assert(isa(grid,'cell'),'The grid should be a cell.');
                assert(~isempty(grid),'The grid should at least contain one model.');
                assert(range(reshape(cellfun(@(x) size(x,1),grid),[numel(grid) 1]))==0 && range(reshape(cellfun(@(x) size(x,2),grid),[numel(grid) 1]))==0,'All systems of your grid should have the same dimensions.');
                assert(isa(params,'cell') && ((size(params,1) == ndims(grid)) || (size(params,1) == 1)),'The parameter space should be represented by a cell with as many rows as parameters.');
                assert(all(cellfun(@(x) isnumeric(x) && (size(x,1) == 1 || size(x,2) == 1),params(:,2))),'Parameter values should be stored in numeric vectors.');
                assert(all(cellfun(@(x) ischar(x),params(:,1))),'Parameter labels should be strings.');

                isfrd = cellfun(@(x) isa(x,'FRDmod'),grid);
                islti = cellfun(@(x) isa(x,'AbstractLTImod'),grid);
                assert(all(isfrd(:)) || all(islti(:)),'Gridded models should contain FRD or LTI models.');

                self.grid_ = grid;
                self.params_ = params;
            else
                self = gridme(varargin{:});
            end
        end
        
        function self = subs(self,params)
        % Eliminates one or multiple dimensions of the grid (only retain one
        % value for a parameter).
        %
        % Parameters:
        %  params: \c cell with names of the scheduling parameters in its first column and the sampled values in the second column @type cell
        % 
        % Return values:
        %  self : the new grid of models @type Gridmod
        
            assert(all(ismember(fieldnames(params),self.params_(:,1))),'At least one of the parameters you provided was not in the grid.');
            assert(all(structfun(@(x) isnumeric(x) && isscalar(x),params)),'You can only substitute parameters by numeric scalars.');
            names = fieldnames(params);
            
            cmd = cell(nparameters(self),1);
            cmd(:) = {':'};
            for i = 1:length(names)
                this_param = params.(names{i});
                par_number = find(ismember(self.params_(:,1),names{i}));
                assert(ismembertol(this_param,self.params_{par_number,2},1e-8),['The value ' num2str(this_param) ' is not included in the grid values of ' names{i} '.']);
                cmd(par_number) = {find(ismember(this_param,self.params_{par_number,2}))};
            end
            
            substr = substruct('()',cmd);
            self.grid_ = subsref(self.grid_,substr);
            self.params_(ismember(self.params_(:,1),names),:) = [];
        end
        
        function inp = nin(self)
        % Returns the number of inputs of the models in the grid.
        % 
        % Return values:
        %  inp : number of inputs @type double
            inp = size(subsref(subsref(self.grid_,substruct('()',num2cell(ones(nparameters(self),1)))),substruct('{}',{1})),2);
        end
        
        function outp = nout(self)
        % Returns the number of outputs of the models in the grid.
        % 
        % Return values:
        %  outp : number of outputs @type double
            outp = size(subsref(subsref(self.grid_,substruct('()',num2cell(ones(nparameters(self),1)))),substruct('{}',{1})),1);
        end
        
        function s = gridsize(self)
        % Returns the grid dimensions.
        % 
        % Return values:
        %  s : size of the grid @type double
            s = size(self.grid_);
        end
        
        function n = nparameters(self)
        % Returns the number of parameters that are sampled in the grid.
        % 
        % Return values:
        %  n : number of parameters @type double
            n = size(self.params_,1);
        end
        
        function self = lft(self,other,nu,ny)
        % Calculates the star interconnection of each model in the grid
        % with another model. 
        % Overloading MATLAB's \c lft() for LCToolbox models. 
        %
        % Either \c self or \c other is supposed to be a Gridmod object. 
        % Two Gridmod objects cannot be combined. 
        % 
        % Parameters:
        %  self : any Model object @type Model
        %  other : any Model object @type Model
        %  nu : first \c nu outputs of \c other are connected to last \c nu
        %  inputs of \c self @type double
        %  ny : last \c ny outputs of \c self are connected to first \c ny
        %  inputs of \c other @type double
        %
        % Return values:
        %  self: the new grid of models @type Gridmod
            isgridmod = cellfun(@(x) isa(x,'Gridmod'),{self,other});
            assert(sum(isgridmod) == 1,'Combining different grid of modelss is not supported (yet).');
            
            if isgridmod(1)
                grid = gridfun(self,@(x)lft(x,other,nu,ny));
            elseif isgridmod(2)
                grid = gridfun(other,@(x)lft(self,x,nu,ny));
            end
            self = Gridmod(grid,self.params_);
        end
        
        function self = blkdiag(varargin)
        % Calculates the block diagonal of every model in the grid with other models (stacking inputs and outputs in the order that is provided). 
        % Overloading MATLAB's \c blkdiag() for LCToolbox models. 
        %
        % Either \c self or \c other is supposed to be a Gridmod object.
        % Multiple Gridmod objects cannot be combined, unless their parameters and grids are equal. 
        % 
        % Parameters:
        %  varargin: list of Model objects, one of them being a Gridmod
        %  object. 
        %
        % Return values:
        %  self: the new grid of models @type Gridmod
            isgridmod = cellfun(@(x) isa(x,'Gridmod'),varargin);
            gridmods = varargin(isgridmod);
            othermods = varargin(~isgridmod);
            if sum(isgridmod)>1
                % combine the gridmods
                for k = 2:length(gridmods)
                    assert(~any(size(gridmods{1}.params_(:))-size(gridmods{2}.params_(:))), 'Multiple Gridmod objects can only be combined if their parameter grids are equal.');
                    assert(all(cellfun(@(x,y)all(eq(x,y)),gridmods{1}.params_(:),gridmods{k}.params_(:))),'Multiple Gridmod objects can only be combined if their parameter grids are equal.');
                end
                grids = cellfun(@(x)x.grid_,gridmods(:),'un',0);
                grid = cellfun(@blkdiag,grids{:},'un',0);
                combinedgridmods = Gridmod(grid,varargin{1}.params_);
                
                % combine the new gridmod with the other models 
                % ! change inputs and outputs to retain original order !
                if isempty(othermods)
                    self = combinedgridmods;
                else
                    self = blkdiag(combinedgridmods,othermods{:});
                    i = 0; no = sum(cellfun(@nout, gridmods)); ni = sum(cellfun(@nin, gridmods));
                    idxoutold{1} = 1:nout(varargin{1});
                    idxinold{1} = 1:nin(varargin{1});
                    for k=2:nargin
                        idxoutold{k} = max(idxoutold{k-1})+(1:nout(varargin{k}));
                        idxinold{k} = max(idxinold{k-1})+(1:nin(varargin{k}));
                    end
                    idxout = []; idxin = [];
                    for k=1:nargin
                        if ~isgridmod(k)
                            idxout = [idxout (idxoutold{k}+no-sum(cellfun(@nout, gridmods(1:i))))];
                            idxin = [idxin (idxinold{k}+ni-sum(cellfun(@nin, gridmods(1:i))))];
                        else
                            i = i+1;
                            idxout = [idxout (1:nout(gridmods{i}))+sum(cellfun(@nout, gridmods(1:i-1)))];
                            idxin = [idxin (1:nin(gridmods{i}))+sum(cellfun(@nin, gridmods(1:i-1)))];
                        end
                    end
                    self = submodel(self,idxout,idxin);
                end
            else
                % check which argument is the gridmod and combine
                idx = find(isgridmod);
                varargin(~isgridmod) = cellfun(@(x)Gridmod(x,varargin{idx(1)}.params_),varargin(~isgridmod),'un',0);
                self = blkdiag(varargin{:});
            end
        end
        
        function stdmodel = std(self)
        % Converts the LCToolbox model to a standard MATLAB array of
        % models.
        %
        % Return values:
        %  stdmodel : a standard MATLAB array of linear time invariant
        %  models. @type numlti 
            stdgrid = gridfun(self,@std);
            stdmodel = cat(3,stdgrid{:});
            stdmodel = reshape(stdmodel,size(stdgrid));
            
            outlist = cell(nparameters(self),1);
            [outlist{:}] = ndgrid(self.params_{:,2});
            t = [self.params_(:,1)'; outlist(:)'];
            stdmodel.SamplingGrid = struct(t{:});
        end
        
        function self = simplify(self)
        % Applies the \c simplify operation to every model that is a part
        % of the grid. 
        %
        % Return values:
        %  self : the new grid of models @type Gridmod
            self.grid_ = gridfun(self,@simplify);
        end
        
        function n = norm(self,varargin)
        % Calculates the norm of every model that is a part of the grid. 
        % Overloading MATLAB's \c norm() for LCToolbox models. 
        %
        % Return values:
        %  self : numerical array with the norms of every model in the grid
        %  @type double
            n = gridfun(self,@(x)norm(x,varargin{:}));
        end
        
        function el = element(self,varargin)
        % Returns a specific model of the grid. 
        % 
        % Parameters:
        %  varargin: indices of the model to be extracted
        %
        % Return values:
        %  self : the new grid of models @type Gridmod
            assert(length(varargin)==nparameters(self));
            el = self.grid_{varargin{:}};
        end
    end
    
    methods (Access=protected)
        function self = submodel(self,idxout,idxin)
        % Only keeps the specified inputs and outputs of every model of the
        % grid. 
        %
        % Parameters: 
        % idxout : indices specifying which outputs to keep @type double
        % idxin : indices specifying which inputs to keep @type double
        %
        % Return values:
        %  self : the new grid of models @type Gridmod
            self.grid_ = gridfun(self,@(x)submodel(x,idxout,idxin));
        end
    end
    
    methods (Access=private)
        function varargout = gridfun(self,f)
        % Applies the function \c f to every model in the grid. 
        %
        % Parameters:
        %  f : the function to be applied @type function_handle
        %
        % Return values:
        %  self : the new grid of models @type Gridmod
            [varargout{1:nargout}] = cellfun(f,self.grid_,'un',0);
        end
    end
end

