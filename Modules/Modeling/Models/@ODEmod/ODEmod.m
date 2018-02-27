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

classdef (InferiorClasses = {?ss,?tf,?zpk}) ODEmod < Model & AbstractLPVmod
    %ODEMOD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f;      % E*dx/dt = @(x,u) f(x,u) 
        g;      % y = @(x,u) g(x,u)
        E;      % for improper systems
        x0;     % initial state
        Ts = 0; % sampling time
    end
    
    properties (Access=private)
        nu_ = 0;
        ny_ = 0;
        nx_ = 0;
        g_conn = @(x,uf,p) zeros(0,1); % yc = g_conn(x,uf,p)
    end
    
    methods
        
        [y,t,x] = sim(varargin)
        
        varargout = impulse(varargin)
        
        varargout = simfixed(varargin)
        
        varargout = step(varargin)
            
        function self = ODEmod(varargin)
            % ODEMOD (nonlinear) general ODE model
            % Possible input combinations
            % 1) ODE: (p can be omitted in case of systems which do not depend on a parameter 
            %   - f(x,u,p): state equation so that E*dot(x) = f(x,u,p)
            %   - g(x,u,p): output equation so that y = g(x,u,p)
            %   - E(x,u,p): mass matrix of the ODE (when omitted, E=I)
            %   - p: cell of SchedulingParameters
            %   - size: [nin, nout, nx]
            % 2) Linear model (LPVLFT,LPVDSS,LTIDSS,LTILFT)
            if isa(varargin{1},'LTIDSSmod')
                f_ = @(x,u,p) varargin{1}.a*x + varargin{1}.b*u;
                g_ = @(x,u,p) varargin{1}.c*x + varargin{1}.d*u;
                E_ = @(x,u,p) varargin{1}.e;
                p_ = varargin{1}.parameters();
                name_ = varargin{1}.name;
                siz = [varargin{1}.nout(),varargin{1}.nin(),varargin{1}.nx];
                x0_ = varargin{1}.x0;
                Ts_ = varargin{1}.Ts;
            elseif isa(varargin{1},'LTILFTmod')
                varargin{1} = simplify(varargin{1});
                f_ = @(x,u,p) varargin{1}.a*x + varargin{1}.b*u;
                g_ = @(x,u,p) varargin{1}.c*x + varargin{1}.d*u;
                E_ = @(x,u,p) varargin{1}.e;
                p_ = varargin{1}.parameters();
                name_ = varargin{1}.name;
                siz = [varargin{1}.nout(),varargin{1}.nin(),varargin{1}.nx];
                x0_ = varargin{1}.x0;
                Ts_ = varargin{1}.Ts;
            elseif isa(varargin{1},'LPVDSSmod')
                p_ = varargin{1}.arguments();
                f_ = @(x,u,p) safeeval(varargin{1}.a,p,p_)*x + safeeval(varargin{1}.b,p,p_)*u;
                g_ = @(x,u,p) safeeval(varargin{1}.c,p,p_)*x + safeeval(varargin{1}.d,p,p_)*u;
                E_ = @(x,u,p) safeeval(varargin{1}.e,p,p_);
                p_ = varargin{1}.parameters();
                name_ = varargin{1}.name;
                siz = [varargin{1}.nout(),varargin{1}.nin(),varargin{1}.nx];
                x0_ = varargin{1}.x0;
                Ts_ = varargin{1}.Ts;
            elseif isa(varargin{1},'LPVLFTmod')
                p_ = varargin{1}.arguments();
                S_ = @(p) simplify(evalme(varargin{1},p,p_));
                f_ = @(x,u,p) A(S_(p))*x + B(S_(p))*u;
                g_ = @(x,u,p) C(S_(p))*x + D(S_(p))*u;
                E_ = @(x,u,p) e(S_(p));
                p_ = varargin{1}.parameters();
                name_ = varargin{1}.name;
                siz = [varargin{1}.nout(),varargin{1}.nin(),varargin{1}.nx];
                x0_ = varargin{1}.x0;
                Ts_ = varargin{1}.Ts;
            elseif isa(varargin{1},'ODEmod')
                f_ = varargin{1}.f;
                g_ = varargin{1}.g;
                E_ = varargin{1}.E;
                p_ = varargin{1}.parameters();
                name_ = varargin{1}.name;
                siz = [varargin{1}.nout(),varargin{1}.nin(),varargin{1}.nx];
                x0_ = varargin{1}.x0;
                Ts_ = varargin{1}.Ts;
            else
                f_ = varargin{1};
                g_ = varargin{2};
                varargin = varargin(3:end);
                
                isfun = cellfun(@(x)isa(x,'function_handle'),varargin);
                iscel = cellfun(@iscell,varargin);
                isparam = cellfun(@(x)isa(x,'SchedulingParameter'),varargin);
                isnum = cellfun(@isnumeric,varargin);

                if any(iscel)
                    assert(sum(iscel)==1,'ODEmod expects one cell argument only.');
                    p_ = varargin{iscel};
                elseif any(isparam)
                    assert(sum(isparam)==1,'ODEmod expects one SchedulingParameter argument only. Otherwise, use cell');
                    p_ = varargin(isparam);
                else
                    p_ = {};
                end
                
                if nargin(f_) == 2, f_ = @(x,u,p) f_(x,u); end
                if nargin(g_) == 2, g_ = @(x,u,p) g_(x,u); end

                if any(isnum)
                    assert(sum(isnum)<3,'ODEmod expects at most 2 numerical inputs: the size [ny,nu,nx] and, if discrete, the sampling time Ts.');
                    issiz3 = cellfun(@(x) length(x)==3, varargin(isnum));
                    issiz1 = cellfun(@(x) length(x)==1, varargin(isnum));
                    comb = sum([issiz3 ; issiz1],2);
                    assert(max(comb)==1,'The size should have three entries [ny,nu,nx] and the sampling time Ts can only be scalar.');
                    numargs = varargin(isnum); 
                    if any(issiz3); siz = numargs{issiz3}; else; [siz(1),siz(2),siz(3)] = ODEmod.sizeguess(f_,g_,p_); end
                    if any(issiz1); Ts_ = numargs{issiz1}; else; Ts_ = 0; end
                else
                    Ts_ = 0;
                    [siz(1),siz(2),siz(3)] = ODEmod.sizeguess(f_,g_,p_);
                end

                if any(isfun)
                    assert(sum(isfun)==1,'ODEmod expects one algebraic function argument only.');
                    E_ = varargin{isfun};
                else
                    E_ = [];
                end
                
                name_ = [];
                x0_ = zeros(siz(3), 1);
            end
            
            self@AbstractLPVmod(p_);
            self.f = f_;
            self.g = g_;
            self.Ts = Ts_;
            
            self.ny_ = siz(1);
            self.nu_ = siz(2);
            self.nx_ = siz(3);
            
            if isempty(E_), E_ = @(x,u,p) eye(self.nx_,self.nx_);
            elseif nargin(E_) == 2, E_ = @(x,u,p) E_(x,u); end
            self.E = E_;
            self.x0 = x0_;
            self.name = name_;
        end
        
        function blk = blkdiag(varargin)
            function e = blkdiagc(c)
                e = c{1};
                for j = 2:length(c)
                    ck = c{j};
                    e = @(x,u,p) blkdiag(e(x,u,p),ck(x,u,p));
                end
            end
                
            if all(cellfun(@(x)isa(x,'ODEmod'),varargin))
                p = mergeparameters(varargin{:});
                sp = cellfun(@(x)SchedulingParameter.ismember(parameters(x),p),varargin,'un',0);
                
                ix = 0; iu = 0; iy = 0;
                f_ = cell(size(varargin));
                g_ = cell(size(varargin));
                g_conn_ = cell(size(varargin));
                E_ = cell(size(varargin));
                for k = 1:nargin
                    s = varargin{k};
                    f_{k} = @(x,u,p) s.f(x(ix+(1:s.nx)),u(iu+(1:s.nin())),p(sp{k}));
                    g_{k} = @(x,u,p) s.g(x(ix+(1:s.nx)),u(iu+(1:s.nin())),p(sp{k}));
                    g_conn_{k} = @(x,u,p) s.g_conn(x(ix+(1:s.nx)),u(iu+(1:s.nin())),p(sp{k}));
                    E_{k} = @(x,u,p) s.E(x(ix+(1:s.nx)),u(iu+(1:s.nin())),p(sp{k}));
                    ix = ix + s.nx;
                    iu = iu + s.nin();
                    iy = iy + s.nout();
                end
                
                f_ = @(x,u,p) cell2mat(cellfun(@(c)c(x,u,p),f_','un',0));
                g_ = @(x,u,p) cell2mat(cellfun(@(c)c(x,u,p),g_','un',0));
                g_conn_ = @(x,u,p) cell2mat(cellfun(@(c)c(x,u,p),g_conn_','un',0));
                E_ = blkdiagc(E_);
                blk = ODEmod(f_,g_,E_,p,[iy,iu,ix]);
                blk.g_conn = g_conn_;
                
                % Combine initial states
                x0_ = reshape(cell2mat(cellfun(@(x) x.x0, varargin, 'un', 0)'),[sum(cell2mat(cellfun(@(x) x.nx, varargin, 'un', 0))),1]); 
                blk = blk.setx0(x0_);
                
            else
                blk = varargin{1};
                for i = 2:nargin
                    blk = blkdiag_model(blk,varargin{i});
                end
            end
        end
        
        function self = lft(self,other,nu,ny)  
            function ret = dispret(label,val)
                % for debugging
                fprintf([label ': ']); disp(num2str(val));
                ret = val;
            end
            
            if ~isempty(other) % return other when empty
                if nargin == 2
                    if all(size(other) <= size(self))
                        [nu,ny] = size(other);
                    else
                        [nu,ny] = size(self);
                    end
                end
                
                blk = blkdiag(self,other);

                % Separate connected and unconnected stuff
                s1 = size(self);
                s2 = size(other);
                freein = [1:(s1(2)-nu),s1(2)+((ny+1):s2(2))];
                freeout = [1:(s1(1)-ny),s1(1)+((nu+1):s2(1))];
                connin = [s1(2)-nu+(1:nu),s1(2)+(1:ny)]; % connin(1),connin(2)
                connout = [s1(1)+(1:nu),s1(1)-ny+(1:ny)]; % connout(2),connout(1)

                Duc = zeros(nin(blk),length(connin));
                Duc(sub2ind(size(Duc),connin,1:length(connin))) = 1;
                Duf = zeros(nin(blk),length(freein));
                Duf(sub2ind(size(Duf),freein,1:length(freein))) = 1;
                Syc = zeros(length(connout),nout(blk));
                Syc(sub2ind(size(Syc),1:length(connout),connout)) = 1;
                Syf = zeros(length(freeout),nout(blk));
                Syf(sub2ind(size(Syf),1:length(freeout),freeout)) = 1;
                Sx = [eye(nx(blk)),zeros(nx(blk),length(connin))];
                Suc = [zeros(length(connin),nx(blk)),eye(length(connin))];
                
                % Compute new f,g,E
                f_ = @(x,u,p) vertcat(blk.f(x(1:blk.nx),Duc*x(blk.nx+(1:(nu+ny)))+Duf*u,p),... 
                                  Syc*blk.g(x(1:blk.nx),Duc*x(blk.nx+(1:(nu+ny)))+Duf*u,p)-x(blk.nx+(1:(nu+ny))));
                g_ = @(x,u,p) Syf*blk.g(x(1:blk.nx),Duc*x(blk.nx+(1:(nu+ny)))+Duf*u,p);
                g_conn_ = @(x,uf,p) vertcat(blk.g_conn(x(1:blk.nx), Duf*uf+Duc*Suc*x, p), ...
                                            Syc*blk.g(x(1:blk.nx), Duf*uf+Duc*Suc*x, p));
                E_ = @(x,u,p) blkdiag(blk.E(x(1:blk.nx),Duc*x(blk.nx+(1:(nu+ny)))+Duf*u,p),zeros(nu+ny,nu+ny));
                siz = [blk.nout-ny-nu,blk.nin-ny-nu,blk.nx+ny+nu];
                
                % Add 'internal' states arising from the connections 
                x0_ = [blk.x0 ; NaN*ones(ny+nu,1)]; % NaN will be replaced by the appropriate numerical value at simulation time
                
                % Make model
                self = ODEmod(f_,g_,E_,parameters(blk),siz);
                self = self.setx0(x0_);
                self.g_conn = g_conn_;
            end
        end
        
        function n = nin(self)
            n = self.nu_;
        end
        
        function n = nout(self)
            n = self.ny_;
        end
        
        function n = nx(self)
            n = self.nx_;
        end
        
        function self = std(self)
            warning('No matlab standard for ODE model');
        end
        
        function self = evalme(self,p)
            assert(length(p) == 1,'only one parameter for now');
            f_ = @(x,u) self.f(x,u,p);
            g_ = @(x,u) self.g(x,u,p);
            siz = [self.nout,self.nin,self.nx_];
            if ~isempty(self.E)
                E_ = @(x,u) self.E(x,u,p);
                self = ODEmod(f_,g_,E_,siz);
            else
                self = ODEmod(f_,g_,siz);
            end
        end
        
        function self = setx0(self,x0_)
            assert(length(x0_) == self.nx() && any(size(x0_)==1), ['The number of initial states (' num2str(length(x0_)) ') does not match the order of your system (' num2str(self.nx) ').']);
            self.x0 = x0_(:);
        end
        
        function self = c2d(self,Ts,varargin)
            assert(nargin>1, 'c2d requires at least two input arguments.');
            assert(self.Ts == 0, 'Your system is already discrete!');
            if nargin==2
                assert(isnumeric(Ts) && isscalar(Ts) && Ts>0, 'The sampling time should be a positive number.');
                method = 'fweuler';
            else
                method = varargin{1};
                assert(any(strcmp(varargin{1},{'fweuler'})),'Only ''fweuler'' is supported for the time being.');
            end
            
            switch method
                case 'fweuler'
                    % discretize
                    self.Ts = Ts;
                    self.f = @(x,u,p) self.f(x,u,p)*Ts + x;
                    
                    % initial states remain the same
            end
        end
    end

    methods(Access=protected)
        function sub = submodel(self,idxout,idxin)
            Sy = zeros(length(idxout),nout(self));
            Sy(sub2ind(size(Sy),1:length(idxout),idxout)) = 1;
            Su = zeros(nin(self),length(idxin));
            Su(sub2ind(size(Su),idxin,1:length(idxin))) = 1;
            
            f_ = @(x,u,p) self.f(x,Su*u,p);
            g_ = @(x,u,p) Sy*(self.g(x,Su*u,p));
            E_ = @(x,u,p) self.E(x,Su*u,p);
            sub = ODEmod(f_,g_,E_,self.parameters(),[length(idxout),length(idxin),self.nx_]);
            sub = sub.setx0(self.x0);
        end
    end
    
    methods(Access=private)
        
        function s = siz(self)
            s = [self.ny_,self.nu_,self.nx_];
        end
        
    end
    
    methods(Static,Access=private)
        
        function tfinal = findss(F,x0,t0)
            h = 1e-3;
            Fp = @(xh) xh-x0-h*F(t0,xh);
            xh = fsolve(Fp,x0);
            tau = max(abs(xh./(xh-x0)));
            tfinal = t0+3*tau;
            
%             dotx = F(t0,x0)./x0;
%             tau = max(abs(1./dotx));
%             tfinal = t0+3*tau;
        end
        
        function [ny_,nu_,nx_] = sizeguess(f,g,p)
            nx_ = -1;
            for total = 0:10
                if nx_~= -1, break; end
                for NX = 1:(total-1)
                    if nx_~=-1, break; end
                    x0 = zeros(NX,1);
                    u0 = zeros(total-NX,1);
                    p0 = zeros(size(p,1),1);
                    try
                        dx0 = f(x0,u0,p0);
                        assert(length(dx0) == length(x0));
                        y0 = g(x0,u0,p0);
                        ny_ = length(y0);
                        nu_ = length(u0);
                        nx_ = length(x0);
                    catch
                        %nothing
                    end
                end 
            end
            assert(nx_~=-1,['Could not guess the model size. ' ...
                'This is normal if the number of states + number of inputs ' ...
                'exceeds the maximum value (10). In that case, set the number ' ...
                'of inputs and outputs manually. This is also possible if the number '...
                'of states is 0. Try pass the size along with f,g,(p).']);
        end
    end

end