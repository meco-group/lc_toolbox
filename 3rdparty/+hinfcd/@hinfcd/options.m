function s = options(obj)
% OPTIONS Returns the standard settings and allows the user to easily change his/her preferred options in the structure
    
    % synthesis
    s.synthesis.gammasolver = 'lmilab';     % 'lmilab', 'cvx', 'yalmip'
    s.synthesis.minrank = false;            % true, false
    s.synthesis.lyapunovsolver = 'lmilab';  % 'lmilab', 'cvx', 'yalmip'
    s.synthesis.relaxation = 1.01;          % any real number >= 1
    s.synthesis.zerotol = 0;                % any real number >= 0
    s.synthesis.formulation = 'statespace'; % 'statespace' or 'descriptor'
    
    % reconstruction
    s.reconstruction.reducedorder = false;  % true, false
    s.reconstruction.solver = 'basiclmi';   % 'lmilab', 'cvx', 'yalmip', 'basiclmi'
    s.reconstruction.method = 'transform';  % 'actual', 'transform'
    s.reconstruction.zerotol = 1e-8;        % any real number >= 0
    s.reconstruction.maxiter = 50;          % any integer number >= 0
    s.reconstruction.singvaltol = 0;        % any real number >= 0
    
    % LMILAB interface
    s.lmilab = [1e-8 500 0 0 0];           % LMILAB options, see the mincx() and feasp() documentation
    
    % YALMIP interface
    try
        s.yalmip = sdpsettings();           % YALMIP options, see the YALMIP documentation
    catch 
        s.yalmip = struct(); 
    end
    
    % CVX interface
    s.cvx.solver = '';                      % solver for CVX, see the CVX documentation 
    s.cvx.solver_settings = {''};           % solver settings for CVX, see the CVX documentation
    s.cvx.precision = 'default';            % precision settings for CVX, see the CVX documentation
    
end