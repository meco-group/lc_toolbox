import splines.*
try 
    BSplineBasis([-1 0 1],1);
catch
    error('OptiSpline was not found. Make sure to <a href="matlab:web(''https://github.com/meco-group/optispline/releases/tag/v0.1'',''-browser'')">download the toolbox</a> and to put it on your MATLAB path. Be careful: do NOT add the subfolders of the OptiSpline package to your path.');
end