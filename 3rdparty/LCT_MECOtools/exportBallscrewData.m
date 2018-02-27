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

function exportBallscrewData(TimeSignal)
%EXPORTBALLSCREWDATA Exports TimeSignal into a text file which can be
%parsed by the ball screw interface. 

   [time,val] = signal(TimeSignal);
   data = [time, val];

   [file,path] = uiputfile('inputdata.txt','Export data for ball screw set-up');
   path = fullfile(path, file);

   dlmwrite(path,data,'delimiter','\t');
   
   % !! You should make sure yourself that you assign the right reference
   % values to the right control input for identification !! 
   
end