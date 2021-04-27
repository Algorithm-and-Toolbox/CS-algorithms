function v = ompsver(history)
%OMPSVER OMPSBox version information.
%  OMPSVER displays the current Sparse OMP toolbox version information.
%
%  OMPSVER('history') also displays history information about the previous
%  versions of the Sparse OMP toolbox and their change logs.
%
%  V = OMPSVER returns the version number of the current Sparse OMP
%  toolbox, and does not display any information.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


ver = 1;

if (nargout>0)
  if (nargin>0)
    error('Invalid number of parameters.');
  end
  v = ver;
  
else
  
  if (nargin==0 || (nargin==1 && strcmpi(history,'history')))
    
    disp(' ');
    disp('--------------------------------------');
    printf('     Sparse OMP Toolbox version %d        ',ver);
    disp('--------------------------------------');
    disp(' ');
    
  else
    error('Unknown parameters.');
  end
  
  if (nargin>0)
    
    
    disp(' ');
    disp(' ');
    disp('Sparse OMP Toolbox version 1, 18.10.09');
    disp('---------------------------------------');
    disp(' ');
    disp('First public release.');
    disp(' ');
    disp(' ');
    
    
  end
end
