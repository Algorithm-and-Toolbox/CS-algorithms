function v = ksvdsver(history)
%KSVDSVER Sparse K-SVD toolbox version information.
%  KSVDSVER displays the current Sparse KSVD toolbox version information.
%
%  KSVDSVER('history') also displays history information about the previous
%  versions of the Sparse KSVD toolbox and their change logs.
%
%  V = KSVDSVER returns the version number of the current KSVD toolbox, and
%  does not display any information.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


ver = 11;

if (nargout>0)
  if (nargin>0)
    error('Invalid number of parameters.');
  end
  v = ver;
  
else
  
  if (nargin==0 || (nargin==1 && strcmpi(history,'history')))
    
    disp(' ');
    disp('-----------------------------------------');
    printf('    Sparse KSVD Toolbox version %d       ',ver);
    disp('-----------------------------------------');
    disp(' ');
    
  else
    error('Unknown parameters.');
  end
  
  if (nargin>0)

    disp(' ');
    disp(' ');
    disp('Sparse K-SVD Toolbox version 11, 18.10.09');
    disp('------------------------------------------');
    disp(' ');
    disp('First public release.');
    disp(' ');
    disp(' ');
    
   end
end

