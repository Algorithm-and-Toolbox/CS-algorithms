function make
%MAKE Build the OMPSBox package.
%  MAKE compiles all OMPSBox MEX functions, using Matlab's default MEX
%  compiler. If the MEX compiler has not been set-up before, please run
%
%    mex -setup
%
%  before using this MAKE file.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% detect platform 

compstr = computer;
is64bit = strcmp(compstr(end-1:end),'64');


% compilation parameters

compile_params = cell(0);
if (is64bit)
  compile_params{1} = '-largeArrayDims';
end


% Compile files %

dictsources = {'mexutils.c','myblas.c','sparsedict.c'};
ompsources = {'mexutils.c','ompscore.c','omputils.c','myblas.c','ompprof.c','sparsedict.c'};

disp('Compiling ompsmex...');
mex('ompsmex.c', ompsources{:},compile_params{:});

disp('Compiling omps2mex...');
mex('omps2mex.c', ompsources{:},compile_params{:});

disp('Compiling dictsepmex...');
mex('dictsepmex.c', dictsources{:},compile_params{:});

disp('Compiling dicttsepmex...');
mex('dicttsepmex.c', dictsources{:},compile_params{:});
