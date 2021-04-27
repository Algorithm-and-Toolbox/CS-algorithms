function [A,Gamma,err,gerr] = ksvds(params,varargin)
%KSVDS Sparse K-SVD dictionary training.
%  [A,GAMMA] = KSVDS(PARAMS) runs the sparse K-SVD dictionary training
%  algorithm on the specified set of signals, returning the sparse
%  dictionary representation matrix A, and the signal representation
%  matrix GAMMA.
%
%  KSVDS has two modes of operation: sparsity-based and error-based. For
%  sparsity-based minimization, the optimization problem is given by
%
%      min  |X-B*A*GAMMA|_F^2      s.t.  |Gamma_i|_0 <= T
%    A,Gamma
%
%  where B is the base dictionary, X is the set of training signals,
%  Gamma_i is the i-th column of Gamma, and T is the target sparsity. For
%  error-based minimization, the optimization problem is given by
%
%      min  |Gamma|_0      s.t.  |X_i - B*A*Gamma_i|_2 <= EPSILON
%    A,Gamma
%
%  where X_i is the i-th training signal, and EPSILON is the target error.
%
%  [A,GAMMA,ERR] = KSVDS(PARAMS) also returns the target function values
%  after each algorithm iteration. For sparsity-constrained minimization,
%  the returned values are given by
%
%  ERR(A,GAMMA) = RMSE(X,B*A*GAMMA) = sqrt( |X-B*A*GAMMA|_F^2 / numel(X) ).
%
%  For error-constrained minimization, the returned values are given by
%
%    ERR(A,GAMMA) = mean{ |Gamma_i|_0 } = |Gamma|_0 / size(X,2) .
%
%  Error computation slightly increases function runtime.
%
%  [A,GAMMA,ERR,GERR] = KSVDS(PARAMS) computes the target function values
%  on the specified set of test signals as well, usually for the purpose of
%  validation (testing the generalization of the dictionary). This requires
%  that the field 'testdata' be present in PARAMS (see below). The length
%  of ERR and GERR is identical.
%
%  [...] = KSVDS(...,VERBOSE) where VERBOSE is a character string,
%  specifies messages to be printed during the training iterations. VERBOSE
%  should contain one or more of the characters 'i', 'r' and 't', each of
%  which corresponds to a certain piece of information:
%           
%    i - iteration number
%    r - number of replaced atoms
%    t - target function value (and its value on the test data if provided)
%
%  Specifying either 'r', 't' or both, also implies 'i' automatically. For
%  example, KSVDS(PARAMS,'tr') prints the iteration number, number of
%  replaced atoms, and target function value, at the end of each iteration.
%  The default value for VERBOSE is 't'. Specifying VERBOSE='' invokes
%  silent mode, and cancels all messages.
%
%  [...] = KSVDS(...,MSGDELTA) specifies additional messages to be printed
%  within each iteration. MSGDELTA should be a positive number representing
%  the interval in seconds between messages. A zero or negative value
%  indicates no such messages (default). Note that specifying VERBOSE=''
%  causes KSVDS to run in silent mode, ignoring the value of MSGDELTA.
%
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'data' - Training data.
%      A matrix containing the training signals as its columns.
%
%    'basedict' - Separable base dictionary.
%      Specifies the base dictionary over which the sparse dictionary
%      should be represented. KSVDS supports separable base dictionaries of
%      either 2 or 3 dimensions, specified as a 2 or 3 element cell array.
%      The base dictionary for the i-th dimension is given in basedict{i},
%      as a matrix of size Ni X Ki. The size of the full base dictionary is
%      prod{Ni} X prod{Ki}, with prod{Ni} = size(data,1).
%
%    'Tdata' / 'Edata' - Sparse coding target.
%      Specifies the number of coefficients (Tdata) or the target error in
%      L2-norm (Edata) for coding each signal. If only one is present, that
%      value is used. If both are present, Tdata is used, unless the field
%      'codemode' is specified (below).
%
%    'Tdict' - Sparsity of trained atoms.
%      Specifies the number of coefficients used to represent each trained
%      atom over the base dictionary.
%
%    'initA' / 'dictsize' - Initial dictionary / no. of atoms to train.
%      At least one of these two should be present in PARAMS.
%
%      'dictsize' specifies the number of dictionary atoms to train. If it
%      is specified without the parameter 'initA', the dictionary is
%      initialized by sparse-coding dictsize randomly selected training
%      signals over the base dictionary.
%
%      'initA' specifies the initial representation matrix A for the
%      training. It should be either a sparse matrix of size KxM, where K
%      is the number of atoms in the base dictionary, or an index vector of
%      length M, specifying the indices of the examples to sparse-code as
%      initial atoms. If only 'initA' is specified, dictsize is set to M.
%      If 'dictsize' and 'initA' are both present, M must be >= dictsize,
%      and in this case the first dictsize columns from initA are used for
%      initialization.
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'testdata' - Validation data.
%      If present, specifies data on which to compute generalization error.
%      Should be a matrix containing the validation signals as its columns.
%
%    'iternum' - Number of training iterations.
%      Specifies the number of Sparse K-SVD iterations to perform. If not
%      specified, the default is 10.
%
%    'memusage' - Memory usage.
%      This parameter controls memory usage of the function. 'memusage'
%      should be one of the strings 'high', 'normal' (default) or 'low'
%      ('high' and 'normal' are equivalent in this function). When set to
%      'high' or 'normal', the fastest implementation of OMPS is used,
%      which involves precomputing the Gram matrix G=D'*D. This increases
%      speed but also requires more memory. When set to 'low', G is not
%      computed, and the slower version of OMPS is used. The 'low' setting
%      should only be used when the trained dictionary is highly redundant
%      and memory resources are low, as it will dramatically increase
%      runtime. See function OMPS for more information.
%
%    'codemode' - Sparse-coding target mode.
%      Specifies whether the 'Tdata' or 'Edata' fields should be used for
%      the sparse-coding stopping criterion. This is useful when both
%      fields are present in PARAMS. 'codemode' should be one of the
%      strings 'sparsity' or 'error'. If it is not present, and both fields
%      are specified, sparsity-based coding takes place.
%
%
%  Optional fields in PARAMS - advanced:
%  -------------------------------------
%
%    'maxatoms' - Maximal number of atoms in signal representation.
%      When error-based sparse coding is used, this parameter can be used
%      to specify a hard limit on the number of atoms in each signal
%      representation (see parameter 'maxatoms' in OMPS2 for more details).
%
%    'muthresh' - Mutual incoherence threshold.
%      This parameter can be used to control the mutual incoherence of the
%      trained dictionary, and is typically between 0.9 and 1. At the end
%      of each iteration, the trained dictionary is "cleaned" by discarding
%      atoms with correlation > muthresh. The default value for muthresh is
%      0.99. Specifying a value of 1 or higher cancels this type of
%      cleaning completely. Note: the trained dictionary is not guaranteed
%      to have a mutual incoherence less than muthresh. However, a method
%      to track this is using the VERBOSE parameter to print the number of
%      replaced atoms each iteration; when this number drops near zero, it
%      is more likely that the mutual incoherence of the dictionary is
%      below muthresh.
%
%
%   Summary of all fields in PARAMS:
%   --------------------------------
%
%   Required:
%     'data'                   training data
%     'basedict'               separable base dictionary
%     'Tdata' / 'Edata'        sparse-coding target
%     'Tdict'                  sparsity of trained atoms
%     'initA' / 'dictsize'     initial A matrix / dictionary size
%
%   Optional (default values in parentheses):
%     'testdata'               validation data (none)
%     'iternum'                number of training iterations (10)
%     'memusage'               'low, 'normal' or 'high' ('normal')
%     'codemode'               'sparsity' or 'error' ('sparsity')
%     'maxatoms'               max # of atoms in error sparse-coding (none)
%     'muthresh'               mutual incoherence threshold (0.99)
%
%
%  References:
%  [1] R. Rubinstein, M. Zibulevsky, and M. Elad, "Learning Sparse
%      Dictionaries for Sparse Signal Approximation", Technical Report -
%      CS, Technion, June 2009.
%
%  See also KSVDSDENOISE.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


global CODE_SPARSITY CODE_ERROR codemode
global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams

CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%%%% parse input parameters %%%%%


data = params.data;
Tdict = params.Tdict;

basedict = params.basedict;
ndims = numel(basedict);
if (~iscell(basedict))
  error('KSVDS expected a cell array for basedict.');
end
if (ndims<2 || ndims>3)
  error('KSVDS only accepts 2- and 3-D separable base dictionaries.');
end

% base dictionary size
N = prod(cellfun(@(z) size(z,1), basedict));
K = prod(cellfun(@(z) size(z,2), basedict));
if (N ~= size(data,1))
  error('Base dictionary size is incompatible with signal size.');
end

ompparams = {'checkdict','off'};

% coding mode %

if (isfield(params,'codemode'))
  switch lower(params.codemode)
    case 'sparsity'
      codemode = CODE_SPARSITY;
      thresh = params.Tdata;
    case 'error'
      codemode = CODE_ERROR;
      thresh = params.Edata;
    otherwise
      error('Invalid coding mode specified');
  end
elseif (isfield(params,'Tdata'))
  codemode = CODE_SPARSITY;
  thresh = params.Tdata;
elseif (isfield(params,'Edata'))
  codemode = CODE_ERROR;
  thresh = params.Edata;

else
  error('Data sparse-coding target not specified');
end


% max number of atoms %

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
  ompparams{end+1} = 'maxatoms';
  ompparams{end+1} = params.maxatoms;
end


% memory usage %

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% iteration count %

if (isfield(params,'iternum'))
  iternum = params.iternum;
else
  iternum = 10;
end


% omp function %

if (codemode == CODE_SPARSITY)
  ompfunc = @omps;
else
  ompfunc = @omps2;
end


% status messages %

printiter = 0;
printreplaced = 0;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
  if (ischar(varargin{i}))
    verbose = varargin{i};
  elseif (isnumeric(varargin{i}))
    msgdelta = varargin{i};
  else
    error('Invalid call syntax');
  end
end

for i = 1:length(verbose)
  switch lower(verbose(i))
    case 'i'
      printiter = 1;
    case 'r'
      printiter = 1;
      printreplaced = 1;
    case 't'
      printiter = 1;
      printerr = 1;
      if (isfield(params,'testdata'))
        printgerr = 1;
      end
  end
end

if (msgdelta<=0 || isempty(verbose))
  msgdelta = -1; 
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr = (nargout>=3 || printerr);


% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
  testdata = params.testdata;
  if (nargout>=4 || printgerr)
    testgen = 1;
  end
end


% mutual incoherence limit %

if (isfield(params,'muthresh'))
  muthresh = params.muthresh;
else
  muthresh = 0.99;
end
if (muthresh < 0)
  error('invalid muthresh value, must be non-negative');
end


% precomputations %

Ik = speye(K);
if (memusage >= MEM_NORMAL)
  baseG = dicttsep(basedict,Ik,dictsep(basedict,Ik,Ik));
else
  baseG = [];
end


% determine dictionary size %

if (isfield(params,'initA'))
  if (any(size(params.initA)==1) && all(iswhole(params.initA(:))))
    dictsize = length(params.initA);
  else
    dictsize = size(params.initA,2);
  end
end
if (isfield(params,'dictsize'))    % this superceedes the size determined by initA
  dictsize = params.dictsize;
end

if (size(data,2) < dictsize)
  error('Number of training signals is smaller than number of atoms to train');
end


% verify base dictionary normalization %

for i = 1:ndims
  if (any(abs(sum(basedict{i}.*basedict{i})-1) > 1e-2))
    error('Base dictionary atoms are not normalized to unit length');
  end
end


% initialize the dictionary %

if (isfield(params,'initA'))
  if (any(size(params.initA)==1) && all(iswhole(params.initA(:))))
    data_ids = params.initA(1:dictsize);
    A = omps(basedict, Ik, data(:,data_ids), baseG, Tdict, 'checkdict', 'off');
  else
    if (size(params.initA,1)~=K || size(params.initA,2)<dictsize)
      error('Invalid initial matrix A');
    end
    A = sparse(params.initA(:,1:dictsize));
  end
else
  data_ids = find(colnorms_squared(data)>1e-6);   % ensure no zero data elements are chosen
  perm = randperm(length(data_ids));
  A = omps(basedict, Ik, data(:,data_ids(perm(1:dictsize))), baseG, Tdict, 'checkdict', 'off');
end


% normalize the initial dictionary %

A = normdictsep(basedict,A);


err = zeros(1,iternum);
gerr = zeros(1,iternum);

if (codemode == CODE_SPARSITY)
  errstr = 'RMSE';
else
  errstr = 'mean atomnum';
end



%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%


for iter = 1:iternum
  
  G = [];
  if (memusage >= MEM_NORMAL)
    G = A'*baseG*A;
  end
  
  
  %%%%%  sparse coding  %%%%%
  
  Gamma = sparsecode(data,basedict,A,G,thresh);
  
  
  %%%%%  dictionary update  %%%%%
  
  replaced_atoms = zeros(1,dictsize);  % mark each atom replaced by optimize_atom
  
  unused_sigs = 1:size(data,2);  % tracks the signals that were used to replace "dead" atoms.
                                 % makes sure the same signal is not selected twice
  
  p = randperm(dictsize); 
  tid = timerinit('updating atoms', dictsize);
  for j = 1:dictsize
    [A(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(data,basedict,A,p(j),Gamma,baseG,Ik,Tdict,unused_sigs,replaced_atoms);
    Gamma(p(j),data_indices) = gamma_j;
    if (msgdelta>0)
      timereta(tid, j, msgdelta);
    end
  end
  if (msgdelta>0)
    printf('updating atoms: iteration %d / %d', dictsize, dictsize);
  end
  
  
  %%%%%  compute error  %%%%%
  
  if (comperr)
    err(iter) = compute_err(basedict,A,Gamma,data);
  end
  if (testgen)
    if (memusage >= MEM_NORMAL)
      G = A'*baseG*A;
    end
    GammaG = sparsecode(testdata,basedict,A,G,thresh);
    gerr(iter) = compute_err(basedict,A,GammaG,testdata);
  end
  
  
  %%%%%  clear dictionary  %%%%%
  
  [A,cleared_atoms] = cleardict(basedict,A,Gamma,data,baseG,Ik,Tdict,muthresh,unused_sigs,replaced_atoms);
  
  
  %%%%%  print info  %%%%%
  
  info = sprintf('Iteration %d / %d complete', iter, iternum);
  if (printerr)
    info = sprintf('%s, %s = %.4g', info, errstr, err(iter));
  end
  if (printgerr)
    info = sprintf('%s, test %s = %.4g', info, errstr, gerr(iter));
  end
  if (printreplaced)
    info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
  end
  
  if (printiter)
    disp(info);
    if (msgdelta>0), disp(' '); end
  end
  
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            optimize_atom             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [atom,gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(X,basedict,A,j,Gamma,baseG,Ik,Tdict,unused_sigs,replaced_atoms)

% data samples which use the atom, and the corresponding nonzero
% coefficients in Gamma
[gamma_j, data_indices] = sprow(Gamma, j);

if (length(data_indices) < 1)
  maxsignals = 5000;
  perm = randperm(length(unused_sigs));
  perm = perm(1:min(maxsignals,end));
  E = sum((X(:,unused_sigs(perm)) - dictsep(basedict,A,Gamma(:,unused_sigs(perm)))).^2);
  [d,i] = max(E);
  atom = omps(basedict,Ik,X(:,unused_sigs(perm(i))),baseG,Tdict,'checkdict','off');
  d = dictsep(basedict,Ik,atom);
  atom = atom./norm(d);
  gamma_j = zeros(size(gamma_j));
  unused_sigs = unused_sigs([1:perm(i)-1,perm(i)+1:end]);
  replaced_atoms(j) = 1;
  return;
end

smallGamma = Gamma(:,data_indices);
Aj = A(:,j);

E_Gamma_j = collincomb(X,data_indices,gamma_j) - dictsep(basedict,A,sparse(smallGamma*gamma_j')) + dictsep(basedict,Aj,sparse(gamma_j*gamma_j'));
atom = omps(basedict,Ik,E_Gamma_j,baseG,Tdict,'checkdict','off');

d = dictsep(basedict,Ik,atom);
atom = atom./norm(d);
d = d./norm(d);

gamma_j = rowlincomb(d,X,1:size(X,1),data_indices) - dicttsep(basedict,A,d)'*smallGamma + dicttsep(basedict,Aj,d)*gamma_j;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sparsecode               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gamma = sparsecode(data,basedict,A,G,thresh)

global ompfunc ompparams

Gamma = ompfunc(basedict,A,data,G,thresh,ompparams{:});

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             compute_err              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err = compute_err(basedict,A,Gamma,data)
  
global CODE_SPARSITY codemode

if (codemode == CODE_SPARSITY)
  err = sqrt(sum(reperror2(data,basedict,A,Gamma))/numel(data));
else
  err = nnz(Gamma)/size(data,2);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           cleardict                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,cleared_atoms] = cleardict(basedict,A,Gamma,X,baseG,Ik,Tdict,muthresh,unused_sigs,replaced_atoms)

use_thresh = 4;  % at least this number of samples must use the atom to be kept

dictsize = size(A,2);

% compute error in blocks to conserve memory
err = zeros(1,size(X,2));
blocks = [1:3000:size(X,2) size(X,2)+1];
for i = 1:length(blocks)-1
  err(blocks(i):blocks(i+1)-1) = sum((X(:,blocks(i):blocks(i+1)-1) - dictsep(basedict,A,Gamma(:,blocks(i):blocks(i+1)-1))).^2);
end

cleared_atoms = 0;
usecount = sum(abs(Gamma)>1e-7, 2);

for j = 1:dictsize
  
  % compute G(:,j)
  if (numel(baseG)>0)
    Gj = A'*(baseG*A(:,j));
  else
    Gj = dicttsep(basedict,A,dictsep(basedict,Ik,A(:,j)));
  end
  Gj(j) = 0;
  
  % replace atom
  if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
    [y,i] = max(err(unused_sigs));    
    atom = omps(basedict, Ik, X(:,unused_sigs(i)), baseG, Tdict, 'checkdict', 'off');
    d = dictsep(basedict, Ik, atom);
    atom = atom./norm(d);
    A(:,j) = atom;
    unused_sigs = unused_sigs([1:i-1,i+1:end]);
    cleared_atoms = cleared_atoms+1;
  end
  
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,basedict,A,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  err2(blockids) = sum((X(:,blockids) - dictsep(basedict,A,Gamma(:,blockids))).^2);
end

end


function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  Y(blockids) = sum(X(:,blockids).^2);
end

end

