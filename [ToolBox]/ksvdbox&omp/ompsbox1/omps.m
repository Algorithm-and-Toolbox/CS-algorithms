function gamma = omps(varargin)
%OMPS Sparsity-constrained Sparse Orthogonal Matching Pursuit.
%  GAMMA = OMPS(Bsep,A,X,G,T) solves the optimization problem
%
%       min  |X - B*A*GAMMA|_2     s.t.  |GAMMA|_0 <= T
%      gamma
%
%  for each of the signals in X, using Batch Orthogonal Matching Pursuit.
%  In the above, Bsep and A are the representation of the sparse-seperable
%  dictionary D = B*A (see 'help sparsedict' for more information), X is a
%  matrix containing column signals, T is the # of non-zeros in each signal
%  representation, and G is the Gramm matrix D'*D = A'*B'*B*A. The output
%  GAMMA is a matrix containing the sparse representations as its columns.
%
%  GAMMA = OMPS(Bsep,A,X,[],T) performs the same operation, but without the
%  matrix G, using OMP-Cholesky. This call produces the same output as
%  Batch-OMP, but is significantly slower. Using this syntax is only
%  recommended when available memory is too small to store G.
%
%  GAMMA = OMPS(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies additional
%  parameters for OMPS. Available parameters are:
%
%    'gammamode' - Specifies the representation mode for GAMMA. Can be
%                  either 'full' or 'sparse', corresponding to a full or
%                  sparse matrix, respectively. By default, GAMMA is
%                  returned as a sparse matrix.
%    'messages'  - Specifies whether progress messages should be displayed.
%                  When positive, this is the number of seconds between
%                  status prints. When negative, indicates that no messages
%                  should be displayed (this is the default).
%    'checkdict' - Specifies whether dictionary normalization should be
%                  verified. When set to 'on' (default) the dictionary
%                  atoms are verified to be of unit L2-norm. Setting this
%                  parameter to 'off' disables verification and accelerates
%                  function performance. Note that an unnormalized
%                  dictionary will produce invalid results.
%    'profile'   - Can be either 'on' or 'off'. When 'on', profiling
%                  information is displayed at the end of the funciton
%                  execution.
%
%
%  Summary of OMPS versions:
%
%     version             |   speed     |   memory
%  -----------------------------------------------------
%   OMPS(Bsep,A,X,G,T)    |  very fast  |  moderate
%   OMPS(Bsep,A,X,[],T)   |  moderate   |  small
%  -----------------------------------------------------
%
%
%  References:
%  [1] R. Rubinstein, M. Zibulevsky, and M. Elad, "Learning Sparse
%      Dictionaries for Sparse Signal Approximation", Technical Report -
%      CS, Technion, June 2009.
%  [2] M. Elad, R. Rubinstein, and M. Zibulevsky, "Efficient Implementation
%      of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit",
%      Technical Report - CS, Technion, April 2008.
%
%  See also OMPS2, SPARSEDICT.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% default options

sparse_gamma = 1;
msgdelta = -1;
checkdict = 1;
profile = 0;


% determine number of parameters

paramnum = 1;
while (paramnum<=nargin && ~ischar(varargin{paramnum}))
  paramnum = paramnum+1;
end
paramnum = paramnum-1;


% parse options

for i = paramnum+1:2:length(varargin)
  paramname = varargin{i};
  paramval = varargin{i+1};

  switch lower(paramname)

    case 'gammamode'
      if (strcmpi(paramval,'sparse'))
        sparse_gamma = 1;
      elseif (strcmpi(paramval,'full'))
        sparse_gamma = 0;
      else
        error('Invalid GAMMA mode');
      end

    case 'messages'
      msgdelta = paramval;
      
    case 'checkdict'
      if (strcmpi(paramval,'on'))
        checkdict = 1;
      elseif (strcmpi(paramval,'off'))
        checkdict = 0;
      else
        error('Invalid checkdict option');
      end

    case 'profile'
      if (strcmpi(paramval,'on'))
        profile = 1;
      elseif (strcmpi(paramval,'off'))
        profile = 0;
      else
        error('Invalid profile mode');
      end

    otherwise
      error(['Unknown option: ' paramname]);
  end
  
end


% get parameters

if (paramnum==5)
  B = varargin{1};
  A = varargin{2};
  X = varargin{3};
  G = varargin{4};
  T = varargin{5};
else
  error('Invalid number of parameters');
end


% verify base dictionary normalization

for i = 1:length(B)
  atomnorms = sum(B{i}.*B{i});
  if (any(abs(atomnorms-1) > 1e-2))
    error('Atoms of B{%d} are not normalized to unit length', i);
  end
end


% verify dictionary normalization

if (checkdict)
  if (isempty(G))
    atomnorms = zeros(1,size(A,2));
    gamma = speye(size(A,2));
    for i = 1:size(A,2)
      atomnorms(i) = sum(dictsep(B,A,gamma(:,i)).^2);
    end
  else
    atomnorms = diag(G);
  end
  if (any(abs(atomnorms-1) > 1e-2))
    error('Dictionary atoms are not normalized to unit length');
  end
end


% omp

gamma = ompsmex(B,A,X,G,T,sparse_gamma,msgdelta,profile);
