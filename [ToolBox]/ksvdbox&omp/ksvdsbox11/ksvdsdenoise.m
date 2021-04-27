function [y,A,nz] = ksvdsdenoise(params,msgdelta)
%KSVDSDENOISE Sparse K-SVD denoising.
%  [Y,A] = KSVDSDENOISE(PARAMS) denoises the specified 2-D or 3-D signal
%  using Sparse K-SVD denoising. Y is the denoised signal and A is the
%  sparse representation of the trained dictionary.
%
%  [Y,A] = KSVDSDENOISE(PARAMS,MSGDELTA) specifies the frequency of message
%  printing during the process. MSGDELTA should be a positive number
%  representing the interval in seconds between messages. A zero or
%  negative value cancels all messages. Default is MSGDELTA=5.
%
%  [Y,A,NZ] = KSVDSDENOISE(...) also returns the average number of non-zero
%  coefficients in the representations of the denoised blocks.
%
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'x' - Noisy signal.
%      The 2-D or 3-D signal to denoise. Should be of type double, and (for
%      PSNR computations) with values within [0,1] (to specify a different
%      range, see parameter 'maxval' below).
%
%    'blocksize' - Size of block.
%      Indicates the size of the blocks to operate on. Should be either an
%      array of the form [N1 N2 ... Np], where p=2,3 is the number of
%      dimensions of x, or simply a scalar N, representing the square block
%      [N N ... N]. See parameter 'stepsize' below to specify the amount of
%      overlap between the blocks.
%
%    'dictsize' - Size of dictionary to train.
%      Specifies the number of dictionary atoms to train by Sparse K-SVD.
%
%    'basedict' - Separable base dictionary.
%      Specifies the base dictionary over which the sparse dictionary
%      should be represented. KSVDSDENOISE supports separable base
%      dictionaries of either 2 or 3 dimensions, specified as a 2 or 3
%      element cell array. The base dictionary for the i-th dimension is
%      given in basedict{i}. See 'help sparsedict' for more information.
%
%    'Tdict' - Sparsity of trained atoms.
%      Specifies the number of coefficients used to represent each trained
%      atom over the base dictionary.
%
%    'psnr' / 'sigma' - Noise power.
%      Specifies the noise power in dB (psnr) or the noise standard
%      deviation (sigma), used to determine the target error for
%      sparse-coding each block. If both fields are present, sigma is used
%      unless the field 'noisemode' is specified (below). When specifying
%      the noise power in psnr, make sure to set the 'maxval' parameter
%      as well (below) if the signal values are not within [0,1].
%
%    'trainnum' - Number of training blocks.
%      Specifies the number of training blocks to extract from the noisy
%      signal for Sparse K-SVD training.
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'initA' - Initial sparse dictionary representation matrix.
%      Specifies the initial representation matrix A for the Sparse K-SVD
%      training. Should be either a sparse matrix of size MxL, where M is
%      the number of atoms in the base dictionary, or the string 'data' to
%      initialize by sparse-coding random signal blocks. When a sparse
%      matrix is specified for 'initA', L must be >= dictsize, and in this
%      case the dictionary is initialized using the first dictsize columns
%      from initA. By default, initdict='data'.
%
%    'stepsize' -  Interval between neighboring blocks.
%      Specifies the interval (in pixels/voxels) between neighboring blocks
%      to denoise in the OMP denoising step. By default, all overlapping
%      blocks are denoised and averaged. This can be changed by specifying
%      an alternate stepsize, as an array of the form [S1 S2 ... Sp] (where
%      p=2,3 is the number of dimensions of x). This sets the distance
%      between neighboring blocks to be Si in the i-th dimension. Stepsize
%      can also be a scalar S, corresponding to the step size [S S ... S].
%      Each Si must be >= 1, and, to ensure coverage of the entire noisy
%      signal, size(x,i)-Ni should be a multiple of Si for all i. The
%      default stepsize is 1.
%
%    'iternum' - Number of Sparse K-SVD iterations.
%      Specifies the number of Sparse K-SVD training iterations to perform.
%      If not specified, the default is 10.
%
%    'maxval' - Maximal intensity value.
%      Specifies the range of the signal values. Used to determine the
%      noise standard deviation when the noise power is specified in psnr.
%      By default, the signal values are assumed to be within [0,1]. When
%      'maxval' is specified, this range changes to [0,maxval].
%
%    'memusage' - Memory usage.
%      This parameter controls memory usage of the function. 'memusage'
%      should be one of the strings 'high', 'normal' (default) or 'low'.
%      When 'memusage' is specified, both KSVDS and OMPSDENOISE are invoked
%      using this memusage setting. Note that specifying 'low' will
%      significantly increase runtime.
%
%
%  Optional fields in PARAMS - advanced:
%  -------------------------------------
%
%    'noisemode' - Noise power mode.
%      Specifies whether the 'psnr' or 'sigma' fields should be used to
%      determine the noise power. This is useful when both fields are
%      present in PARAMS. 'noisemode' should be one of the string 'psnr' or
%      'sigma'. If it is not present, and both fields are specified,
%      'sigma' is used.
%
%    'gain' - Noise gain.
%      A positive value (usually near 1) controlling the target error for
%      sparse-coding each block. When gain=1, the target error is precisely
%      the value derived from the psnr/sigma fields. When gain is different
%      from 1, the target error is multiplied by this value. The default
%      value is gain = 1.15.
%
%    'lambda' - Weight of the input signal.
%      Specifies the relative weight attributed to the noisy input signal
%      in determining the output. The default value is 0.1*(maxval/sigma),
%      where sigma is the standard deviation of the noise. See function
%      OMPSDENOISE for more information.
%
%    'maxatoms' - Maximal number of atoms.
%      This parameter can be used to specify a hard limit on the number of
%      atoms used to sparse-code each block. Default value is
%      prod(blocksize)/2, i.e. half the number of samples in a block. See
%      function OMPS2 for more information.
%
%
%   Summary of all fields in PARAMS:
%   --------------------------------
%
%   Required:
%     'x'                      signal to denoise
%     'blocksize'              size of block to process
%     'dictsize'               size of dictionary to train
%     'basedict'               separable base dictionary
%     'Tdict'                  sparsity of trained atoms
%     'psnr' / 'sigma'         noise power in dB / standard deviation
%     'trainnum'               number of training signals
%
%   Optional (default values in parentheses):
%     'initA'                  initial representation matrix A ('data')
%     'stepsize'               distance between neighboring blocks (1)
%     'iternum'                number of training iterations (10)
%     'maxval'                 maximal intensity value (1)
%     'memusage'               'low, 'normal' or 'high' ('normal')
%     'noisemode'              'psnr' or 'sigma' ('sigma')
%     'gain'                   noise gain (1.15)
%     'lambda'                 weight of input signal (0.1*maxval/sigma)
%     'maxatoms'               max # of atoms per block (prod(blocksize)/2)
%
%
%  References:
%  [1] M. Elad and M. Aharon, "Image Denoising via Sparse and Redundant
%      representations over Learned Dictionaries", the IEEE Trans. on Image
%      Processing, Vol. 15, no. 12, pp. 3736-3745, December 2006.
%  [2] R. Rubinstein, M. Zibulevsky, and M. Elad, "Learning Sparse
%      Dictionaries for Sparse Signal Approximation", Technical Report -
%      CS, Technion, June 2009.
%
%  See also KSVDS, OMPSDENOISE, OMPS2, SPARSEDICT.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


%%%%% parse input parameters %%%%%

x = params.x;
blocksize = params.blocksize;
trainnum = params.trainnum;
dictsize = params.dictsize;

p = ndims(x);
if (p<2 || p>3)
  error('KSVDSDENOISE only supports 2-D and 3-D signals.');
end


% blocksize %
if (numel(blocksize)==1)
  blocksize = ones(1,p)*blocksize;
end


% maxval %
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval = 1;
  params.maxval = maxval;
end


% gain %
if (isfield(params,'gain'))
  gain = params.gain;
else
  gain = 1.15;
  params.gain = gain;
end


% msgdelta %
if (nargin<2)
  msgdelta = 5;
end

verbose = 't';
if (msgdelta <= 0)
  verbose='';
  msgdelta = -1;
end


% initial dictionary %

if (~isfield(params,'initA'))
  params.initA = 'data';
end

if (isfield(params,'initA') && ischar(params.initA))
  if (strcmpi(params.initA,'data'))
    params = rmfield(params,'initA');    % causes initialization using random examples
  else
    error('Invalid initial matrix A specified.');
  end
end

if (isfield(params,'initA'))
  params.initA = params.initA(:,1:dictsize);
end


% noise mode %
if (isfield(params,'noisemode'))
  switch lower(params.noisemode)
    case 'psnr'
      sigma = maxval / 10^(params.psnr/20);
    case 'sigma'
      sigma = params.sigma;
    otherwise
      error('Invalid noise mode specified');
  end
elseif (isfield(params,'sigma'))
  sigma = params.sigma;
elseif (isfield(params,'psnr'))
  sigma = maxval / 10^(params.psnr/20);
else
  error('Noise strength not specified');
end

params.Edata = sqrt(prod(blocksize)) * sigma * gain;   % target error for omp
params.codemode = 'error';

params.sigma = sigma;
params.noisemode = 'sigma';


% make sure test data is not present in params
if (isfield(params,'testdata'))
  params = rmfield(params,'testdata');
end


%%%% create training data %%%

ids = cell(p,1);
[ids{:}] = reggrid(size(x)-blocksize+1, trainnum, 'eqdist');
params.data = sampgrid(x,blocksize,ids{:});

% remove dc in blocks to conserve memory %
blocksize = 2000;
for i = 1:blocksize:size(params.data,2)
  blockids = i : min(i+blocksize-1,size(params.data,2));
  params.data(:,blockids) = remove_dc(params.data(:,blockids),'columns');
end


%%%%% KSVDS training %%%%%

if (msgdelta>0)
  disp('Sparse K-SVD training...');
end
A = ksvds(params,verbose,msgdelta);


%%%%%  denoise the signal  %%%%%

if (~isfield(params,'lambda'))
  params.lambda = maxval/(10*sigma);
end

params.A = A;

if (msgdelta>0)
  disp('OMPS denoising...');
end

% call the appropriate ompdenoise function
if (p==2)
  [y,nz] = ompsdenoise2(params,msgdelta);
elseif (p==3)
  [y,nz] = ompsdenoise3(params,msgdelta);
else
  [y,nz] = ompsdenoise(params,msgdelta);
end

end
