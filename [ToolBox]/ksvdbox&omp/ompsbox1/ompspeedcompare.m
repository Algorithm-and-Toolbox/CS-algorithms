function ompspeedcompare
%OMPSPEEDCOMPARE Speed comparison of standard vs. sparse OMP.
%  OMPSPEEDCOMPARE invokes standard and sparse OMP on a random set of
%  signals, and compares their runtimes. The sparse dictionary assumes a
%  separable base dictionary. The function automatically selects the number
%  of signals for the test based on the speed of the system.
%
%  To run the test, type OMPSPEEDCOMPARE from the Matlab prompt.
%
%  See also SPARSEDICT.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


n = 8; k = 10;  % size of each 1-D dictionary (n x k)
ndims = 3;      % number of dimensions (2 or 3)
                % size of base dictionary is (n^ndims) x (k^ndims)
m = 1000;       % number of atoms in sparse dictionary
Tdict = 16;     % sparsity of dictionary atoms
Tdata = 8;      % sparsity of test signals


% random separable base dictionary %

Bsep = cell(1,ndims);

for i = 1:ndims
  Bsep{i} = rand(n,k);
  Bsep{i} = Bsep{i}*spdiags(1./sqrt(sum(Bsep{i}.^2))',0,size(Bsep{i},2),size(Bsep{i},2));  % normalize
end


% random matrix A %

A = spalloc(k^ndims, m, m*Tdict);
for i = 1:m
  p = randperm(k^ndims);
  A(p(1:Tdict),i) = rand(Tdict,1);
end

A = normdictsep(Bsep, A);


% the full dictionary %

D = dictsep(Bsep, A, speye(size(A,2)));


% check for ompbox installation %

if (exist('omp','file') && exist('ompver','file'))
  ompbox = 1;
else
  disp('Note: OMPBox not detected. Speed comparison to standard OMP will not be performed.');
  ompbox = 0;
end


% select signal number according to computer speed %

if (ompbox)
  x = randn(size(D,1),20);
  tic; omp(D,x,[],Tdata,'messages',-1); t=toc;
  signum = ceil(20/(t/20));     % approximately 20 seconds of OMP-Cholesky
else
  x = randn(size(D,1),300);
  tic; omps(Bsep,A,x,[],Tdata,'messages',-1); t=toc;
  signum = ceil(5/(t/300));     % approximately 5 seconds of sparse OMP-Cholesky
end


% generate random signals %

X = randn(size(D,1),signum);


% run OMP  %

if (ompbox)
  printf('\nRunning OMP-Cholesky...');
  tic; omp(D,X,[],Tdata,'messages',4); t1=toc;

  printf('\nRunning Batch-OMP...');
  tic; omp(D,X,D'*D,Tdata,'messages',1); t2=toc;
end

printf('\nRunning Sparse OMP-Cholesky...');
tic; omps(Bsep,A,X,[],Tdata,'messages',1); t3=toc;

printf('\nRunning Sparse Batch-OMP...');
tic; omps(Bsep,A,X,dicttsep(Bsep,A,dictsep(Bsep,A,speye(size(A,2)))),Tdata,'messages',1); t4=toc;


% display summary  %

printf('\n\nSpeed summary for %d signals, dictionary size %d x %d', signum, size(D,1), size(D,2));
printf('Base dictionary: %d-D separable, atom sparsity %d\n', ndims, Tdict);
printf('Call syntax            Algorithm                Total time');
printf('--------------------------------------------------------------');
if (ompbox)
  printf('OMP(D,X,[],T)          OMP-Cholesky             %5.2f seconds', t1);
  printf('OMP(D,X,G,T)           Batch-OMP                %5.2f seconds', t2);
end
printf('OMPS(Bsep,A,X,[],T)    Sparse OMP-Cholesky      %5.2f seconds', t3);
printf('OMPS(Bsep,A,X,G,T)     Sparse Batch-OMP         %5.2f seconds\n', t4);

