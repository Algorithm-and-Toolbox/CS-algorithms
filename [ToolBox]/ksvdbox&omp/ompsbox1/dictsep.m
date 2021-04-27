function x = dictsep(Bsep,A,Gamma)
%DICTSEP Sparse-separable dictionary forward operator.
%  X = DICTSEP(Bsep,A,GAMMA) computes the dictionary operator X=B*A*GAMMA
%  for the sparse dictionary D=B*A and the sparse matrix GAMMA. The cell
%  array Bsep and the sparse matrix A represent the sparse dictionary, see
%  'help sparsedict' for more information.
%
%  Note: In general, the explicit Matlab call
%  
%    X = kron(kron(Bsep{3},Bsep{2}),Bsep{1})*A*GAMMA
%
%  will be faster than the call
%
%    X = dictsep(Bsep,A,GAMMA)
%
%  due to the sparsity of GAMMA and the highly optimized Matlab
%  implementation (note that this is not the case for DICTTSEP, which is
%  faster than its Matlab code equivalent). Thus, DICTSEP is mostly useful
%  for cases where memory constraints prohibit the computation of the full
%  base dictionary B.
%
%  See also DICTTSEP, SPARSEDICT.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


x = dictsepmex(Bsep,A,Gamma);
