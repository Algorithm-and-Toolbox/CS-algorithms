function y = dicttsep(Bsep,A,x)
%DICTTSEP Sparse-separable dictionary transpose operator.
%  Y = DICTTSEP(Bsep,A,X) computes the dictionary operator Y=(B*A)'*X for
%  the sparse dictionary D=B*A and the matrix X. The cell array Bsep and
%  the sparse matrix A represent the sparse dictionary, see 'help
%  sparsedict' for more information.
%
%  Note: the call
%
%    Y = dicttsep(Bsep,A,X)
%
%  is equivalent to the explicit call
%  
%    Y = A'*B'*X
%
%  where B=kron(kron(Bsep{3},Bsep{2}),Bsep{1}). However, the first syntax
%  is significantly faster.
%
%  See also DICTSEP, SPARSEDICT.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


y = dicttsepmex(Bsep,A,x);
