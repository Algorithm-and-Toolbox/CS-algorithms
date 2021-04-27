%This is the Matlab interface to the Sparse OMP MEX implementation.
%The function is not for independent use, only through omps.m.


%OMPSMEX Matlab interface to the Sparse OMP MEX implementation.
%  GAMMA = OMPSMEX(Bsep,A,X,G,L,SPARSE_G,MSGDELTA,PROFILE) invokes the OMP
%  MEX function according to the specified parameters. The parameter G is
%  not required, and if not specified, should be passed as [].
%
%  L - the target sparsity.
%  SPARSE_G - returns a sparse GAMMA when nonzero, full GAMMA when zero.
%  MSGDELTA - the delay in secs between messages. Zero means no messages.
%  PROFILE - nonzero means that profiling information should be printed.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
