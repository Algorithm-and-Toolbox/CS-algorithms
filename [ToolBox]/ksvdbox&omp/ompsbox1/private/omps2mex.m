%This is the Matlab interface to the Sparse OMP2 MEX implementation.
%The function is not for independent use, only through omps2.m.


%OMPS2MEX Matlab interface to the Sparse OMP2 MEX implementation.
%  GAMMA = OMPS2MEX(Bsep,A,X,G,EPSILON,SPARSE_G,MSGDELTA,MAXATOMS,PROFILE)
%  invokes the OMP2 MEX function according to the specified parameters. The
%  parameter G is not required, and if not specified, [] should be passed.
%
%  EPSILON - the target error.
%  SPARSE_G - returns a sparse GAMMA when nonzero, full GAMMA when zero.
%  MSGDELTA - the delay in secs between messages. Zero means no messages.
%  MAXATOMS - the max number of atoms per signal, negative for no max.
%  PROFILE - nonzero means that profiling information should be printed.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
