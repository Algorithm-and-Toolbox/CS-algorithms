%This is the Matlab interface to the MEX implementation of the sparse
%dictionary forward operator. The function is not for independent use, only
%through dictsep.m.


%DICTSEPMEX Matlab interface to the sparse dictionary forward operator.
%  X = DICTSEPMEX(Bsep,A,GAMMA) applies the sparse dictionary (Bsep,A) to
%  the specified sparse matrix GAMMA, producing X = B*A*GAMMA.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
