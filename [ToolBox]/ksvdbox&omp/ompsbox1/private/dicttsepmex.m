%This is the Matlab interface to the MEX implementation of the sparse
%dictionary transpose operator. The function is not for independent use,
%only through dicttsep.m.


%DICTTSEPMEX Matlab interface to the sparse dictionary transpose operator
%  Y = DICTTSEPMEX(Bsep,A,X) applies the transpose of the sparse dictionary
%  (Bsep,A) to the specified matrix X, producing Y = (B*A)'*X.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
