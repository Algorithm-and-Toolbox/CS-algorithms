%
%Overview of sparse dictionaries and sparse OMP.
%
%  Introduction to sparse dictionaries:
%
%  A sparse dictionary is a dictionary which has a sparse structure of the
%  form D=B*A, where B is a fixed base dictionary and A is a sparse matrix.
%  The base dictionary B typically has a fast algorithmic implementation -
%  specifically, faster than explicit matrix multiplication - which makes
%  the sparse dictionary very efficient to apply. Sparse dictionaries
%  combine the efficiency of fast transforms with the adaptability of an
%  explicit matrix representation, and provide a simple and powerful
%  dictionary structure. For more information, see the reference below.
%
%  Sparse dictionaries in OMPSBox:
%
%  OMPSBox supports sparse dictionaries with a separable base dictionary of
%  either 2 or 3 dimensions. This means that if X=B*A*GAMMA, where X is a
%  matrix of column signals, B*A is the sparse dictionary, and GAMMA is a
%  matrix of sparse representations, then each column in X represents a
%  vectorized 2-D or 3-D signal of size N1 X N2 or N1 X N2 X N3,
%  respectively, and similarly, each column in B represents a vectorized
%  2-D or 3-D atom of the same size. 
%
%  A sparse dictionary in OMPSBox is represented as a pair (Bsep,A). The
%  base dictionary B is specified by the 2- or 3-element cell array Bsep,
%  where each Bsep{i} is the base dictionary for the i-th dimension, and is
%  a matrix of size Ni X Ki. The full base dictionary is explicitly given
%  by B = kron(Bsep{2},Bsep{1}) or B = kron(kron(Bsep{3},Bsep{2}),Bsep{1}),
%  respectively, and is of size prod{Ni} X prod{Ki}. The matrix A is a
%  sparse matrix of size prod{Ki} X M , containing the sparse
%  representations of the dictionary atoms over B.
%
%  References:
%  [1] R. Rubinstein, M. Zibulevsky, and M. Elad, "Learning Sparse
%      Dictionaries for Sparse Signal Approximation", Technical Report -
%      CS, Technion, June 2009.
%
%
%  See also DICTSEP, DICTTSEP, NORMDICTSEP, OMPS, OMPS2.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


disp('Enter ''help sparsedict'' for an overview of sparse dictionaries and OMPSBox.');
