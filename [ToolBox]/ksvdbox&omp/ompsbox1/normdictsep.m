function Aout = normdictsep(Bsep,A)
%NORMDICTSEP Sparse-separable dictionary normalization.
%  Aout = NORMDICTSEP(Bsep,A) normalizes the columns of A such that the
%  sparse-separable dictioanry (Bsep,Aout) has unit-energy atoms.
%
%  See also SPARSEDICT.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2009


atomnorms = zeros(size(A,2),1);
gamma = speye(size(A,2));

for i = 1:size(A,2)
  atomnorms(i) = sqrt(sum(dictsep(Bsep,A,gamma(:,i)).^2));
end

Aout = A*spdiags(1./atomnorms,0,size(A,2),size(A,2));

end
