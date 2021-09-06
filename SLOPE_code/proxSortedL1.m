function x = proxSortedL1(y,lambda)

% Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

% This file is part of SLOPE Toolbox version 1.0.
%
%    The SLOPE Toolbox is free software: you can redistribute it
%    and/or  modify it under the terms of the GNU General Public License
%    as published by the Free Software Foundation, either version 3 of
%    the License, or (at your option) any later version.
%
%    The SLOPE Toolbox is distributed in the hope that it will
%    be useful, but WITHOUT ANY WARRANTY; without even the implied
%    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%    See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the SLOPE Toolbox. If not, see
%    <http://www.gnu.org/licenses/>.

% Normalization
lambda = lambda(:);
y       = y(:);
sgn     = sign(y); % Returns phase for complex numbers
[y,idx] = sort(abs(y),'descend');

% Simplify the problem
k = find(y > lambda,1,'last');

% Compute solution and re-normalize
n = numel(y);
x = zeros(n,1);

if (~isempty(k))
   v1 = y(1:k);
   v2 = lambda(1:k);
   v = proxSortedL1Mex(v1,v2);
   x(idx(1:k)) = v;
end

% Restore signs
x = sgn .* x;
