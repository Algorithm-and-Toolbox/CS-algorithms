proxSortedL1 <- function(x,lambda)
{
   # Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

   # This file is part of SLOPE Toolbox version 1.0.
   #
   #    The SLOPE Toolbox is free software: you can redistribute it
   #    and/or  modify it under the terms of the GNU General Public License
   #    as published by the Free Software Foundation, either version 3 of
   #    the License, or (at your option) any later version.
   #
   #    The SLOPE Toolbox is distributed in the hope that it will
   #    be useful, but WITHOUT ANY WARRANTY; without even the implied
   #    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   #    See the GNU General Public License for more details.
   #
   #    You should have received a copy of the GNU General Public License
   #    along with the SLOPE Toolbox. If not, see
   #    <http://www.gnu.org/licenses/>.

   # Normalize x
   if (is.complex(x))
   {
      sign = Arg(x)
      x = Mod(x)
   }
   else
   {
     sign = sign(x)
     x = abs(x)
   }
   s = sort(x, decreasing=TRUE, index.return=TRUE)
   s$ix <- s$ix - 1
   n <- length(x)
   
   # Apply prox and inverse ordering
   result <- .C("cproxSortedL1", as.double(s$x), as.double(lambda),
                as.double(vector("double",n)),as.integer(n),
                as.integer(s$ix))[[3]]
   
   # Restore sign or phase
   result <- result * sign
   
   return(result)
}
