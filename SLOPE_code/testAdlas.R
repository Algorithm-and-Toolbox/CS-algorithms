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

dyn.load('cproxSortedL1.so')
source('proxSortedL1.R')
source('Adlas.R')

A      <- as.matrix(read.table('data/testData_A.txt'     ,header=FALSE,colClasses="numeric"))
b      <- as.matrix(read.table('data/testData_b.txt'     ,header=FALSE,colClasses="numeric"))
lambda <- as.matrix(read.table('data/testData_lambda.txt',header=FALSE,colClasses="numeric"))

colnames(A) <- NULL
colnames(b) <- NULL
colnames(lambda) <- NULL

info <- Adlas(A,b,lambda)

print(info$x)
