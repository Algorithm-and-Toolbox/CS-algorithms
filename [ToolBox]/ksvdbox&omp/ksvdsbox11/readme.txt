
KSVDSBox v11 README
October 18, 2009



KSVDSBox installation:
----------------------

1. Make sure OMPSBox v1 is installed prior to installing this package.
2. Unpack the contents of the compressed file to a new directory, named e.g. "ksvdsbox".
3. If you have not done so before, configure Matlab's MEX compiler by entering
    >> mex -setup
   prior to using MAKE. For optimal performance, it is recommended that you select a compiler
   that performs optimizations. For instance, in Windows, MS Visual Studio is preferred to Lcc.
4. Within Matlab, navigate to the KSVDSBox directory, and then to the "private" directory within it,
   and enter MAKE to run the compilation script.
5. Add the KSVDSBox package directory to the Matlab path (you can use the ADDPATH command for this).
   Do not add the private directory to the path.


KSVDSBox quick start:
---------------------

1. Enter "ksvdsdenoisedemo" at the Matlab command prompt to run a demonstration of Sparse K-SVD
   denoising.
2. For a complete list of functions in the package, enter
 >> help ksvdsbox
This assumes the package was installed to a directory named "ksvdsbox". If not, replace ksvdsbox
in the above with the (unqualified) name of the KSVDSBox installation directory.

