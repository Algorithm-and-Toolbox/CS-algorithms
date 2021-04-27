
OMPSBox v1 README
October 18, 2009



OMPSBox installation:
---------------------

1. Unpack the contents of the compressed file to a new directory, named e.g. "ompsbox".
2. If you have not done so before, configure Matlab's MEX compiler by entering
    >> mex -setup
   prior to using MAKE. For optimal performance, it is recommended that you select a compiler
   that performs optimizations. For instance, in Windows, MS Visual Studio is preferred to Lcc.
3. Within Matlab, navigate to the OMPSBox directory, and then to the "private" directory within it,
   and enter MAKE to run the compilation script.
4. Add the OMPSBox package directory to the Matlab path (you can use the ADDPATH command for this).
   Do not add the private directory to the path.


OMPSBox quick start:
--------------------

1. Enter "help sparsedict" at the Matlab prompt for an overview of sparse dictionaries in OMPSBox.
2. Enter "ompspeedcompare" to test the OMPSBox installation and compare the speed of sparse and standard OMP.
3. For a complete list of functions in the package, enter
 >> help ompsbox
This assumes the package was installed to a directory named "ompsbox". If not, replace ompsbox
in the above with the (unqualified) name of the OMPSBox installation directory.

Also see faq.txt for some frequently asked questions about the package.
