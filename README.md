# CS-algorithms
An Archive of Reconstruction Algorithms for Compressive Sensing and Compressive Imaging

### Archived

- **ISTA&FISAT**: A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems (*SIAM J. Img. Sci. 2, 1*) by Amir Beck and Marc Teboulle, January 2009.[[pdf]](https://dl.acm.org/doi/10.1137/080716542) [[doi]](https://doi.org/10.1137/080716542) [[code_FISTA_matlab]](https://github.com/tiepvupsu/FISTA) [[code_FISTA_python]](https://github.com/JeanKossaifi/FISTA) [[code_denoise_matlab\]](https://github.com/sandeepbanik/Image-denoise-and-TV-solver)[[code_ista_python\]](https://github.com/seunghwanyoo/ista_lasso)
- **TwIST**: Two-Step Iterative Shrinkage/Thresholding Algorithms for Image Restoration) in *IEEE Transactions on Image processing* 2007 by J. Bioucas-Dias and [M. Figueiredo](http://www.lx.it.pt/~mtf/). [[code](http://www.lx.it.pt/~bioucas/code.htm)]
- **GPSR**: Gradient Projection  for Sparse  Reconstruction in *IEEE Journal of Selected Topics in Signal Processing* 2007 by [Mário Figueiredo](http://www.lx.it.pt/~mtf), [Robert D. Nowak](http://www.ece.wisc.edu/~nowak/) and  [Stephen J. Wright](http://www.cs.wisc.edu/~swright/). [[code](http://www.lx.it.pt/~mtf/GPSR/)]
- **OMP_SP_IHT_CoSaMP_GBP_IRLS**: A collection of "CS Recovery Algorithms containing OMP, SP, IHT, CoSaMP, GBP and IRLS" provided by [Chengfu Huo](http://home.ustc.edu.cn/~roy)
- **TVAL3** : TV Minimization by Augmented Lagrangian and Alternating Direction Algorithms, by Chengbo Li and Yin Zhang [[code]](https://www.caam.rice.edu/~optimization/L1/TVAL3/)
- **3DCS**: Three-dimensional (3D) compressive sensing algorithms: towards real-time volumetric imaging, by [Yang Liu](https://github.com/liuyang12) [[code]](https://github.com/liuyang12/3DCS)



### ToolBox

- **ksvdbox&omp**: [Toolboxs]((http://www.cs.technion.ac.il/~ronrubin/software.html#)) for ksvd and omp by [Ron Rubinstein](http://www.cs.technion.ac.il/~ronrubin/software.html#)

  - [OMP-Box v10](./[ToolBox]/ksvdbox&omp/ompbox10/) Implementation of the Batch-OMP and OMP-Cholesky algorithms for quick sparse-coding of large sets of signals [[ref](http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf)].

  - [OMPS-Box v1](./[ToolBox]/ksvdbox&omp/ompsbox1/) Implementation of the Batch-OMP and OMP-Cholesky algorithms for sparse dictionaries [[ref](http://www.cs.technion.ac.il/~ronrubin/Publications/sparsedict.pdf)].

    [KSVD-Box v13](./[ToolBox]/ksvdbox&omp/ksvdbox13/) Implementation of the K-SVD and Approximate K-SVD dictionary training algorithms, and the K-SVD Denoising algorithm [[ref](http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf)]. Requires OMP-Box v10.

    [KSVDS-Box v11](./[ToolBox]/ksvdbox&omp/ksvdsbox11/) Implementation of the sparse K-SVD dictionary training algorithm and the sparse K-SVD Denoising algorithm [[ref](http://www.cs.technion.ac.il/~ronrubin/Publications/sparsedict.pdf)]. Requires OMPS-Box v1. The package is also available without the demo volumes (less recommended) at at [KSVDS-Box v11-min](http://www.cs.technion.ac.il/~ronrubin/Software/ksvdsbox11-min.zip).
  
- **Dictionary Learning Tools for Matlab** : A [toolbox](https://www.ux.uis.no/~karlsk/dle/) for dictionary learning by [Karl Skretting](https://www.ux.uis.no/~karlsk/), University of Stavanger.



### Others

- **IRtools**: **MATLAB package of iterative regularization methods and large-scale test problems.** [[github]](https://github.com/jnagy1/IRtools)
- **SPORCO**: SPORCO is a **Python package for solving optimisation problems with sparsity-inducing regularisation**. These consist primarily of sparse coding and dictionary learning problems, including convolutional sparse coding and dictionary learning, but there is also support for other problems such as Total Variation regularisation and Robust PCA. The optimisation algorithms in the current version are based on the Alternating Direction Method of Multipliers (ADMM) or on the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA).[[github]](https://github.com/bwohlberg/sporco)
- **3DCS**: Three-dimensional (3D) compressive sensing algorithms: towards real-time volumetric imaging, provided by  [Yang Liu](https://github.com/liuyang12).  [[github]](https://github.com/liuyang12/3DCS)
- **Reproducible Deep Compressive Sensing**: Collection of source code for deep learning-based compressive sensing (DCS). Links for source code, pdf, doi are available.Collection of reproducible deep learning for compressive sensing, provided by [Thuong Nguyen](https://github.com/ngcthuong). [[github]](https://github.com/ngcthuong/Reproducible-Deep-Compressive-Sensing)



### Note

1. Please be aware that the codes containing in this repository is downloaded from the internet, and may not be updated to the latest version. For the latest version, please visit the original code website provided above for the updates.
2. Some codes may not contain source information as they are directly collected from 3rd party website.
3. Some algorithms are listed above without providing codes in this repository for their large memory size, but the code link is given above.