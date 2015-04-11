# RobustFuzzyCMeans

Implementation of Robust Fuzzy C-means algorithm as presented in Dzung Pham "Spatial Models for Fuzzy Clustering", CVIU, 2001. 

```
INPUTS:


data        - NumSamples X NumFeaturesPerSample
NumClusters - number of clusters to make
ImgDim      - 2x1 or 3x1 vector specifying the dimensions of the image
                being segmented, to help find the 4- or 6- or 8-nbrhood for
                each voxel
options     - a structure containing the following optional fields
    options.ExpntQ       - exponent for the matrix U     (Default: 2)
    options.beta    - trade-off between spatial smoothing and FCM objective function
    options.MaxIter - max. number of iterations     (Default: 100)
    options.tol     - tolerance/threshold for min improvement per iteration (Default: 1e-5)
    options.verbose - logical variable, whether to display the info per each iteration (Default: 1)



OUTPUTS:


centers     - centers of the different clusters ( NumClusters x 1)
U           - membership functions ( NumClusters x NumSamples )
ObjFun      - value of objective function ( MaxIter x 1)
```

* Implemented by Pradeep Reddy Raamana. 
 * Comments and bug reports are welcome. 
 * No warranty whatsoever, implied or otherwise.
