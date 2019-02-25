
# Resultant Based Incremental Recovery of Camera Pose from Pairwise Matches

Yoni Kasten, Meirav Galun, Ronen Basri (WACV'19)


MATLAB Mex Version 1.0 (2019-02-25)
Copyright (C) Yoni Kasten, Weizmann Institute, 2019.
Licensed for noncommercial research use only.


## Background

This code retrieves a camera pose given 2D correspondences to images with known camera poses.

For more information see:
```
@article{kasten2019resultant,
  title={Resultant Based Incremental Recovery of Camera Pose from Pairwise Matches},
  author={Kasten, Yoni and Galun, Meirav and Basri, Ronen},
  journal={arXiv preprint arXiv:1901.09364},
  year={2019}
}
```

[[arXiv]](https://arxiv.org/abs/1901.09364)
Please cite these paper if you use this code in an academic publication.


## Installation

mex/C++ code:
In order to use the code it is necessary to compile the mex functions.
We supply compiled versions for Windows and Linux.
For compilation, run compileMex.m for compiling ResultantExtractor and extractLinearEquationsT.
Note that the compilation time is long.




## Use

We supply sample configurations for the general case (configurations.mat) and for the 4+2 case - 4 points correspondences come from the same with known camera (configurations4_2.mat). 
We also provide two demo scripts that show how to use each type of configuration
```
exampleGeneral.m
example4_2.m
```

core functions:
 - getRtResultant.m - solves the pose of the unknown camera given 6 correspondences to 6 known cameras where up to 3 correspondences come from the same known camera.
 - getRtResultant_4_2.m - solves the pose of the unknown camera given 6 correspondences to known cameras where 4 correspondences come from the same known camera. 


## License

   This software is provided under the provisions of the Lesser GNU Public License (LGPL). 
   see: http://www.gnu.org/copyleft/lesser.html.

   This software can be used only for research purposes, you should cite
   the aforementioned papers in any resulting publication.

   The Software is provided "as is", without warranty of any kind.




## Version History


* Version 1.0 (2019-02-25)
   Initial Release
