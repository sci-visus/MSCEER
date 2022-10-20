# MSCEER

### Changes 10/20/2022

***Please NOTE!***

*This is a complete swap of the repository. The initial upload ended up not being the active development path (with improvements and bug fixes, etc.). Given the irreconcilable code divergence between the SVN and Github versions of MSCEER, I wiped the original github repository and replaced with the SVN development version.*

*A lot of stuff is different. This new repository has a ton of projects that either are incomplete, reference other projects that are not available, or were abandoned -- its not sanitized!!*

**However, if you want to compute Morse-Smale complexes, you've come to the right place. Below you will find some descriptions of how to build the more useful parts.**

## Overview

Philosophically, the code is broken into "compute discrete gradient" and then "live analysis", where the first will write a .grad file. Typically this computation is the bulk of the work, but fortunately this only needs to be done once. The live analysis then will read the image, precomputed discrete grad, and construct a MS complex, build a simplification hierarchy -- and then extract the structure you need for analysis at the persistence threshold desired. 

The code is structured into projects that produce small command line applications. For instance, there are apps to build the discrete gradient many different ways:
- steepest_lstar *build discrete gradient using relaxed version of Robin's lower stars algorithm, by sliding a template over the grid -- the fastest of the methods. This ignores full lexicographic order, using just lowest then highest vertex to order cells, attempting to match 0-1 pair direction to reduce 2-manifold "wrinkles".*
- steepest_lstar_tiled *a tiled version of above that may be faster and lower memory (but still in development)*
- steepest *build discrete gradient using relaxed version of Robin's lower stars algorithm.*
- ondemandaccurate *first numerically integrate to get a souce/destination field, then apply constrainted Robin's steepest descent to build gradient that "respects" this. User can choose which asc/dsc manfolds should be computed accurately.*
- ondemandaccurate2d *same as above for 2d examples.*
- convergent *get "accuaracy" through region-growing, using the integrated probability of reaching a critical point to get accurate separatrices.*

After a discrete gradient is computed, to compute the MS complex or any features, use one of the "extractXXX" projects. These really are examples of how to use in your own C++ projects. "extractmsc" probably has the template for what most people want to do - load a discrete gradient, compute the MS complex, build a simplification hierarchy, and iterate over features at simplification thresholds.

## Acknowledgement

**I love feedback**
If you use MSCEER please message me! I will be more than happy to get you up and running!

MSCEER is open source (BSD3). We ask that if you use it in your research to cite:

A. Gyulassy, P. Bremer and V. Pascucci, "Shared-Memory Parallel Computation of Morse-Smale Complexes with Improved Accuracy," in IEEE Transactions on Visualization and Computer Graphics.
doi: 10.1109/TVCG.2018.2864848
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8440824&isnumber=4359476	
	
The IEEE TVCG paper can be found here:
 http://sci.utah.edu/~jediati/Documentation/ParallelAccuratGeom_Gyulassy2018.pdf

And the talk presented at Vis 2018 can be found here:
 http://sci.utah.edu/~jediati/Documentation/IEEEVis2018Talk_Gyulassy.pptx

--------------------------------------
Build
--------------------------------------

The project can be built using CMake 3.11

do the usual cmake thing... 
just activate the projects you want to build.


NOTE: I have tested things on windows and linux is in progress, not mac -- but there are basically no dependencies (unless you enable some of the "extractXXX" projects). Any errors in the compile are typically because Visual Studio lets you get away with leaving off the "this->" and un-instantiated template code is never compiled.  



