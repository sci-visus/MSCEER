MSCEER
=================================

--------------------------------------
Overview
--------------------------------------

MSCEER - This is the initial upload of Morse-Smale complex computation with on-demand accuracy. 

The library is under construction - expect massive changes as I bring more of my code base into MSCEER

#### What exits now:
* MSCEER library - has grids, gradients, morse-smale complexes, etc.
* Discrete gradient compuataion with on-demand accuracy command line
* A test for python access to MSCEER

#### What is coming:
* More functionality for MSCEER
	* Bring in algorithms for alternative discrete gradient compuataion methods
	* A lot more complex ways of querying the MS complex, including geometry extraction
	* Dataflows for complex feature-based pipelines
* Morse stand-alone command line topology apps
	* Additional discrete gradient computation methods
	* command-line fixed-pipeline feature extraction from discrete gradients (e.g. nodes, arcs, etc)
* 2D and 3D visualization tool based on GLUT (probably, although it might be GLFW)
* Documentation: yes, sorry, nothing to see yet, and also, things may be named poorly
* Python wrapping
	
--------------------------------------
Build
--------------------------------------

The project can be build using CMake 3.11

do the usual cmake thing... 

what this builds:
MSCEER.lib
OnDemandAccurate2D.exe
OnDemandAccurate3d.exe
(Lots more to come) 

Dependencies:
* Swig (although i will make this optional)

NOTE: I have tested things on windows and linux, not mac. 

--------------------------------------
Examples included
--------------------------------------

In the ExampleData directory there is a sample file - GaussianBox_50x50x50.raw

Input files to the ondemand accurate discrete gradient are just binary files of size XxYxZxSizeof(float32)

MSCEER indexing has X going fastest, then Y then Z. 

Running OnDemandAccurate2D.exe, or OnDemandAccurate3d.exe with no arguments outputs the usage. 



Th