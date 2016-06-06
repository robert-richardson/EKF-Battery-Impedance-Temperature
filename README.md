Temperature Estimation using an EKF with Impedance Measurements
===============================================================

<img src="Temperature Estimation using Impedance with Matlab - Github/Figures/GraphicalAbstract.png" width="400" align="right">

This repository contains the Matlab source code to implement an Dual Extended Kalman Filter (DEKF) for battery temperature monitoring using single frequency impedance measurements and a 1-D thermal model.

This code accompanies the following paper:
Richardson, Robert R., and David A. Howey. "Sensorless battery internal temperature estimation using a kalman filter with impedance measurement." Sustainable Energy, IEEE Transactions on 6.4 (2015): 1190-1199. [Publisher copy][1] and [Open-acces pre-print][2]

I would ask that you cite this paper if you want to use this code for your own research.
This code has been developed at the Department of Engineering Science of the University of Oxford. 
For information about our lithium-ion battery research, visit the [Howey Research Group][3] website.
If you are interested in our energy research, check out our research group website [Energy and Power Group][4].

Feel free to email the authors with any questions:  
[Robert Richardson](mailto:robert.richardson@eng.ox.ac.uk)  
[David Howey](mailto:david.howey@eng.ox.ac.uk) 


Requirements
============
You will need MATLAB to run this code. This code has been developed and 
tested in MATLAB R2015b and should work with later versions. 
Although it has not been tested with earlier MATLAB releases, it should 
also work with no or minor modifications.

 
Installation
============

##Option 1 - Downloading a .zip file##
Download a .zip file of the code at:

[https://github.com/robert-richardson/EKF-Battery-Impedance-Temperature/archive/master.zip][5]

Then, unzip the folder in a chosen directory on your computer.

##Option 2 - Cloning the repository with Git##
To clone the repository, you will first need to have [Git][6] installed on 
your computer. Then, navigate to the directory where you want to clone the 
repository in a terminal, and type:
```
git clone https://github.com/robert-richardson/EKF-Battery-Impedance-Temperature.git
```
The folder containing all the files should appear in your chosen directory.


Getting started
===============

Execution of the [MainScript.m](Temperature Estimation using Impedance with Matlab - Github/MainScript.m) file runs the simulation.
The simulation comprises a Dual Extended Kalman Filter for estimation of battery internal temperature distribution with unknown  convection coefficient using single frequency Electrochemical Impedance Spectroscopy (EIS) measurements as measurement input. A 3500 sec drive cycle is simulated, and the estimated core and surface temperatures are compared with 'ground truth' thermocouple measurements. The experimental data was obtained using a 26650 A123 lithium iron phosphate cell.

The MainScript file can also be read using Matlab Publish (rather than simply running the file) to generate a PDF or html file summarizing the file code and outputs. A published pdf resulting from this process is already available in the html subfolder.


License
=======

This open-source MATLAB code is published under the BSD 3-clause License,
please read the `LICENSE.txt` file for more information.


[1]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=7097077&tag=1
[2]: http://arxiv.org/abs/1501.06160
[3]: http://users.ox.ac.uk/~engs1053/
[4]: http://epg.eng.ox.ac.uk/
[5]: https://github.com/robert-richardson/EKF-Battery-Impedance-Temperature/archive/master.zip
[6]: https://git-scm.com/


