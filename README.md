 # cytokine-network-inference #
A method for the inference of cytokine interaction networks: IL23 model example code
- - - -
* To run the example code, first download and install the lastest version of Data2Dynamics: https://github.com/Data2Dynamics/d2d, see also https://github.com/Data2Dynamics/d2d/wiki/Installation
* Then download the cytokine-network-inference example code, contained in the folder 'Jansen2021'
* Add the folder 'Jansen2021' to the Data2Dynamics folder 'arFramework3'
* Start MATLAB and manually add the folder 'arFramework3' to the path.
* Run Setup_IL23_model.m (Takes about 15 minutes to run on a Macbook Pro laptop with a 2.3 GHz Intel Core i5 processor.)

* Or, using the command line:
- git clone https://github.com/Data2Dynamics/d2d
- cd Data2Dynamics/d2d/arFramework3
- git clone https://github.com/Joanneke-Jansen/cytokine-network-inference
- /Applications/MATLAB_R2018b.app/bin/matlab -nodisplay (or, for mac: /Applications/MATLAB_R2018b.app/bin/ -nodisplay)
- cd cytokine-network-inference/
- addpath(genpath('../..'))
- Setup_IL23_model
