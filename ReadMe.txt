The rights of the package are reserved by the compressive sensing lab, Carnegie Mellon University

Author: zhuo hui

C++ code: 
 The package contains the code to capture the Image based on stereo settins and given pointGrey Fly GS3 cameras. 
 In order to adopt to other types of cameras, you can change corresponding lines in sl_depth.cpp. 
 The sl_depth contains more functions which can generate the real time depth map and visulize the results. 
 Moreover, we have also provided the calibration code based on Catech Tool box. 
 Note that the package is based on opencv lib and FLY lib (if you want to utilize the cameras), pls download them first when
 you want to utilize the packpage. 

 
Matlab code: 
 The matlab code attached simulate the same functions as shown in C++ code, and all the codes are based on the generated images. 
 	
 Functions included: 
 1) Coding the structured light images based on grey code.
 2) Decode codemap and calcualte the disparity
 The "runscripts.m" gave the demo....	

 
 If you publish corresponding results based on the code, please cite our submitted work in ECCV, thanks.  
	
 Have fun!
	