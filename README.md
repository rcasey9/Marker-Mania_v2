This package only works on Windows. It requires MATLAB (2022b or later) and Vicon Nexus.

Follow these steps for processing your MoCap data:

1) Format your data so that all trials are grouped into folders that contain their respective markerset, Static, and FJC (if you use FJCs) 

*** your static trial must include "Static", "static", or "STATIC" somewhere in its name. No other trials should contain these substrings
     
2) Label your static trial(s), set the region of interest to just 2 total frames, and clear events outside of the region of interest. 
 Run the "static skeleton calibration" and "Set Autolabel Pose" pipelines on your static trial(s).
    
3) If your data quality necessitates Functional Joint Calibrations, label your FJC trial(s) and run the "Functional Joint Calibration" pipeline on them.

4) (OPTIONAL) At this point you may want to check that your setup in Vicon up to this point is sufficient. To test this, run a "Reconstruct And Label" pipeline in Vicon on one of the trials you intend to process. Ensure that it completes without any errors and that Vicon's labeling does an okay job before moving on to MATLAB.

5) Run the Setup.m script. It may prompt you to input the path to the Nexus.exe executable on your pc if it has trouble finding it. This is mandatory before moving on to the next step
 
6) Input the path(s) to your data at the top of the Main.m script and hit run. It can take ~24 hours to process an entire subject depending on the capacity of your computer.
  
7) Input the same path to your data at the top of the Cleanup.m script
 
8) Parse all of the trials in the /Failed/ folder(s) and unlabel any drifiting or out of place markers. Some trials in the /Failed/ folder may not have any issues present.
  
9) Run the Cleanup.m script. All of your processed data should now be in the /Finished/ folder.
