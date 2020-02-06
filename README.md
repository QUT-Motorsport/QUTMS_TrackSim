# TrackSim
This repo is used to generate a track bounded with cones using GPS data.

This is done using MATLAB.

## To generate map
Place Excel list of latitudes and longitudes in same location as Matlab Code. 

If necessary modify Matlab code to target Excel file. 

Alternatively, modify code to set acceleration_test or skidpad_test to true. 

Change any other relevant parameters in user modifiable variables section, e.g. track width and code. 

Confirm code has produced expected track result with Matlab plots. 

The code will output a .sdf and .config file for the placement of cones for the simulator track data in same location as code.

## To use in Gazebo
The subfolders (e.g. lakeside_test_track) go into `~/.gazebo/models`. 

For example, a properly cloned path would be: `~/.gazebo/models/lakeside_test_track/model.sdf`.

Then you should be able to just insert them in Gazebo or spawn them using .world file. 
