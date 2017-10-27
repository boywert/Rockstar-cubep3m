# Rockstar-cubep3m

This is Rockstar code with the extension to run on CubeP3M simulations.

There are 2 file formats (FILE_FORMAT) added to the reading routines.

	- “CUBEP3M” : for old format of CubeP3M (*xv*.dat and *PID*.dat)
	- “ZCUBEP3M” : for new zip format (*zip[0-3]*.dat and *PID*.dat)
	
The required parameters are;

	- CUBEP3M_NDIM  : number of domains per dimension.
	- CUBEP3M_NCDIM : (only for ZCUBEP3M)  
	- CUBEP3M_NP : number of particles per dimension, i.e. the total number of particle in the simulation box is CUBEP3M_NP^3
	- CUBEP3M_PID : = 1 to read particle ID data from *PID* files, = 0 otherwise.

The required format of either CUBEP3M or ZCUBEP3M can be used in the format;

	- FILENAME = “[directory]/<snap>xvPID<block>.dat” (for CUBEP3M), 
  	
  	- FILENAME = “[directory]/node<snap>/<snap>xvPID<block>.dat” (for ZCUBEP3M)

where "<snap>xvPID<block>.dat” will be changed to the correct filename format in the reading routine. The string “xvPID” is the string to be searched and replaced. It can be changed to any strings to suit the application.
