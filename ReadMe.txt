************** ReadMe ***************
Software requirements
     	1. Fiji, with the smb-plugins added (https://github.com/SingleMolecule/smb-			plugins)
	2. MATLAB 2018a, Mathworks

Analysis steps
	1. Background correction and single-particle tracking. Trajectories from 		individual acquisitions are saved to Results Table (txt files). 
	*** Use Fiji ***
	
	2. For a single experiments, combine all trajectories from multiple acquisitions 	collected at a specified time-lapse intervals.
	*** Use MATLAB: batch_residence_time.m ***
	
	3. From multiple experiments, combine all trajectories and perform bootstrapping 	analysis and global fitting.
	*** MATLAB: bootstrap_global_Fit.m and globalFit.m ***

Functions and outputs of batch_residence_time.m
	
	1. Read single particle tracking data that were saved in Results Table
	file (txt). Each Results Table contain all trajectories detected in ONE
	acquisition, at ONE time-lapse interval
	
	2. For each time-lapse interval, assign each trajectory a unique ID. Append
	trajectories, and their lifetimes, from all acquisitions and save to 
	LifetimeInFrames.txt file. This is critical for the bootstrapping
	analysis (using bootstrap_global_Fit.m code) where a proportion of the
	trajectories (e.g. 80%) are randomly drawn.

	3. From LifetimeInFrames.txt file, construct the cumulative residence
	time distribution (crtd) in frames from a single experiment. Save the 
	crt to csv file.
	
	4. Generate crtd and keff*tautl plots

	5. Generate tautl.txt file, containing the specified time-lapse intervals

Functions and outputs of bootstrap_global_Fit.m and globalFit.m (the two codes are to be placed in the same folders)

	1. Prompt users to provide directories where data from individual experiments are 	stored, and an extra directory to store outputs from these codes.

	2. Read taut.txt file to determine at which time-lapse interval the data were 		acquired.

	3. For data acquired at the same time-lapse interval, combine all trajectories and 	their lengths (in LifetimeInFrames.txt file)
	
	4. Perform n rounds of bootstrap analysis (defined by n_repeat variable). Each 		round, randomize trajectory entries and keep 80% of the trajectories.

	5. From bootstrapped samples, construct cumulative residence time distributions 	(crtds) and global fitting.

	6. Results are output as:
		kb kb_err koff_1 koff_1_err amp_1 amp_1_err koff_2 koff_2_err amp_2
		where 'kb' is the photobleaching rate 'koff_1' and 'koff_2' are off rates
		and 'amp' represents fractional population dissociating according to
		the indicated off rate. 'err' represents standard deviation of the
		bootstrap distribution