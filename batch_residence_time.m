function batch_residence_time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function "batch_residence_time" helps to build cumulative residence
% time distributions from a single interval imaging. This function is to 
% be used following the peak detection and single-particle tracking (SPT). 
%
% %%%%%%%%%%%%%%%%%%%%%%%
% Prior to this code:
% Peak detection and SPT were performed in Fiji, using Single-Molecule
% Biophysics plugin and the Results Table were saved as 'txt' files.
% The data are organized in the following format (see Example Data):
%   1. Acquisition files, their corresponding bright field images and
%   Results Tables for each interval are stored in a sub-folder.
%   2. A parental folder contains all sub-folders.
%
%  %%%%%%%%%%%%%%%%%%%%%%
% The function of this code:
% 1. Read single particle tracking data that were saved in Results Table
% file (txt). Each Results Table contain all trajectories detected in ONE
% acquisition, at ONE time-laspe interval
%
% 2. For each time-lapse interval, assign each trajectory a unique ID. 
% Append trajectories, and their lifetimes, from all acquisitions and save 
% to LifetimeInFrames.txt file. This is critical for the bootstrapping
% analysis (using bootstrap_global_Fit.m code) where a proportion of the
% trajectories (e.g. 80%) are randomly drawn.
%
% 3. From LifetimeInFrames.txt file, construct the cumulative residence
% time distribution (crtd) in frames from a single experiment. Save the 
% crt to csv file.
%
% 4. Generate crtd and keff*tautl plots
%
% 5. Generate tautl.txt file, containing the specified time-lapse intervals
%
%%%%%%%%%%%%%%%%%%%%%%%
% Please check entries in lines: 48, 52, 59, 61, 64 and 66. Run this code
% first.
%
%%%%%%%%%%%%%%%%%%%%%%%
% % Author: Harshad Ghodke & Han Ho
% %  Date of Last edit: 28 Feb 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
w = warning ('off','all');
%%%%%Interval imaging settings
tau_tl = [0.1; 0.2; 0.5; 1; 2; 4; 8; 10];           %vector containing timelapse times

%%%%%Filehandling - filepath and file names
%%%Starting search path for folder containing 'Analysis' folder
dirpath = 'F:\NER projects\Example data\UvrAYPet dmfd\20180821_091738_415\Analysis';       % path to folder containing Analysis folder

%%%Results table file name descriptors
f_Results = 'Results';
f_PF = '_Peakfit';
f_PT = '_PT';
% Starting frame of the single-molecule phase
f1 = '51';
% Last frame of the single-molecule phase
f2 = '150';
% Step size (in pixel) used in linking trajectories in single particle
% tracking (Fiji)
f_step = '_3';          %minimum distance between two peaks. Setting in Peak Fitter in Fiji/SingleMoleculeBiophysics plugins
f_f1_f2 = strcat('_',f1,'_',f2);
thresh = '_7';         %threshold used for analyzing peaks in peak fitter
ext = '.txt';

%%%Filepath for data. Point to 'Analysis' folder that contains results files.
[filepath] = uigetdir(dirpath);
files = dir(filepath);
subdir = [files.isdir];
subFolders = files(subdir); %list all the sub-folders contained in the parent folder,
                            % including two hidden sub-folders '.' and '..'
%%
% Read Results Table and generate 'LifetimeInFrames.txt' files
for tautl = 3:length(subFolders)
    % obtain directory of the current sub-folder    
    subfolder_path = fullfile(filepath,subFolders(tautl).name);
    cd(subfolder_path);  
    filelist = dir(subfolder_path);       %list of all files in the chosen folder
    
    %Get number of acquisition files, based on the number of bright field
    % images
    n_acq = 0;                  %number of acquisitions
    for i=1:length(filelist)      
        if strfind(filelist(i).name, 'BF')>0
            n_acq = n_acq+1;
        end
    end
    %%
    %%%Read results tables for all acquisitions in the experiment
    for i=1:n_acq
        peak_loc_PT =       [];
        % results_filename: Name of the Results Table generated using Fiji
        % in the previous step of peak localization and tracking.
        results_filename = strcat(f_Results, f_f1_f2,'_', num2str(i), thresh, f_PF, f_PT, f_step, ext);
      
        try
            disp(results_filename);
            if exist(results_filename,'file')==2
                % read Results Table
                [col_PT, peak_loc_PT_temp] = readResultsTable(results_filename);                
                peak_loc_PT(:,1) = peak_loc_PT_temp(:,find(strcmpi('baseline', col_PT)));
                peak_loc_PT(:,2) = peak_loc_PT_temp(:,find(strcmpi('height', col_PT)));
                peak_loc_PT(:,3) = peak_loc_PT_temp(:,find(strcmpi('x', col_PT)));
                peak_loc_PT(:,4) = peak_loc_PT_temp(:,find(strcmpi('y', col_PT)));
                peak_loc_PT(:,5) = peak_loc_PT_temp(:,find(strcmpi('sigma_x', col_PT)));
                peak_loc_PT(:,6) = peak_loc_PT_temp(:,find(strcmpi('sigma_y', col_PT)));
                peak_loc_PT(:,7) = peak_loc_PT_temp(:,find(strcmpi('fwhm_x', col_PT)));
                peak_loc_PT(:,8) = peak_loc_PT_temp(:,find(strcmpi('fwhm_y', col_PT)));
                peak_loc_PT(:,9) = peak_loc_PT_temp(:,find(strcmpi('fwhm', col_PT)));
                peak_loc_PT(:,10) = peak_loc_PT_temp(:,find(strcmpi('error_baseline', col_PT)));
                peak_loc_PT(:,11) = peak_loc_PT_temp(:,find(strcmpi('error_height', col_PT)));
                peak_loc_PT(:,12) = peak_loc_PT_temp(:,find(strcmpi('error_x', col_PT)));
                peak_loc_PT(:,13) = peak_loc_PT_temp(:,find(strcmpi('error_y', col_PT)));
                peak_loc_PT(:,14) = peak_loc_PT_temp(:,find(strcmpi('error_sigma_x', col_PT)));
                peak_loc_PT(:,15) = peak_loc_PT_temp(:,find(strcmpi('error_sigma_y', col_PT)));
                peak_loc_PT(:,16) = peak_loc_PT_temp(:,find(strcmpi('error_fwhm_x', col_PT)));
                peak_loc_PT(:,17) = peak_loc_PT_temp(:,find(strcmpi('error_fwhm_y', col_PT)));
                peak_loc_PT(:,18) = peak_loc_PT_temp(:,find(strcmpi('error_fwhm', col_PT)));
                peak_loc_PT(:,19) = peak_loc_PT_temp(:,find(strcmpi('slice', col_PT)));
                peak_loc_PT(:,20) = peak_loc_PT_temp(:,find(strcmpi('trajectory', col_PT)));
                peak_loc_PT(:,21) = peak_loc_PT_temp(:,find(strcmpi('step_size', col_PT)));
                peak_loc_PT(:,22) = peak_loc_PT_temp(:,find(strcmpi('trajectory_length', col_PT)));
                
                % Filter result table                                
                peak_loc_PT(peak_loc_PT(:,1)<0,:) = [];                                
                peak_loc_PT(peak_loc_PT(:,2)<0,:) = [];
%%
% To facilitate bootstrapping analysis in the next step, each trajectory is
% identified by a unique ID, so that certain proportion of trajectories
% (e.g. 80%) can be drawn in a random manner.
%
% Generate a set of files with trajectory name and length of trajectory
                % Get unique trajectory labels
                C = [];   % empty array to store trajectory label and lifetime in frames                
                [C(:,1),ia,~] = unique(peak_loc_PT(:,20));  % get unique trajectory labels
                C(:,2) = peak_loc_PT(ia(:,1),22);           % identify lifetimes
                % residence time distribution from all acquisitions
                % in the same interval is appended to the LifetimeInFrames.txt
                % file
                if ~isempty(C)
                   dlmwrite (strcat(subfolder_path, filesep,'LifetimeInFrames',thresh, f_f1_f2, f_step, '.txt'),C,'-append');
                end
                disp('LifeTimeinFrames written');
                
            else
                disp('Could not write peak locations file');
            end
        catch ME
            disp(ME)
        end
    end
end

cum_traj_f1 = (1:16)'; % lifetime in frames
cum_traj_f1(:,2:(length(subFolders)-1)) = 0; % columns to store crtd
%%
% export crtd and keff*tautl plot
for tautl = 3:length(subFolders)
    subfolder_path = fullfile(filepath,subFolders(tautl).name);
    cd(subfolder_path);                                              
%  Open file containing lifetime data in frames from all acquistions
    [C] = dlmread(strcat(subfolder_path, filesep,'LifetimeInFrames',thresh, f_f1_f2, f_step, '.txt'));
% residence time distribution from all acquisitions in the same interval
    d1 = C(:,2);  % Obtain histogram in frames
    [n,~] = histcounts(d1);
    n = n';
    k = length(n);
    if length(n) < 16
        n(k:16) = zeros;
    else
        n = n(1:16);
    end
    %  Obtain the cumulative residence time distribution 
    for kk = 1:15
        n(kk) = n(kk) + sum(n(kk+1:16));                     
    end   
    cum_traj_f1(:,tautl-1) = n;   
end
%%
% Plot individual crtd and fit to single exponential model
% using customized function fit_single_exponential
figure(1);
for kk = 1:length(tau_tl)     
     fit_single_exponential(tau_tl(kk,1), cum_traj_f1(:,kk+1),strcat(filepath, filesep,  'SingleExponentialFit1',thresh, f_f1_f2, f_step, '.txt'),'r', 'ro');
end
% Save crtd plot
saveas(gcf, fullfile(filepath,filesep,strcat('Crtd_',thresh,f_f1_f2, f_step,'.fig')));

%   Plot keff*tautl vs tautl 
[ktau1] = dlmread(strcat(filepath, filesep,  'SingleExponentialFit1',thresh, f_f1_f2, f_step, '.txt'));
figure(2);
plot(ktau1(:,1), ktau1(:,10),'r','LineWidth',2);
    
xlabel('\tau_{tl}(s)','Fontsize',20);
ylabel('k_{eff} \tau_{tl}','Fontsize',20);
saveas(gcf, fullfile(filepath,filesep,strcat('kefftautlvstautl plot',thresh,f_f1_f2, f_step,'.fig')));
   
% Save crtd to csv file    
cum_traj_f1 = [[0 tau_tl'];cum_traj_f1];
cum_traj_filename1 = strcat(filepath, filesep, 'Cumulative residence time1', thresh,f_f1_f2, f_step,'.csv');
dlmwrite(cum_traj_filename1, cum_traj_f1);
    
%cd 'C:\Research\van Oijen Lab\Analysis\Stroboscopic fluorescence measurement analysis' %Matlab working directory
end

function fit_single_exponential(tautl, y, filename_str, distcolor, modelcolor)
                    FitParams = [];
                    t = (1:10)'*tautl;
                    hold on;
%                     Plot the dwell time distribution
                    plot(t(2:10,1), y(2:10,1), distcolor,'LineWidth',2);
                    
%                    Fit the dwell time distribution to a single exponential fit
                    options = fitoptions('exp1', 'Lower', [0 -Inf 0 -Inf]);
                    [expFit, gof] = fit(t(2:10,1), y(2:10,1),'exp1',options);
                    FitParams(1,1) = tautl;
                    FitParams(1,2) = expFit.a;                  %amplitude of a*exp(bx)
                    FitParams(1,3) = expFit.b;                  % -keff
                    FitParams(1,4) = gof.rsquare;
                    ci = confint(expFit, 0.95);
                    FitParams(1,5) = ci(1,1);                  %95 percent confidence intervals for a
                    FitParams(1,6) = ci(1,2);                  %95 percent confidence intervals for a
                    FitParams(1,7) = ci(2,1);                  %95 percent confidence intervals for b
                    FitParams(1, 8) = ci(2,2);                 %95 percent confidence intervals for b
                    FitParams(1,9) = -expFit.b;                 %keff
                    FitParams(1,10) = tautl*FitParams(1,9);       %keff*tau_tl
                    
                    ybar1 = (expFit.a).*exp((expFit.b).*t);      %model prediction from fit  
%                     plot model prediction
                    plot(t(2:end,1), ybar1(2:end,1), modelcolor,'LineWidth',2);
                    xlabel('time(s)','Fontsize',20);
                    ylabel('Occurrence','Fontsize',20);
                    
                    dlmwrite(filename_str,FitParams,'-append'); 
end

function [sub,fls] = subfolder(CurrPath,sub,fls)
%------------------------------------------------
tmp = dir(CurrPath);
tmp = tmp(~ismember({tmp.name},{'.' '..'}));
for i = {tmp([tmp.isdir]).name}
   sub{end+1} = [CurrPath '\' i{:}];
   if nargin==2
      sub = subfolder(sub{end},sub);
   else
      tmp = dir(sub{end});
      fls{end+1} = {tmp(~[tmp.isdir]).name};
      [sub, fls] = subfolder(sub{end},sub,fls);
   end% if
end% if
end
function [columns, table] = readResultsTable(filename)
    disp(strcat('Reading file ', filename));
    fid = fopen(filename);

    line = fgetl(fid);
    
    % determine column names
    columns = regexp(line, '(\t|,)', 'split');
    %columns = regexp(line, '\t', 'split');
    
    % put all values into table matrix
    maxRows = 1000000;
    table = zeros(maxRows, length(columns));
    
    line = fgetl(fid);
    row = 0;
    
    while ischar(line) && row < maxRows
        %values = regexp(line, '\t|\n', 'split');
        values = regexp(line, '(,|\t|\n)', 'split');
        row = row + 1;
        
        for i = 1:length(values)
            table(row,i) = str2double(values(i));
        end
        
        line = fgetl(fid);
    end
    
    fclose(fid);
    
    table = table(1:row,:);
    if isempty(table)
        disp('Could not read data from file');
    end
end