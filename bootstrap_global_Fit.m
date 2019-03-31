function bootstrap_global_Fit()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function "bootstrap_global_Fit" helps to pool trajectories from
% multiple experiments and generate bootstrapped samples from the pooled
% population. The bootstrapped samples are globally fitting, and keff*tautl
% plot is generated.
%
% Require globalFit.m to be present in the same folder
% %%%%%%%%%%%%%%%%%%%%%%%
% Prior to this code: Run batch_residence_time.m
% The batch_residence_time.m code generates two types of files that are
% required for this step. 
%
% 1. tautl.txt file containing the time-lapse interval
%
% 2. LifetimeInFrames.txt containing all trajectories (unique ID) and their
% length. When combining, the trajectories can be randomized and only a
% portion of the total combined trajectoris is selected to form the
% bootstrapped sample.
%
%  %%%%%%%%%%%%%%%%%%%%%%
% Readme: When prompted, provide the number of datasets that need to be
% bootstrapped. Provide path to the 'Analysis' folders for each folder.
% For example: 
% '...\Example data\UvrAYPet dmfd\20180821_091738_415\Analysis'
% '...\Example data\UvrAYPet dmfd\20180821_133304_274\Analysis'
% '...\Example data\UvrAYPet dmfd\20180821_174915_379\Analysis'
% Finally, provide the path to the folder where the analysis files should be
% saved. For example: '...\Example data\UvrAYPet dmfd'
%
% Results are output as:
% kb kb_err koff_1 koff_1_err amp_1 amp_1_err koff_2 koff_2_err amp_2
% where 'kb' is the photobleaching rate 'koff_1' and 'koff_2' are off rates
% and 'amp' represents fractional population dissociating according to
% the indicated off rate. 'err' represents standard deviation of the
% bootstrap distribution
%%%%%%%%%%%%%%%%%%%%%%%
% Before running this code, please check entries in lines: 
% 51, 53, 55, 59 and 71 to 134.
%
%%%%%%%%%%%%%%%%%%%%%%%
% % Author: Harshad Ghodke & Han Ho
% %  Date of Last edit: 28 Feb 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
%Select the folders containing the results
dir = 'F:\NER projects\Example data\UvrAYPet dmfd';
answer = inputdlg('Enter the number of datasets');
answer = cell2mat(answer);
thresh = '_7';         %threshold used for peakfitting
% Starting frame of the single-molecule phase
f1 = '51';
% Last frame of the single-molecule phase
f2 = '150';
f_f1_f2 = strcat('_',f1,'_',f2);
% Step size (in pixel) used in linking trajectories in single particle
% tracking (Fiji)
f_step = '_3';
tint = 0.1; %camera integration time (exposure time)
n_repeats = 10; %number of bootstrap samples to be drawn
% Select directories that need to be combined
C = cell(str2num(answer),1);
i=1;
while i<=str2num(answer)
    C{i,1}  = uigetdir(dir);
    i = i + 1;
end
[save_file_path] = uigetdir(dir);

% Specify time-lapse intervals.
d(1).name = 'd1';
d(1).tautl = 0.1;
d(1).data = [];

d(2).name = 'd2';
d(2).tautl = 0.2;
d(2).data = [];

d(3).name = 'd3';
d(3).tautl = 0.3;
d(3).data = [];

d(4).name = 'd4';
d(4).tautl = 0.4;
d(4).data = [];

d(5).name = 'd5';
d(5).tautl = 0.5;
d(5).data = [];

d(6).name = 'd6';
d(6).tautl = 0.6;
d(6).data = [];

d(7).name = 'd7';
d(7).tautl = 0.7;
d(7).data = [];

d(8).name = 'd8';
d(8).tautl = 0.8;
d(8).data = [];

d(9).name = 'd9';
d(9).tautl = 1;
d(9).data = [];

d(10).name = 'd10';
d(10).tautl = 2;
d(10).data = [];

d(11).name = 'd11';
d(11).tautl = 3;
d(11).data = [];

d(12).name = 'd12';
d(12).tautl = 4;
d(12).data = [];

d(13).name = 'd13';
d(13).tautl = 5;
d(13).data = [];

d(14).name = 'd14';
d(14).tautl = 6;
d(14).data = [];

d(15).name = 'd15';
d(15).tautl = 8;
d(15).data = [];

d(16).name = 'd16';
d(16).tautl = 10;
d(16).data = [];

for i = 1:str2num(answer) %Check every directory
    for j = 1:16          %identify the number of time-lapse times and iterate over those
        if exist(fullfile(C{i,1},filesep,num2str(j)), 'dir') == 7
            data1 = [];  %Empty matrix with trajectory label and lifetimes           
            tautl = dlmread(fullfile(C{i,1},filesep,num2str(j),filesep,'tautl.txt'));
            
            if exist(fullfile(C{i,1},filesep,num2str(j),filesep,strcat('LifetimeInFrames',thresh, f_f1_f2, f_step, '.txt')),'file') == 2
                [data1] = dlmread(fullfile(C{i,1},filesep,num2str(j),filesep,strcat('LifetimeInFrames',thresh, f_f1_f2, f_step, '.txt')));
                data1 = data1(:,2);
%                 disp(strcat('Reading',fullfile(C{i,1},filesep,num2str(j),filesep,strcat('LifetimeInFrames',thresh, f_f3_f2, f_step, '.txt'))));
                disp(strcat('Concatenating data from ',fullfile(C{i,1},filesep,num2str(j))));
            else
                disp('Could not find lifetime file');
            end
            
            if d(1).tautl == tautl
                d(1).data = [d(1).data; data1];
            elseif d(2).tautl == tautl
                d(2).data = [d(2).data; data1];
            elseif d(3).tautl == tautl
                d(3).data = [d(3).data; data1];
            elseif d(4).tautl == tautl
                d(4).data = [d(4).data; data1];
            elseif d(5).tautl == tautl
                d(5).data = [d(5).data; data1];
            elseif d(6).tautl == tautl
                d(6).data = [d(6).data; data1];
            elseif d(7).tautl == tautl
                d(7).data = [d(7).data; data1];
            elseif d(8).tautl == tautl
                d(8).data = [d(8).data; data1];
            elseif d(9).tautl == tautl
                d(9).data = [d(9).data; data1];
            elseif d(10).tautl == tautl
                d(10).data = [d(10).data; data1];
            elseif d(11).tautl == tautl
                d(11).data = [d(11).data; data1];
            elseif d(12).tautl == tautl
                d(12).data = [d(12).data; data1];
            elseif d(13).tautl == tautl
                d(13).data = [d(13).data; data1];
            elseif d(14).tautl == tautl
                d(14).data = [d(14).data; data1];
            elseif d(15).tautl == tautl
                d(15).data = [d(15).data; data1];
            elseif d(16).tautl == tautl
                d(16).data = [d(16).data; data1];
            else
                disp('Cannot concatenate data!');
            end
        else
            disp(strcat('no time point found', num2str(j)));
        end
    end
end
%% Plot CRTDs from the pooled data
figure(1);
frequency = [];
for kk = 1:16          
        clear a
        a = d(kk).data;          %Extract vector with lifetimes in frames               
        if numel(d(kk).data)>2         
            tautl_temp = d(kk).tautl;
            %Y_temp is a column vector containing counts as a function of frames;
            [~,Y_temp] = fit_single_exponential2(tautl_temp, a);
            X_temp = tautl_temp*(1:9);
            frequency = [frequency X_temp' Y_temp];
            plot_single_exponential(X_temp', Y_temp,'bo', 'b');            
        end                    
end
saveas(gcf, fullfile(save_file_path,filesep,strcat('CRTD.fig')));
dlmwrite(strcat(save_file_path, filesep, 'CRTD.csv'), frequency);
%% Bootstrap and Global Fitting
fitting_message = [];
ktautl = zeros(n_repeats, 16); %Vector containing keff from fits
% n_repeats is the number of bootstrapped samples to be drawn
for j = 1:n_repeats
    disp(strcat('Bootstrapping sample ',num2str(j)));    
    X = []; Y = [];
    for kk = 1:16          
        clear a
        a = d(kk).data;          %Extract vector with lifetimes in frames
        a = a(randperm(length(a)));  %randomize entries in a
        %     Choose a subset of a to estimate means
        n_bootstrap = ceil(0.8*length(a));               %Select 80 percent of the population
        
        if numel(d(kk).data)>2         
            [ktautl(j,kk),Y_temp] = fit_single_exponential2(d(kk).tautl, datasample(a,n_bootstrap));
            X_temp = d(kk).tautl*(1:9);
            X = [X;X_temp];
            Y = [Y; Y_temp']; 
        end            
    end
    
    %% Fit to mono-exponential model    
    koff1_initial   = 1;
    koff1_lb        = 1e-5;
    n_para          = 2; %two parameters: kb and koff
    [p1(j,:),output1] = globalFit(n_para,X,Y,tint,koff1_initial,koff1_lb);
    %% Fit to bi-exponential model
    n_para = 4; %four parameters: kb, koff1, B1 and koff2 
    [p2(j,:), output2] = globalFit(n_para,X,Y,tint,koff1_initial,koff1_lb);
    koff1_initial = (ktautl(j,end)-ktautl(j,1))/(X(end,1)-X(1,1));
    if p2(j,2) < 1e-3 %|| p2(j,2) <= koff1_initial/2
        disp('guess for koff1 in 2 koff fitting');
        try
            while p2(j,2) < 1e-3              
                %koff1_initial = koff1_initial*3;
                koff1_lb = koff1_lb*10;
                [p2(j,:), output2] = globalFit(n_para,X,Y,tint,koff1_initial*10,koff1_lb);        
            end
        catch
        end
    end
    mes1_temp = strcat('Repeat ', num2str(j), {' '},...
            '1 koff',{' '}, output1);
    mes2_temp = strcat('Repeat ', num2str(j), {' '},...
            '2 koff',{' '}, output2);    
    fitting_message = [fitting_message; mes1_temp; mes2_temp];
end

%% Global Fitting summary
%Fiter koff < 0.001
p1(p1(:,2)<=0.001,:) = [];
p2(p2(:,2)<=0.001,:) = [];
%1 koff model
kb_1 = mean(p1(:,1));
kb_sd_1 = std(p1(:,1));

koff_1 = mean(p1(:,2));
koff_sd_1 = std(p1(:,2));

%2 koff model
kb_2 = mean(p2(:,1));
kb_sd_2 = std(p2(:,1));

koff1_2 = mean(p2(:,2));
koff1_sd_2 = std(p2(:,2));
B1_2 = mean(p2(:,3));
B1_sd_2 = std(p2(:,3));

koff2_2 = mean(p2(:,4));
koff2_sd_2 = std(p2(:,4));
B2_2 = 1 - B1_2;

k_all_1 = [kb_1 kb_sd_1 koff_1 koff_sd_1 1 0 0 0 0];
k_all_2 = [kb_2 kb_sd_2 koff1_2 koff1_sd_2 B1_2 B1_sd_2 koff2_2 koff2_sd_2 B2_2];
k_all   = [k_all_1; k_all_2];

filename1 = strcat(save_file_path, filesep, 'koff_1.csv');
dlmwrite(filename1, p1);
filename2 = strcat(save_file_path, filesep, 'koff_2.csv');
dlmwrite(filename2, p2);
filename3 = strcat(save_file_path, filesep, 'koff_mean.csv');
dlmwrite(filename3, k_all);
%% keff*tautl plot
M = zeros(16,3);
for kk = 1:16
    if ktautl(1,kk)> 0
        M(kk,1) = d(kk).tautl;
        M(kk,2) = mean(ktautl(:,kk));
        M(kk,3) = std(ktautl(:,kk));
    end
    
end

figure(2);
M( ~any(M,2), : ) = [];
dt = cellstr(datetime('now','TimeZone','local','Format','yyyymmddHHMMSS'));
dlmwrite(strcat(save_file_path, filesep, 'Kefftautlvskefftautl_bootstrap_',dt{1,1},'.csv'),M);
shadedErrorBar(M(:,1),M(:,2),M(:,3),'-g',1)
saveas(gcf, fullfile(save_file_path,filesep,strcat('kefftautlvstautl bootstrap plot',dt{1,1},'.fig')));

filename4 = strcat(save_file_path, filesep, 'fitting_mess.csv');
T = cell2table(fitting_message);
writetable(T,filename4);
end

function plot_single_exponential(x, y, distcolor, modelcolor)
                    FitParams = [];
                    hold on;
%                     Plot the dwell time distribution
                    plot(x, y, distcolor,'LineWidth',2);
                    
%                    Fit the dwell time distribution to a single exponential fit
                    options = fitoptions('exp1', 'Lower', [0 -Inf 0 -Inf]);
                    [expFit, gof] = fit(x, y,'exp1',options);
                    FitParams(1,1) = x(1);
                    FitParams(1,2) = expFit.a;                  %amplitude of a*exp(bx)
                    FitParams(1,3) = expFit.b;                  % -keff
                    FitParams(1,4) = gof.rsquare;
                    ci = confint(expFit, 0.95);
                    FitParams(1,5) = ci(1,1);                  %95 percent confidence intervals for a
                    FitParams(1,6) = ci(1,2);                  %95 percent confidence intervals for a
                    FitParams(1,7) = ci(2,1);                  %95 percent confidence intervals for b
                    FitParams(1, 8) = ci(2,2);                 %95 percent confidence intervals for b
                    FitParams(1,9) = -expFit.b;                 %keff
                    FitParams(1,10) = x(1)*FitParams(1,9);       %keff*tau_tl
                    
                    ybar1 = (expFit.a).*exp((expFit.b).*x);      %model prediction from fit  
%                     plot model prediction
                    plot(x, ybar1, modelcolor,'LineWidth',2);
                    xlabel('time(s)','Fontsize',20);
                    ylabel('Occurrence','Fontsize',20);                   
                    

end
function [ktautl,dist] = fit_single_exponential2(tautl, y)

hold on;
%                 Obtain probability distribution from trajectory lifetimes                   
                    x_f = (1:16)';                % lifetime in frames
                    [n,~] = hist(y, x_f);           % y is in frames
                    n = n';
                     
%                 Obtain dwell time distribution in frames
                    for jj = 1:14
                         n(jj, 1) = n(jj,1) + sum(n(jj+1:15,1),1);                     
                    end        
                    
                    FitParams = [];
                    t = (1:10)'*tautl;
                    hold on;
                    dist = n(2:10,1);
                
%                    Fit the dwell time distribution to a single exponential fit
                    options = fitoptions('exp1', 'Lower', [0 -Inf 0 -Inf]);
                    [expFit, gof] = fit(t(2:10,1), n(2:10,1),'exp1',options);
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
                    ktautl = FitParams(1,10);
                    
                    if FitParams(1,4)<0.95                  %Only select good fits
                        ktautl = 0;
                    end
end
function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
% function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
%
% Purpose
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image.
%
% Outputs
% H - a structure of handles to the generated plot objects.
%
%
% Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});
%
% Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1);
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1);
% hold off
%
%
% Rob Campbell - November 2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error checking
error(nargchk(3,5,nargin))


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar)
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:)';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:)';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end

%Set default options
defaultProps={'-k'};
if nargin<4, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end

if nargin<5, transparent=0; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot to get the parameters of the line
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. There we have the option of choosing alpha or a
% de-saturated solid colour for the patch surface .

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=0.15; %How de-saturated or transparent to make patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end


%Calculate the error bars
uE=y+errBar(1,:);
lE=y-errBar(2,:);


%Add the patch error bar
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
    'edgecolor','none',...
    'facealpha',faceAlpha);


%Make pretty edges around the patch.
H.edge(1)=plot(x,lE,'-','color',edgeColor);
H.edge(2)=plot(x,uE,'-','color',edgeColor);

%Now replace the line (this avoids having to bugger about with z coordinates)
delete(H.mainLine)
H.mainLine=plot(x,y,lineProps{:});


if ~holdStatus, hold off, end


if nargout==1
    varargout{1}=H;
end
end
