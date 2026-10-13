
clear all

addpath(fullfile(fileparts(mfilename('fullpath')),'ShorelineS_functions')); %Add shorelineS model functions (folder next to this script)
S=struct;
start = '2000-01-01'; S.reftime = start ;
end_simul = '2024-12-30'; S.endofsimulation = end_simul;
ref_time = datetime(S.reftime);
coastSat_start = datetime('1999-08-01');
n_d = days(datetime(end_simul) - datetime(start));
dy = days(ref_time - coastSat_start); %days since 1 august 1999, this parameter defines the initial coastline of the observed data for comparison
site = 2; %embayment number
beach = 2; %number identifying the beach
interv = 30; %Time interval of storage output (in days)
file_out = '30_sept_int_w'; %folder that will hold the model run output
%plot_waves 
[x, y, x_initial, y_initial, g] = initial_grid([2 3], dy, n_d); %creating the shoreline initial position in coordinates (UTM) and the model grid
g = {g};
%S.xyout=g; %Doing this directly defining the model grid, The model will make a projection of the 'coastline position', 'wave conditions' and 'transports' on this grid. xc,yc are the modeled coordinades and zg is the cross-shore change in relation to the grid
S.xmc = x_initial; %input initial coastline coordinates
S.ymc = y_initial;
S.wvcfile='input_25_contour.nc';
S.interpolationmethod='alongshoremapping';
S.ddeep = 25;
S.dnearshore = 25; 
%S.qscal = 0.5; %default = 1 %calibration factor
S.ds0=75; % initial space step [m]   
S.trform = 'KAMP' ;
S.boundaryconditionstart='periodic';
S.boundaryconditionend='periodic';
S.d=10; % active profile height [m],default = 10 %it is possible to use spacially varying d
S.dt = 1/(365*8); %number of time step. unit = years 
S.tc=0; %fixed time step
S.storageinterval=interv;  % Time interval of storage of output file ('output.mat'; [day])
S.outputdir=file_out;

%%% High angle AND SPIT FORMATION corrections%%
S.suppresshighangle = 0; %If wave direction exceeds a threshold (critical angle), the model forces it back to that threshold instead of letting it go higher.
%S.twopoints=1; %spread the sediment at high-angle transitions also over the second downdrift cell
%S.maxangle=50; %default = 60
%S.relaxationlength = 
%S.spitwidth = 0 ; %default = 50m %Critical spit width that if spit below this value overwash is triggered
%S.spitheadwidth = 20; %default = 200m %Width used for the upwind correction

%video settings
S.plotvisible=1;
S.plotDIR=30;
%S.plotUPW=1;
S.plotHS=30;
S.video=1;
% Plot (and store a video frame) every 240 time steps = 30 days at dt = 3 h.
% The ShorelineS default (plotinterval=1) redraws the figure and stores a
% full-size video frame in memory at every 3-hourly step (~73,000 frames for
% 2000-2024), which dominates the run time and memory use.
% Set to 1 to restore the original behaviour.
S.plotinterval=8*30;
S.fignryear=1;

[S]=ShorelineS(S);
project_output(x, y,file_out,dy,interv, x_initial, y_initial,n_d);       

%%% Add River input as nourishments record%%%
%S.nourish=1;
%S.LDBnourish='C:\Users\mpul348\OneDrive - The University of Auckland\Documents\Doctoral level reseacrh\ShorelineS_Material\shorelines\shorelines-master\01_theoretical_beach\nourishment_data_365_days.nor';
%%%Add headland as groynes##
%S.struct=1;
%S.x_hard=[350,350,400,400];
%S.y_hard=[50,200,200,50];
%S.dirspr=10;
%%%%

%%%Add the approach for bypassing test%%%
%S.revet=1;
%S.x_revet=[300,500,700];
%S.y_revet=[100,150,100];
%S.crit_width=1;
%S.diffraction=1;
%%%%%
%S.Hso=2;                                                                   % wave height [m]
%S.phiw0=345;% with respect to north.                                                               % deep water wave angle [°N]
%S.spread=5;                                                               % time step [year] -> use automatic timestep if S.dt==0 || S.tc==1                                                               % switch for using hard structures
%S.twopoints=1;
%S.seaslope=0.01;
%S.landslope=0.05;
%S.seamin=-3;
%S.landmax=3;
%S.tide_interaction=false;
%S.dx=200;
%S.dy=200;
%S.zdeep=-10;
%S.zshallow=-2;
%S.slope=.002;
%S.smoothfac=0.0;
%S.ld=1000;
%S.growth=0;

                                    

