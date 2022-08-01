close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----           SYNAPTUS TOOLBOX TEST: psm_test3D           -----')
disp('-----           Phase shift migration, 3D dataset           -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core'),fullfile(toolboxPath,'misc'));

%% Load data
disp('Loading data');
load(fullfile(toolboxPath,'datasets','PlaneScan3D_PlexiAluFBH.mat'),...
    'fs','ptxy','tDelay','fLow','fHigh','thick','xStep','yStep','cc');

%% Focus using PSM
disp('Processing data'); 
tic
[im,xIm,yIm,zIm] = psm(ptxy,fs,tDelay,cc,thick,fLow,fHigh,[xStep yStep]);
toc 

%% Plot raw and focused image
test3D_plot

%% Add blank line (nicer formatting for test text output).
disp(' ')