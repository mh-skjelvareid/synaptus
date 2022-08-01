% Note: The dataset used here is a full-matrix recording of a set of steel
% pins submerged in water. The pins are examples of point-like scatterers
% and are used to test the full matrix version of the PSM algorithm in a
% single-layer scenario.

%% Close / clear
close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----   SYNAPTUS TEST: array_psm_testFullMatrix_1layer.m    -----')
disp('-----          PSM algorithm adapted for arrays             -----')
disp('-----   Full matrix dataset from linear array, 1 layer      -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core'),fullfile(toolboxPath,'misc'));

%% Load test data
disp('Loading data')
load(fullfile(toolboxPath,'datasets','ArraySteelPins.mat'),'data','MWB')
mwb = MWB*1e-6;         % Measurement window begin

%% Parameters
cc = 1480;              % Water sound speed
fLow = 1e6;             % Lower cutoff frequency
fHigh = 7e6;            % Higher cutoff frequency
fs = 50e6;              % Sampling frequency
thick = 0.05;           % Thickness of water layer (distance from transducer to end of ROI)
aPitch = 0.001;         % Array element pitch
pulseDelay = 0.7e-6;    % Approximate

%% Convert cell to matrix
tmpCell = cell(1,1,length(data));
tmpCell(1,1,:) = data;
ptxw = cell2mat(tmpCell);
[nT,nX,nW] = size(ptxw);

%% Zero-pad according to measurement window begin
ptxw = [zeros(round(mwb*fs),nX,nW); ptxw];

%% Plot example raw data: Received field when middle element is excited
tPlot = ((0:(nT-1))/fs + mwb);
xPlot = (0:(nX-1))*aPitch;
figure
colormap(gray(256))
imagesc(xPlot*1e3,tPlot*1e6,logImage(ptxw(:,:,round(nW/2))))
caxis([-50 0])
colorbar
xlabel('x [mm]')
ylabel('t [us]')
title('Data received by all elements when middle element is excited')

%% Focus
disp('Processing data')
tic 
txDelay = eye(nX) - 1;  % Zero transmit delay, only one element at a time
[im,xIm,zIm] = array_psm(ptxw,fs,txDelay,cc,thick,fLow,fHigh,aPitch,...
    'pulseDelay',pulseDelay);
toc

%%
figure
colormap(gray(256))
imagesc(xIm*1e3,zIm{1}*1e3,logImage(im{1}))
colorbar
caxis([-40 0])
xlabel('x [mm]')
ylabel('z [mm]')
title({'Focused image (dB) based on all transmitted waves'})

%% Add blank line (nicer formatting for test text output).
disp(' ')