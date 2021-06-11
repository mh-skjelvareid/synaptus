% Note: The dataset used here is a full-matrix recording of a copper block
% with side-drilled holes immersed in water. This script is used to test
% the ability of the algorithm to work on multilayer structures.

%% Close / clear
close all
clearvars
clc

%% Load test data
load('../../datasets/ArrayCuBlockStationary.mat','data','MWB')
mwb = MWB*1e-6;

%% Parameters
xrStep = 0.001;         % Array element pitch (assume all elements are used)
xsStep = 0.001;         % Change to 0.002 to use only every second transducer (for example)
xsrOffset = 0.00;       % Used to simulate offset of sending elements (first sending element not equal to first receiving)
cc = [1480 4660];       % Water/copper sound speed
h = 1;                  % No matched filtering
fLow = 1e6;             % Lower cutoff frequency
fHigh = 7e6;            % Higher cutoff frequency
fs = 50e6;              % Sampling frequency
thick = [0.0235 0.04];           % Thickness of water layer (distance from transducer to end of ROI)

%% Convert cell to matrix
tmpCell = cell(1,1,length(data));
tmpCell(1,1,:) = data;
ptx = cell2mat(tmpCell);
[nT,nXr,nXs] = size(ptx);

%% X and Y axis vectors
tPlot = ((0:(nT-1))/fs + mwb);
xPlot = (0:(nXr-1))*xrStep;

%% Plot example raw data: Received field when middle element is excited
figure
colormap(gray(256))
imagesc(xPlot*1e3,tPlot*1e6,logImage(ptx(:,:,round(nXs/2))))
set(gca,'CLim',[-50 0])
colorbar
xlabel('x [mm]')
ylabel('t [us]')
title('Data received by all elements when middle element is excited')

%% Focus using Multilayer Multistatic Omega-K for 2D
tic
[im,xIm,zIm] = msmlok2d(ptx,fs,mwb,cc,thick,fLow,fHigh,...
    xsStep,xrStep,xsrOffset);
toc

%% Plot focused image
figure
colormap(gray(256))
imagesc(xIm*1000,zIm*1000,logImage(im))
set(gca,'CLim',[-65 -10])
colorbar
xlabel('X pos. [mm]')
ylabel('Z pos. [mm]')
title('Focused image of copper block, showing side-drilled holes')