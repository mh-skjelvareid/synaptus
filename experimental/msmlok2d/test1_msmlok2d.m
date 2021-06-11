% Note: The dataset used here is a full-matrix recording of steel pins in
% water. This script is used for testing the algorithm in the simplest
% possible case: A single medium, with a small number of scatterers.

%% Close / clear
close all
clearvars
clc

%% Load test data
load('../../datasets/ArraySteelPins.mat','MWB','data')

%% Parameters
xrStep = 0.001;         % Array element pitch (assume all elements are used)
xsStep = 0.001;         % Change to 0.002 to use only every second transducer (for example)
xsrOffset = 0.000;      % Used to simulate offset of sending elements (first sending element not equal to first receiving)
cc = 1480;              % Water sound speed
h = 1;                  % No matched filtering
fLow = 1e6;             % Lower cutoff frequency
fHigh = 7e6;            % Higher cutoff frequency
fs = 50e6;              % Sampling frequency
mwb = MWB*1e-6;         % Measurement window begin
thick = 0.05;           % Thickness of water layer (distance from transducer to end of ROI)

%% Convert cell to matrix
tmpCell = cell(1,1,length(data));
tmpCell(1,1,:) = data;
ptx = cell2mat(tmpCell);
ptx = ptx(:,:,(1 + xsrOffset*1e3):(xsStep*1e3):32);
[nT,nXr,nXs] = size(ptx);

%% X and Y axis vectors
tPlot = ((0:(nT-1))/fs + mwb);
xPlot = (0:(nXr-1))*xrStep;

%% Plot example raw data: Received field when middle element is excited
figure
colormap(gray(256))
imagesc(xPlot*1e3,tPlot*1e6,logImage(ptx(:,:,round(nXs/2))))
set(gca,'CLim',[-40 -0])
colorbar;
xlabel('x [mm]')
ylabel('t [us]')
title('Data received by all elements when middle element is excited')

%% Focus using Multilayer Multistatic Omega-K for 2D
tic
[im,xIm,zIm] = msmlok2d(ptx,fs,mwb,cc,thick,fLow,fHigh,xsStep,xrStep,xsrOffset);
toc

%% Plot focused image
figure
colormap(gray(256))
imagesc(xIm*1000,zIm*1000,logImage(im))
set(gca,'CLim',[-40 0])
colorbar
xlabel('X pos. [mm]')
ylabel('Z pos. [mm]')
