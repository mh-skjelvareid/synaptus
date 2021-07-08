close all
clearvars
clc

%% Add path to necessary functions
addpath('../core','../misc')

%% Load data
dataFile = '../datasets/Array_WaterSteelLayers_63Angles.mat';
load(dataFile, 'data','txDelay','steeringAngle','fs','fLow','fHigh',...
    'arrayPitch','cc','thick','waveform2Way')

%% Show raw data examples
tt = (0:(size(data,1)-1))/fs;
xx = (0:(size(data,2)-1))*arrayPitch;
figure
a1 = subplot(1,3,1);
imagesc(xx*1e3,tt*1e6,abs(hilbert(data(:,:,1))))
title({'Raw data envelope,','rightmost angle'})
a2 = subplot(1,3,2);
imagesc(xx*1e3,tt*1e6,abs(hilbert(data(:,:,32))))
title({'Raw data envelope,','center angle'})
a3 = subplot(1,3,3);
imagesc(xx*1e3,tt*1e6,abs(hilbert(data(:,:,63))))
title({'Raw data envelope,','leftmost angle'})
xlabel([a1 a2 a3],'x [mm]')
ylabel([a1 a2 a3],'t [us]')

%% Skip water layer, focus in steel layer
[im,xIm,zIm] = array_psm(data,fs,txDelay,cc,thick,fLow,fHigh,...
    arrayPitch,'wf',waveform2Way,'skipLayers',[1]);

%% Plot focused image of steel layer
figure
imagesc(xIm*1e3,zIm{2}*1e3,abs(im{2}))
xlabel('x [mm]')
ylabel('z [mm]')
title({'Focused image of steel layer,', 'based on all transmitted waves'})
