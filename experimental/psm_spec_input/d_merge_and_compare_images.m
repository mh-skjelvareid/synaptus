close all
clearvars

%% Paths
% Update these paths if re-processing on your local computer
rawDataPath = 'E:\Synaptus\RadarAntarcticaData\Preprocessed\';
focDataPath = 'E:\Synaptus\RadarAntarcticaData\Preprocessed\Images\';

%% Load raw data
rawDataFiles = dir([rawDataPath '*.mat']);
rawDataFiles = {rawDataFiles.name};
nRawFiles = length(rawDataFiles);

Pox_all = [];
for ii = 1:nRawFiles
    load([rawDataPath rawDataFiles{ii}],'Pox','xStep')
    Pox_all = [Pox_all Pox];
end

%% Transform to time domain
ptx_all = ifft(Pox_all,[],1);

%% Make t and x axes
bw = 60e6;          % Passband bandwidth
tStep = 1/bw;
tPlot = (0:(size(ptx_all,1)-1)) * tStep;
xPlot = (0:(size(ptx_all,2)-1)) * xStep;

%% Plot raw data
% Note: logImage() is found in the 'misc' folder in the synaptus toolbox
figure(1)
colormap(flipud(gray))
imagesc(xPlot,tPlot*1e6,logImage(ptx_all))
caxis([-70 0])
colorbar
xlabel('x position [m]')
ylabel('t [us]')
ylim([0 34])
title('Raw data (merged) [dB]')

%% Export plot to file
% exportgraphics(gca,'RawDataMerged.tif')

%% Load focused image data
focDataFiles = dir([focDataPath '*.mat']);
focDataFiles = {focDataFiles.name};
nFocFiles = length(focDataFiles);

im_all = [];
for ii = 1:nFocFiles
    load([focDataPath focDataFiles{ii}],'im','zIm')
    im_all = [im_all im{1}];
end

%% Plot focused image
figure(2)
colormap(flipud(gray))
imagesc(xPlot,zIm{1},logImage(im_all))
caxis([-80 0])
colorbar
xlabel('x position [m]')
ylabel('z depth [m]')
title('Focused image (merged) [dB]')
ylim([0 4000])

%% Export plot to file
% exportgraphics(gca,'FocusedDataMerged.tif')

%% Save focused image as *.mat file
% xIm = xPlot;
% save('FocusedDataMerged','im','xPlot','zIm')