close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----        SYNAPTUS TOOLBOX TEST: cpsm_test3D             -----')
disp('-----    Cylindrical phase shift migration, 3D dataset      -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core'),fullfile(toolboxPath,'misc'));

%% Load data and convert to double precision
disp('Loading data')
load(fullfile(toolboxPath,'datasets','CylScan3D_LeadSpheres.mat'),'ptpz','fs','phiStep','zStep',...
    'tDelay','r0','fLow','fHigh','cc');
ptpz = double(ptpz);

%% Set parameters for CPSM processing, range and transfer function
rStep = (cc/2)/(fHigh-fLow);
rStart = r0 + (tDelay*(cc/2));
rEnd = rStart + (size(ptpz,1)/fs)*(cc/2);

transFunc = 'haun';     % Transfer function, options 'haun', 'gardner', 'exact'

%% Focus image
disp('Processing data')
tic
[im,phiIm,rIm,zIm] = cpsm(ptpz,fs,tDelay,cc,fLow,fHigh,[phiStep zStep],...
    r0,rStep,rStart,rEnd,'transFunc',transFunc);
toc 

%% Raw data envelope
envRaw = zeros(size(ptpz));
for ii = 1:size(ptpz,3)
    envRaw(:,:,ii) = abs(hilbert(ptpz(:,:,ii)));
end

%% Plot Raw data C-scan
figure
imagesc(zIm*1e3,phiIm*(180/pi),squeeze(max(envRaw)))
xlabel('z [mm]')
ylabel('phi [deg]')
title('Raw image C-scan (max. envelope)')

%% Plot focused image C-scan
figure
imagesc(zIm*1e3,phiIm*(180/pi),abs(squeeze(max(im))));
xlabel('z [mm]')
ylabel('phi [deg]')
title('Focused image C-scan (max. envelope)')

%% Add blank line (nicer formatting for test text output).
disp(' ')