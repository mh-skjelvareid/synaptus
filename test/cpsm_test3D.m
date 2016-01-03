close all
clear all

%% Load data, convert from single to double precision
load('../datasets/CylScan3D_LeadSpheres','ptpz','fs','phiStep','zStep',...
    'tDelay','r0','fLow','fHigh','cc');
ptpz = double(ptpz);

%% Set parameters for CPSM processing, range and transfer function
rStep = (cc/2)/(fHigh-fLow);
rStart = r0 + (tDelay*(cc/2));
rEnd = rStart + (size(ptpz,1)/fs)*(cc/2);

transFunc = 'haun';     % Transfer function, options 'haun', 'gardner', 'exact'

%% Focus image
tic;
[im,phiIm,rIm,zIm] = cpsm(ptpz,fs,tDelay,cc,fLow,fHigh,[phiStep zStep],...
    r0,rStep,rStart,rEnd,'transFunc',transFunc);
disp(['Total processing time: ' num2str(toc) ' seconds.'])

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
