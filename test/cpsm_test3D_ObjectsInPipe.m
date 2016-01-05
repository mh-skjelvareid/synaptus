%% Close / clear
close all
clear all

%% Load data
load('../datasets/CylScan3D_ObjectsInPipe','ptpz','fs','phiStep','zStep',...
    'tDelay','r0','fLow','fHigh','cc');

%% Gate out reflection from pipe surface (creates artefacts)
rRaw = r0 + (tDelay+(1:size(ptpz,1))/fs) * (cc/2);
rIndRaw = (rRaw <= 0.147);
ptpz = ptpz(rIndRaw,:,:);

%% Processing parameters
rStep = (cc/2)/(fHigh-fLow);
rStart = r0 + (tDelay*(cc/2));
rEnd = rStart + (size(ptpz,1)/fs)*(cc/2);

transFunc = 'haun';     % Transfer function, options 'haun', 'gardner', 'exact'

%% Focus
tic;
[im,phiIm,rIm,zIm] = cpsm(ptpz,fs,tDelay,cc,fLow,fHigh,[phiStep zStep],r0,...
    rStep,rStart,rEnd,'transFunc',transFunc);
disp(['Total processing time: ' num2str(toc) ' seconds.'])

%% Calculate raw signal envelope
envRaw = zeros(size(ptpz));
for ii = 1:size(ptpz,3)
    envRaw(:,:,ii) = abs(hilbert(ptpz(:,:,ii)));
end

%% Make C-scans (maps of maximum envelope)
cScanRaw = squeeze(max(abs(envRaw)));
cScanIm = squeeze(max(abs(im)));

%% Plot results
figure
imagesc(zIm*1e3,phiIm*(180/pi),logImage(cScanRaw));
colorbar
caxis([-30 0])
xlabel('z [mm]')
ylabel('phi [deg]')
title('Raw image C-scan (maximum envelope amplitude, dB)')

figure
imagesc(zIm*1e3,phiIm*(180/pi),logImage(cScanIm));
colorbar
caxis([-30 0])
xlabel('z [mm]')
ylabel('phi [deg]')
title('Focused image C-scan (maximum envelope amplitude, dB)')
