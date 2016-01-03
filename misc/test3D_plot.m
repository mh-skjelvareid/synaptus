%% Remove water layer (empty)
im = im(2:3);
zIm = zIm(2:3);

%% Plotting parameters
[nT,nX,nY] = size(ptxy);		% Size of B-scan matrix
tPlot = (0:(nT-1))/fs + tDelay;	% Time vector for y axis

sliceInd1 = 13;                 % X position for scatterers in 1. layer
sliceInd2 = 24;                 % X position for scatterers in 2. layer

tInd1  = (tPlot >= 88e-6) & (tPlot <= 106e-6);      % Time gating, 1. layer
tInd2  = (tPlot >= 110e-6) & (tPlot <= 122.5e-6);   % Time gating, 2. layer

zInd1  = (zIm{1} >= 0.068) & (zIm{1} <= 0.092);     % Depth gating, 1. layer
zInd2  = (zIm{2} >= 0.1) & (zIm{2} <= 0.14);        % Depth gating, 1. layer

%% Calculate envelopes and C-scans
envRaw = zeros(size(ptxy));
for ii = 1:size(ptxy,3)
    envRaw(:,:,ii) = abs(hilbert(ptxy(:,:,ii)));
end
envIm1 = abs(im{1});
envIm2 = abs(im{2});

cScanRaw1 = squeeze(max(envRaw(tInd1,:,:)));
cScanRaw2 = squeeze(max(envRaw(tInd2,:,:)));
cScanIm1 = squeeze(max(envIm1(zInd1,:,:)));
cScanIm2 = squeeze(max(envIm2(zInd2,:,:)));

%% Plot
% Raw data plots
figure
subplot(4,2,1)
imagesc(xIm*1e3,tPlot(tInd1)*1e6,envRaw(tInd1,:,sliceInd1))
title('Raw data slice, 1. layer')
ylabel('t [us]')

subplot(4,2,3)
imagesc(xIm*1e3,yIm*1e3,cScanRaw1.')
title('Raw data C-scan, 1. layer')
ylabel('y [mm]')

subplot(4,2,5)
imagesc(xIm*1e3,tPlot(tInd2)*1e6,envRaw(tInd2,:,sliceInd2))
title('Raw data slice, 2. layer')
ylabel('t [us]')

subplot(4,2,7)
imagesc(xIm*1e3,yIm*1e3,cScanRaw2.')
title('Raw data C-scan, 1. layer')
ylabel('y [mm]')
xlabel('x [mm]')

% Focused image plots
subplot(4,2,2)
imagesc(xIm*1e3,zIm{1}(zInd1)*1e3,envIm1(zInd1,:,sliceInd1))
title('Focused image slice, 1. layer')
ylabel('z [mm]')

subplot(4,2,4)
imagesc(xIm*1e3,yIm*1e3,cScanIm1.')
title('Focused image C-scan, 1. layer')
ylabel('y [mm]')

subplot(4,2,6)
imagesc(xIm*1e3,zIm{2}(zInd2)*1e3,envIm2(zInd2,:,sliceInd2))
title('Focused image slice, 2. layer')
ylabel('z [mm]')

subplot(4,2,8)
imagesc(xIm*1e3,yIm*1e3,cScanIm2.')
title('Focused image C-scan, 2. layer')
ylabel('y [mm]')
xlabel('x [mm]')
