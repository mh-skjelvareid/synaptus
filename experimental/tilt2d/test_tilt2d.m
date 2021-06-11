%% Close / clear
close all
clearvars

%% Load data
load('../../datasets/CopperBlock_2degtTilt','fs','ptx','tDelay',...
    'fLow','fHigh','xStep','cc')
[nT,nX] = size(ptx);

%% Parameters
xPlot = (0:(nX-1))*xStep;
tPlot = tDelay + (0:(nT-1))/fs;
zPlot = tPlot(:)*(cc(1)/2);

%% Tilt parameters (could be estimated from data)
pos.x = [0.002, 0.299];
pos.z = [0.0434,0.0549];

aa = (pos.z(2) - pos.z(1))/(pos.x(2) - pos.x(1));
bb = pos.z(1) - aa*pos.x(1);

%% Plot raw data
figure
colormap(gray(256))
imagesc(xPlot*1e3,zPlot*1e3,logImage(ptx))
set(gca,'CLim',[-40 0])
xlabel('X pos. [mm]')
ylabel('Z pos. ref. water [mm]')
title('Original B-scan')
hold on

%% Plot line corresponding to interface
line = aa*xPlot + bb;
plot(xPlot*1e3,line*1e3,'b--','LineWidth',2)
legend('Interface line')

%% Tilt and shift image to surface
[ptxTilt,xStepNew] = tilt2d(ptx,fs,tDelay,cc(1),fLow,fHigh,xStep,aa,bb);

%% Plot tilted image
xPlotNew = (0:(nX-1))*xStepNew;

figure
colormap(gray(256))
imagesc(xPlotNew*1e3,zPlot*1e3,logImage(real(ptxTilt)))
set(gca,'CLim',[-40 0])
xlabel('X pos. [mm]')
ylabel('Z pos. ref. water [mm]')
title('B-scan extrapolated to tilted surface')
hold on

%% Use wavefield at surface as starting point for focusing inside block
tDelay_tilt = 0;
thick_copper = 0.08;
[im,xIm,zIm] = psm(ptxTilt,fs,tDelay_tilt,cc(2),thick_copper,fLow,fHigh,xStepNew);

%% Show focused image
figure
colormap(gray(256))
imagesc(xIm*1e3,zIm{1}*1e3,logImage(im{1}))
set(gca,'CLim',[-40 0])
xlabel('X [mm]')
ylabel('Z [mm]')
title('Focused image of copper block')

