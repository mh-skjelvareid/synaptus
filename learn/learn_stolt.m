%%
% LEARN_STOLT - a simplified version of the Stolt / omega-k algorithm
%
%   This version of the algorithm is included in the toolbox to illustrate the
%   main concepts of the algorithm, without the "code clutter" that arises from
%   making the algorithm as general and efficient as possible. The algorithm
%   is simplified by only considering 2D data from a single-layer medium. Plots
%   of data in time/space-domain and Fourier domain are included to illustrate
%   the process step by step.
%
%   Feel free to change the code, play around and see what happens.
%
%   LICENSE DISCLAIMER:
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   2016-01-03  MHS - Initial version
%
%   Copyright (C) 2016  Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com

close all
clear all

%% Load Data
load('../datasets/LineScan2D_WireTargets.mat', ...
    'ptx',...           % Ultrasound data
    'fs',...            % Sampling frequency
    'xStep',...         % Spatial step size
    'tDelay',...        % Time delay between pulse transmission and measurement
    'cc');              % Wave velocity
[nT,nX] = size(ptx);    % Size of data matrix
fLow = 0.4e6;           % Lower cutoff freq., transducer band
fHigh = 2.5e6;          % Upper cutoff freq., transducer band

%% Static variables related to Z axis
zStart = tDelay*(cc/2);                     % Start of image z range
zEnd = zStart + (nT/fs)*(cc/2);             % End of image z range
dz = (cc/2)/(fHigh-fLow);                   % Z resolution
nZ = length(zStart:dz:zEnd);                % Number of z lines in focused im.

%% Make time/space coodinate variables
tt = tDelay + (0:(nT-1))/fs;                % Time
xx = (0:(nX-1))*xStep;                      % X
zz = zStart + (0:(nZ-1))*dz;                % Z

%% Make spectrum coordinate vectors
nFFTt = 2^nextpow2(nT);                     % FFT size, time
nFFTx = 2^nextpow2(nX);                     % FFT size, x
nFFTz = 2^nextpow2(nT*(fHigh-fLow)/(fs/2));

omega = ((0:(nFFTt-1)) - floor(nFFTt/2))*((2*pi*fs)/nFFTt);     % Freq. vactor
omega = ifftshift(omega);                                       % Shift as fft
omegaBandInd = (omega <= -2*pi*fLow) & (omega >= -2*pi*fHigh);  % Transd. band
omegaBand = omega(omegaBandInd);

kxs = (2*pi)/xStep;                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx); % X-axis wave number vector
kx = ifftshift(kx);                                % Shift as fft

kzs = (2*pi)/dz;                                        % Sampling kz
kz = (2*pi*fLow)/(cc/2) + (0:(nFFTz-1))*(kzs/nFFTz);    % kz vector
kz = -kz(end:-1:1);                     % Correspond to neg. omega

%% Plot raw data
figure
imagesc(xx*1e3,tt*1e6,abs(hilbert(ptx)))
xlabel('x [mm]')
ylabel('tt [us]')
title('Raw data envelope')

%% Fourier transform
Pokx = fft2(ptx,nFFTt,nFFTx);

figure
imagesc(fftshift(kx),fftshift(omega),logImage(fftshift(Pokx)))
caxis([-70 0])
xlabel('kx')
ylabel('omega')
title('2D spectrum')

%% Cut out part of spectrum corresponding to transducer freq. band, neg. omega
Pokx = Pokx(omegaBandInd,:);

figure
imagesc(fftshift(kx),omegaBand,logImage(fftshift(Pokx,2)))
caxis([-70 0])
xlabel('kx')
ylabel('omega')
title('2D spectrum, transducer frequency band only')

%% Resample to higher resolution (improves accuracy in Stolt interpol.)
upSamp = 4;
nOmegaBand = length(omegaBand);
dOmega = ((2*pi*fs)/nFFTt);
Pokx = fft(ifft(Pokx),upSamp*nOmegaBand);
tmp = (0:(nOmegaBand*upSamp-1))*(dOmega/upSamp);
omegaBand = -tmp(end:-1:1) + omegaBand(end);

%% Make 2D spectrum coordinate grids
[OMEGA,KX] = ndgrid(omegaBand,kx);          % 2D grids, omega & kx

%% Time shift according to tDelay (back to t=0)
Pokx = Pokx.*exp(-1i*OMEGA*tDelay);

%% Calculate new spectrum coordinates
[KZ_st,KX_st] = ndgrid(kz,kx);
KK2_st = (KZ_st.^2 + KX_st.^2);                     % KK^2
KK_st = -sqrt(KK2_st .* (KK2_st > 0));              % KK, correspond to neg. om.

Akzkx = 1./(1 + (KX_st.^2)./KZ_st.^2);              % Scale factor
Akzkx(isnan(Akzkx)) = 1;

OMEGA_st = (cc/2)*KK_st;                % New omega coord. for interplation

%% Focus through spectrum resampling (this is where the magic happens)
Pkzkx = complex(zeros(nFFTz,nFFTx));
for jj = 1:nFFTx
     Pkzkx(:,jj) = interp1(omegaBand,Pokx(:,jj),OMEGA_st(:,jj));
end

% Values out af range in interpolation are set to NaN. Change to zero.
nonValid = (isnan(Pkzkx)) | (OMEGA_st < omegaBand(1)) | ...
    (OMEGA_st > omegaBand(end));
Pkzkx(nonValid) = 0;

% Amplitude scaling (consequence of variable change from omega to kz)
Pkzkx = Pkzkx.*Akzkx;

%% Shift back to original start of measurement
Pkzkx = Pkzkx.*exp(1i*KZ_st*zStart);

figure
imagesc(fftshift(kx),kz,logImage(fftshift(Pkzkx,2)))
caxis([-70 0])
xlabel('kx')
ylabel('kz')
title('Spectrum after Stolt resampling')

%% Inverse Fourier transform
im = ifft2(Pkzkx);
im = im(1:nZ,1:nX);

figure
imagesc(xx*1e3,zz*1e3,abs(im))
xlabel(['x [mm]'])
ylabel(['z [mm]'])
title('Focused image')
