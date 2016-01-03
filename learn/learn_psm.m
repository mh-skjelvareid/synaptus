%%
% LEARN_PSM - a simplified version of the PSM algorithm
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

omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);    % Freq. vactor
omega = ifftshift(omega);                                       % Shift as fft
omegaBand = (omega <= -2*pi*fLow) & (omega >= -2*pi*fHigh);     % Transd. band

kxs = (2*pi)/xStep;                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx); % X-axis wave number vector
kx = ifftshift(kx);                                % Shift as fft

[OMEGA,KX] = ndgrid(omega(omegaBand),kx);          % 2D grids, omega & kx

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
Pokx = Pokx(omegaBand,:);

figure
imagesc(fftshift(kx),omega(omegaBand),logImage(fftshift(Pokx,2)))
caxis([-70 0])
xlabel('kx')
ylabel('omega')
title('2D spectrum, transducer frequency band only')

%% Time shift according to tDelay (back to t=0)
Pokx = Pokx.*exp(-1i*OMEGA*tDelay);

%% Calculate phase shift factor
KZ2 = (2/cc)^2*OMEGA.^2 - KX.^2;
realWaveIndex = (KZ2 >= 0);                 % Index of real kz
KZ = sqrt(KZ2.*realWaveIndex);              % Calculate kz
Pokx = Pokx.*realWaveIndex;                 % Mask out evanescent waves
phaseShift = exp(-1i*KZ*dz);                % Phase shift for each z step

%% Focus
% Phase shift to beginning of original measurement
PokxS = Pokx.*exp(-1i*KZ*(zStart-dz));
figure
imagesc(abs(ifft2(Pokx)))
title('Wavefield phase shifted to beginning of original measurement')

% Focus line by line
im = complex(zeros(nZ,nX));
for ii = 1:nZ
    PokxS = PokxS.*phaseShift;              % Phase shift dz for each line
    tmp = ifft(sum(PokxS));                 % Inverse transform
    im(ii,:) = tmp(1,1:nX);                 % Crop to original dimensions

    % Show example plot of wavefield focused at middle of z range
    if (ii == round(nZ/2))
        figure
        imagesc(abs(ifft2(PokxS)))
        title(['Wavefield focused to middle of image z range ' ...
            '(note wrapping in z dim.)'])
    end
end

figure
imagesc(xx*1e3,zz*1e3,abs(im))
xlabel('x [mm]')
ylabel('z [mm]')
title('Focused image')
