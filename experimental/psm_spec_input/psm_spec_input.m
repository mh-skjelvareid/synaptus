function varargout = psm_spec_input(ptxy,omega,kx,tDelay,cc,thick,varargin)
% psm_spec_input - Phase Shift Migration with spectral input
%
%   Usage:
%   [im,...] = psm_spec_input(ptxy,fs,tDelay,cc,thick,fLow,fHigh,sss,...)
%
%   Input parameters:
%   Pokx     -  MxN matrix of spectral data, dimensions (omega,kx).
%   omega   -   vector, length M, with omega coordinates.
%               Omega values are assumed to be positive, corresponding to 
%               the passband of the transducer / antenna / signal. 
%   kx      -   vector, length N, with kx coordinates.
%               Kx values are usually both positive and negative (symmetry
%               around zero), but this is not a requirement.
%   tDelay  -    Time delay between pulse transmission and measurement [s]
%   cc      -    Vector of compressional wave velocities for each layer [m/s]
%   thick   -    Vector of thickness for each layer [m/s]
%                  (set last boundary equal to end of ROI)
%
%   Output parameters (varargout):
%   im          -    Reconstructed image as cell array - one cell per layer
%   xIm         -    x positions for pixels in im
%   zIm         -    z positions for pixels in im
%                       (cell array, different z res. for each layer)
%
%   The derivation of the PSM algorithm is described in sections 3.5 and 5.1
%   the PhD thesis "Synthetic aperture ultrasound imaging with application to
%   interior pipe inspection" (2012) by Martin H. Skjelvareid, University of
%   Tromsø, Norway.
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
%   2021-09-16  MHS - Initial version
%
%   Copyright (C) 2021  Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com


%% Dimension testing 
[nOmega,nKx] = size(Pokx);
if length(omega) ~= nOmega
    error('Length of omega vector must correspond to number of rows in Pokx')
end
if length(kx) ~= nKx
    error('Length of kx vector must correspond to number of columns in Pokx')
end

%% Calculate dependent variables
nL = length(thick);                                 % Number of layers
zIF = cumsum([0; thick(:)]);                        % z position of interfaces

dzl = (cc/2)./(max(omega)/(2*pi));                  % Z resolution in each layer
zOffset = tDelay*(cc(1)/2);                         % Z offset due to tDelay
nZ = ceil(sum([thick(1)-zOffset; thick(2:end).']./dzl(:)))+nL;  % Num. Z depths

%% Create grids
[OMEGA,KX] = ndgrid(omega,kx);


%% Phase shift migration through each layer
im = cell(nL,1);                        % Preallocate image cell structure
zIm = cell(nL,1);                       % Preallocate z axis cell structure
planeCount = 0;                         % Counter for # focused lines / planes

fprintf('Progress: 0123');

for ii = 1:nL
    % Calculate z-axis wave number KZ
    KZ2 = (2/cc(ii))^2*OMEGA.^2 - KX.^2;

    realWaveIndex = (KZ2 >= 0);                 % Index of real kz
    KZ = sqrt(KZ2.*realWaveIndex);              % Calculate kz
    Pokx = Pokx.*realWaveIndex;             % Mask out evanescent waves
    phaseShift = exp(1i*KZ*dzl(ii));           % Phase shift for each z step

    if ii == 1
        PokxShifted = Pokx.*exp(1i*KZ*zOffset);    % Shift to tDelay depth
        nPlanesZ = ceil((thick(1)-zOffset)/dzl(1))+1;   % Calc. # Z depths
        zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zOffset;   % Calc. Z depths
    else
        PokxShifted = Pokx;                         % Copy to local wavef.
        nPlanesZ = ceil(thick(ii)/dzl(ii)) + 1;         % Calc. # Z depths
        zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zIF(ii);   % Calc. Z depths
    end

    iml = zeros(nPlanesZ,nX,nY);                        % Preallocate loc. image

    for jj = 1:nPlanesZ
        Pkxky_t0 = sum(PokxShifted,1);                % Sum (IFFT at t=0)
        imPlane = ifft(Pkxky_t0,nFFTx,2);               % FFT along x direction
        iml(jj,:) = imPlane(1:nX);                      % Insert line in matrix

        PokxShifted = PokxShifted.*phaseShift;      % Phase shift to next z
        planeCount = planeCount + 1;
        if ~mod(planeCount,25)                          % Progress information
            fprintf('\b\b\b\b%3d%%', round((planeCount*100)/nZ));
        end
    end

    % Copy to cell structure
    im{ii} = iml;

    % Migrate ref. wavefield to next interface
    if ii < nL
        Pokx = Pokx.*exp(1i*KZ*thick(ii));
    end
end

% Print progress = 100%
fprintf('\b\b\b\b%3d%%\n', 100);

%% Assign output
varargout{1} = im;                                      % Focused image
varargout{2} = param.xStart + (0:(nX-1))*sss(1);        % X-axis pixel positions
varargout{3} = zIm;                                     % Z-axis pixel positions
