function varargout = psm(ptxy,fs,tDelay,cc,thick,fLow,fHigh,sss,varargin)
% psm - Phase Shift Migration for synthetic aperture focusing
%
%   Usage:
%   [im,...] = psm(ptxy,fs,tDelay,cc,thick,fLow,fHigh,sss,...)
%
%   Input parameters:
%   ptxy    -    2D or 3D ultrasonic data, dimensions (t,x) or (t,x,y)
%   fs      -    Sampling frequency [Hz]
%   tDelay  -    Time delay between pulse transmission and measurement [s]
%   cc      -    Vector of compressional wave velocities for each layer [m/s]
%   thick   -    Vector of thickness for each layer [m/s]
%                  (set last boundary equal to end of ROI)
%   fLow    -    Lower limit of transducer frequency band [Hz]
%   fHigh   -    Higher limit of transducer frequency band [Hz]
%   sss     -    spatial step size(s) (distance between measurements) [m]
%                   2D data: [xStep] - 3D data: [xStep,yStep]
%
%   Optional parameter-value input pairs:
%   tFftMult    -    multiplier for t-axis FFT size. Default: 1
%   xFftMult    -    multiplier for x-axis FFT size. Default: 1
%   yFftMult    -    multiplier for y-axis FFT size. Default: 1
%   hh          -    Transd. imp. res., for matched filt. Default: 1 (no filt.)
%   xStart      -    First x value in "xIm" output [m]
%   yStart      -    First y value "yIm" output [m]
%
%   Output parameters (varargout):
%   im          -    Reconstructed image as cell array - one cell per layer
%   xIm         -    x positions for pixels in im
%   yIm         -    y positions for pixels in im
%   zIm         -    z positions for pixels in im
%                       (cell array, different z res. for each layer)
%
%   Note: The z resolution of the output image for layer l is set to
%       dz = (cc(l)/2)/fHigh.
%   This is done to minimize the required number of z depths for each
%   layer. If this resolution is too coarse, try resampling the image to a
%   higher resolution (see resample, interp)
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
%   2010-04-27  MHS - Initial version
%   2015-12-15  MHS - Adaptation for both 2D and 3D data
%
%   Copyright (C) 2016  Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com


%    A note on code notation:
%    ------------------------
%    In this implementation of PSM, several versions of the discretely
%    sampled wavefield is used. A short explanation of the code notation is
%    included here to make reading of the code easier.
%
%    Code variable          Represents wave field
%    ptxy                   p(t,x,y)
%    Poxy                   P(omega,x,y)
%    Pokxky                 P(omega,k_x,k_y)
%    Pkxky_t0               P(t=0,k_x,k_y)


%% Parse optional input arguments
p = inputParser;                        % Create input parser
p.addParamValue('xFftMult',1);          % Multiplier for FFT size in x dir.
p.addParamValue('yFftMult',1);          % Multiplier for FFT size in y dir.
p.addParamValue('tFftMult',1);          % Multiplier for FFT size in t dir.
p.addParamValue('hh',1)                 % Transducer impulse response
p.addParamValue('xStart',0)             % First x value of scan
p.addParamValue('yStart',0)             % First y value of scan
p.parse(varargin{:});                   % Parse possible parameter-value pairs
param = p.Results;                      % Transfer results to "param" structure
clear p

%% Dimension testing (2D is a special case)
dataIs3D = (ndims(ptxy) == 3);

% Check that step sizes are correctly given
if dataIs3D && length(sss) == 1
    error('The sss input parameter must have 2 elements for 3D data sets')
end

% Get dataset size
if dataIs3D
    [nT,nX,nY] = size(ptxy);
else
    [nT,nX] = size(ptxy);
    nY = 1;
end

%% Calculate dependent variables
nL = length(thick);                                 % Number of layers
zIF = cumsum([0; thick(:)]);                        % z position of interfaces

dzl = (cc/2)./(fHigh);                              % Z resolution in each layer
zOffset = tDelay*(cc(1)/2);                         % Z offset due to tDelay
nZ = ceil(sum([thick(1)-zOffset; thick(2:end).']./dzl(:)))+nL;  % Num. Z depths

nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.hh)-1);     % FFT size time
nFFTx = param.xFftMult * 2^nextpow2(nX);                        % FFT size x
if dataIs3D
    nFFTy = param.yFftMult * 2^nextpow2(nY);                    % FFT size y
end

omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);    % Freq. vactor
omega = ifftshift(omega);                                       % Shift as fft
omegaBandIndex = (omega >= -(2*pi*fHigh)) & (omega <= -(2*pi*fLow));

kxs = (2*pi)/sss(1);                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);

if dataIs3D
    kys = (2*pi)/sss(2);                                % Samp. wavenum., y dir.
    ky = ((0:(nFFTy-1)) - floor(nFFTy/2))*(kys/nFFTy);  % X-ax. wave num. vec.
    ky = ifftshift(ky);
end

%% Fourier transform along time dimension for both scan and imp.response
Poxy = fft(ptxy,nFFTt,1);
HH = fft(param.hh,nFFTt,1);

%% Matched filtering and phase shift (compensation for tDelay time delay)
matchAndDelayVector = conj(HH).*exp(-1i*omega*tDelay);
Poxy = Poxy.*repmat(matchAndDelayVector,[1 nX nY]);

%% Extract frequency band of interest
Poxy = Poxy(omegaBandIndex,:,:);

%% Fourier transform in x and y direction
Pokxky = fft(Poxy,nFFTx,2);
if dataIs3D
    Pokxky = fft(Pokxky,nFFTy,3);
end

%% Create grids
if dataIs3D
    [OMEGA,KX,KY] = ndgrid(omega(omegaBandIndex),kx,ky);
else
    [OMEGA,KX] = ndgrid(omega(omegaBandIndex),kx);
end

%% Phase shift migration through each layer
im = cell(nL,1);                        % Preallocate image cell structure
zIm = cell(nL,1);                       % Preallocate z axis cell structure
planeCount = 0;                         % Counter for # focused lines / planes

fprintf('Progress: 0123');

for ii = 1:nL
    % Calculate z-axis wave number KZ
    if dataIs3D
        KZ2 = (2/cc(ii))^2*OMEGA.^2 - KX.^2 - KY.^2;
    else
        KZ2 = (2/cc(ii))^2*OMEGA.^2 - KX.^2;
    end

    realWaveIndex = (KZ2 >= 0);                 % Index of real kz
    KZ = sqrt(KZ2.*realWaveIndex);              % Calculate kz
    Pokxky = Pokxky.*realWaveIndex;             % Mask out evanescent waves
    phaseShift = exp(-1i*KZ*dzl(ii));           % Phase shift for each z step

    if ii == 1
        PokxkyShifted = Pokxky.*exp(-1i*KZ*zOffset);    % Shift to tDelay depth
        nPlanesZ = ceil((thick(1)-zOffset)/dzl(1))+1;   % Calc. # Z depths
        zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zOffset;   % Calc. Z depths
    else
        PokxkyShifted = Pokxky;                         % Copy to local wavef.
        nPlanesZ = ceil(thick(ii)/dzl(ii)) + 1;         % Calc. # Z depths
        zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zIF(ii);   % Calc. Z depths
    end

    iml = zeros(nPlanesZ,nX,nY);                        % Preallocate loc. image

    for jj = 1:nPlanesZ
        Pkxky_t0 = sum(PokxkyShifted,1);                % Sum (IFFT at t=0)
        imPlane = ifft(Pkxky_t0,nFFTx,2);               % FFT along x direction
        if dataIs3D
            imPlane = ifft(imPlane,nFFTy,3);            % FFT along y direction
            iml(jj,:,:) = imPlane(1,1:nX,1:nY);         % Insert plane in matrix
        else
            iml(jj,:) = imPlane(1:nX);                  % Insert line in matrix
        end

        PokxkyShifted = PokxkyShifted.*phaseShift;      % Phase shift to next z
        planeCount = planeCount + 1;
        if ~mod(planeCount,25)                          % Progress information
            fprintf('\b\b\b\b%3d%%', round((planeCount*100)/nZ));
        end
    end

    % Copy to cell structure
    im{ii} = iml;

    % Migrate ref. wavefield to next interface
    if ii < nL
        Pokxky = Pokxky.*exp(-1i*KZ*thick(ii));
    end
end

% Print progress = 100%
fprintf('\b\b\b\b%3d%%\n', 100);

%% Assign output
varargout{1} = im;                                      % Focused image
varargout{2} = param.xStart + (0:(nX-1))*sss(1);        % X-axis pixel positions

if dataIs3D
    varargout{3} = param.yStart + (0:(nY-1))*sss(2);    % Y-axis pixel positions
    varargout{4} = zIm;                                 % Z-axis pixel positions
else
    varargout{3} = zIm;                                 % Z-axis pixel positions
end
