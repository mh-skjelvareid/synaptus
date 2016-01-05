function varargout = mulok(ptxy,fs,tDelay,cc,thick,fLow,fHigh,sss,varargin)
% MULOK - Multi-Layer Omega-K algorithm for synthetic aperture ultrasound
%
%   Usage:
%   [im,...] = mulok(ptxy,fs,tDelay,cc,thick,fLow,fHigh,xStep,yStep,...)
%
%   Input parameters:
%   ptxy        -   2D or 3D ultrasonic data, dimensions (t,x) or (t,x,y)
%   fs          -   Sampling frequency [Hz]
%   tDelay      -   Time delay between pulse transmission and measurement [s]
%   cc          -   Vector of comp. wave velocities for each layer [m/s]
%   thick       -   Vector of thickness for each layer [m/s]
%                   (set last boundary equal to end of ROI)
%   fLow        -   Lower limit of transducer frequency band [Hz]
%   fHigh       -   Higher limit of transducer frequency band [Hz]
%   sss         -   spatial step size(s) (distance between measurements) [m]
%                  2D data: [xStep] - 3D data: [xStep,yStep]
%
%   Optional parameter-value input pairs:
%   interpol    -   Interpolation method for Stolt resampling of spectrum
%                   'linear' - linear resampling (see interp1)
%                   'chirpz' - faster resampling based on linear approximation
%                              of resampling coordinates and Chirp-Z transform
%   fc          -   Transducer center frequency [Hz]. Used in Chirp-Z trans.
%                   Default: mean([fLow,fHigh])
%   tFftMult    -   multiplier for t-axis FFT size. Default: 1
%   xFftMult    -   multiplier for x-axis FFT size. Default: 1
%   yFftMult    -   multiplier for y-axis FFT size. Default: 1
%   zFftMult    -   multiplier for z-axis FFT size. Default: 1
%   upSamp      -   omega upsamp. factor (preproc. before interp.) Default: 4
%                   Used only for linear interpolation
%   hh          -   Transd. imp. res., for matched filt. Default: 1 (no filt.)
%   xStart      -   First x value in "xIm" output [m]
%   yStart      -   First y value "yIm" output [m]
%
%   Output parameters (varargout):
%   im          -   Reconstructed image as cell array - one cell per layer
%   xIm         -   x positions for pixels in im
%   yIm         -   y positions for pixels in im
%   zIm         -   z positions for pixels in im
%                   (cell array, different z res. for each layer)
%
%   Note: The z resolution of the output image for layer l is set to
%       dz = (cc(l)/2)/fHigh.
%   This is done to minimize the required number of z depths for each
%   layer. If this resolution is too coarse, try resampling the image to a
%   higher resolution (see resample, interp)
%
%   The derivation of the MULOK algorithm is described in section 4.2 in
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
%   2010-04-12  MHS - Initial version
%   2015-12-15  MHS - Adaptation for both 2D and 3D data
%
%   Copyright (C) 2016  Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com

%    A note on code notation:
%    ------------------------
%    What Stolt migration does is essentially to transform the sampled
%    wavefield p(t,z=0,x,y) to p(t=0,z,x,y). To understand how this is done,
%    consult f.ex. Margrave (2003) - "Numerical Methods in Exploration
%    Seismology".
%
%    In the process of doing this transformation, many versions of the
%    wavefield are used. A short explanation of the code notation is
%    included here to make reading of the code easier.
%
%    Code variable            Represents wave field
%    ptxy                    p(t,x,y)
%    Poxy                    P(omega,x,y)
%    Pokxky                    P(omega,k_x,ky)
%    Pkzkxky                    P(k_z,k_x,k_y)
%    pzxy                    p(z,x,y) (t=0)

%% Parse optional input parameters
p = inputParser;                        % Input parsing object
p.addParamValue('xFftMult',1);          % Multiplier for FFT size in x dir.
p.addParamValue('yFftMult',1);          % Multiplier for FFT size in x dir.
p.addParamValue('tFftMult',1);          % Multiplier for FFT size in t dir.
p.addParamValue('zFftMult',1)           % Multiplier for FFT size in z dir.
p.addParamValue('upSamp',4)             % Multiplier for interpolation FFT size
p.addParamValue('hh',1);                % Impulse response
p.addParamValue('xStart',0)             % First x value of scan
p.addParamValue('yStart',0)             % First x value of scan
validator = @(str) any(strcmp(str,{'linear','chirpz'}));
p.addParamValue('interpol','linear',validator);     % Interpolation method
p.addParamValue('fc',mean([fHigh,fLow]));           % Center frequency
p.parse(varargin{:});                               % Parse param.-val. pairs
param = p.Results;                                  % Store results in "param"
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
nL = length(thick);                                % Number of layers
zIF = cumsum([0; thick(:)]);                    % z position of each interface

dzl = (cc/2)./(fHigh-fLow);                     % Z resolution for each layer
zOffset = tDelay*(cc(1)/2);                     % Z offset of meas. window

%% Set FFT sizes
nFFTx = param.xFftMult * 2^nextpow2(nX);        % FFT size in x dimension
nFFTy = param.yFftMult * 2^nextpow2(nY);        % FFT size in y dimension
nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.hh)-1);    % FFT size time

%% Create omega and k vectors
dOmega = ((2*pi*fs)/nFFTt);                             % Omega step size
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*dOmega;       % Omega vector
omega = ifftshift(omega);                               % Shift as fft output

kxs = (2*pi)/sss(1);                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);

if dataIs3D
    kys = (2*pi)/sss(2);                                % Samp. wavenum., y dir.
    ky = ((0:(nFFTy-1)) - floor(nFFTy/2))*(kys/nFFTy);  % X-ax. wave num. vec.
    ky = ifftshift(ky);
end

%% Calculate bandpass mask for omega
omegaBandIndex = (omega >= -(2*pi*fHigh)) & (omega <= -(2*pi*fLow));
nOmegaBand = nnz(omegaBandIndex);
fBand = fs*(nOmegaBand/nFFTt);

%% Fourier transform along time dimension
Poxy = fft(ptxy,nFFTt,1);                               % P(omega,x,y)
clear ptxy                                              % Free up space

%% Matched filtering
if param.hh ~= 1
    HH = fft(param.hh,nFFTt,1);                         % Impulse resp. spectrum
    Poxy = Poxy.*repmat(conj(HH),[1 nX nY]);            % Matched filt.
end

%% Cut out frequency band
Poxy = Poxy(omegaBandIndex,:,:);
omegaBand = omega(omegaBandIndex);

%% Interpolate band to higher resolution (improves accuracy in Stolt interpol.)
if strcmp(param.interpol,'linear') && param.upSamp ~= 1
    Poxy = fft(ifft(Poxy),param.upSamp*nOmegaBand);
    tmp = (0:(nOmegaBand*param.upSamp-1))*(dOmega/param.upSamp);
    omegaBand = -tmp(end:-1:1) + omegaBand(end);
end

%% tDelay compensation
Poxy = Poxy.*repmat(exp(-1i*omegaBand(:)*tDelay),[1 nX nY]);

%% Fourier transform in x and y direction
Pokxky = fft(Poxy,nFFTx,2);
if dataIs3D
    Pokxky = fft(Pokxky,nFFTy,3);
end
clear Poxy

%% Create "PSM grids" used for extrapolation between interfaces
if dataIs3D
    [OMEGA_psm,KX_psm,KY_psm] = ndgrid(omegaBand,kx,ky);
else
    [OMEGA_psm,KX_psm] = ndgrid(omegaBand,kx);
end

%% Extraplolate wavefield to each interface, and image using Stolt transf.
im = cell(nL,1);        % Preallocate cell structure for each layer image
zIm = cell(nL,1);       % Preallocate cell structure for corresponding z axis

for ii = 1:nL
    disp(['Processing layer ' num2str(ii) ' of ' num2str(nL)]);

    if strcmp(param.interpol,'linear')
        %%%%%%%%  REGULAR STOLT RESAMPLING (LINEAR INTERPOLATION) %%%%%%%%
        nFFTz = param.zFftMult * 2^nextpow2(nT*(fHigh-fLow)/(fs/2));
        kzs = (2*pi)/dzl(ii);                   % Sampling kz
        kz = (2*pi*fLow)/(cc(ii)/2) + (0:(nFFTz-1))*(kzs/nFFTz);    % kz vector
        kz = -kz(end:-1:1);                     % Correspond to neg. omega

        % "Stolt grids" (correspond to low resolution P(kz,kx,ky) wavefield)
        if dataIs3D
            [KZ_st,KX_st,KY_st] = ndgrid(kz,kx,ky);
            KK_st = (KZ_st.^2 + KX_st.^2 + KY_st.^2);           % KK^2
            Akzkxky = 1./(1 + (KX_st.^2 + KY_st.^2)./KZ_st.^2); % Scale factor
            clear KX_st KY_st
        else
            [KZ_st,KX_st] = ndgrid(kz,kx);
            KK_st = (KZ_st.^2 + KX_st.^2);                      % KK^2
            Akzkxky = 1./(1 + (KX_st.^2)./KZ_st.^2);            % Scale factor
            clear KX_st
        end
        KK_st = -sqrt(KK_st .* (KK_st > 0));        % Calc. KK with square root
        Akzkxky(isnan(Akzkxky)) = 0;                % Remove single NaN point.
                                                    % TODO: Check NaN - why?

        % Calc. omega values to be interpolated for
        OMEGA_st = (cc(ii)/2)*KK_st;

        % Interpolate for each (kx,ky)
        if dataIs3D
            Pkzkxky = complex(zeros(nFFTz,nFFTx,nFFTy));
            for jj = 1:nFFTx
                for kk = 1:nFFTy
                     Pkzkxky(:,jj,kk) = interp1(omegaBand(:),Pokxky(:,jj,kk),...
                        OMEGA_st(:,jj,kk));
                end
            end
        else
            Pkzkxky = complex(zeros(nFFTz,nFFTx));
            for jj = 1:nFFTx
                 Pkzkxky(:,jj) = interp1(omegaBand(:),Pokxky(:,jj),...
                    OMEGA_st(:,jj));
            end
        end

        % Values out af range in interpolation are set to NaN. Change to zero.
        nonValid = (isnan(Pkzkxky)) | (OMEGA_st < omegaBand(1)) | ...
            (OMEGA_st > omegaBand(end));
        Pkzkxky(nonValid) = 0;

        % Amplitude scaling (because of variable change from omega to kz)
        Pkzkxky = Pkzkxky.*Akzkxky;

        % If first layer, shift according to measurement z offset
        if ii == 1
            Pkzkxky = Pkzkxky.*exp(1i*KZ_st*zOffset);
        end

        % Inverse transform
        pzxy = ifftn(Pkzkxky);

    else
        %%%%%%%%  CHIRP-Z APPROXIMATION TO STOLT INTERPOLATION %%%%%%%%
        Ptkxky = ifft(Pokxky,[],1);
        kzc = (2*pi*param.fc)/(cc(ii)/2);

        % Interpolate for each kx
        if dataIs3D
            Pkzkxky = complex(zeros(nOmegaBand,nFFTx,nFFTy));   % Preallocate
            for jj = 1:nFFTx
                for kk = 1:nFFTy
                    % Calc K (scale factor), A (spec. offset) and W (spec. step)
                    K = sqrt(1 + (kx(jj)^2+ky(kk)^2)/kzc^2);
                    A = -2*pi*((fLow/fBand)*(1/K-1) + (param.fc/fBand)*(K-1/K));
                    W = 2*pi*(1/K)*(fs/fBand)/nFFTt;
                    Pkzkxky(:,jj,kk) = qczt(Ptkxky(:,jj,kk),nOmegaBand,W,A);
                end
            end
        else % Data is 2D
            Pkzkxky = complex(zeros(nOmegaBand,nFFTx));         % Preallocate
            for jj = 1:nFFTx
                % Calc K (scale factor), A (spec. offset) and W (spec. step)
                K = sqrt(1 + (kx(jj)/kzc)^2);
                A = -2*pi*((fLow/fBand)*(1/K-1) + (param.fc/fBand)*(K-1/K));
                W = 2*pi*(1/K)*(fs/fBand)/nFFTt;
                Pkzkxky(:,jj) = qczt(Ptkxky(:,jj),nOmegaBand,W,A);
            end
        end % if dataIs3D

        % If first layer, shift according to measurement z offset
        if ii == 1
            kz = omegaBand/(cc(ii)/2);
            if dataIs3D
                Pkzkxky = Pkzkxky.*exp(-1i*repmat(kz,1,nFFTx,nFFTy)*zOffset);
            else
                Pkzkxky = Pkzkxky.*exp(-1i*repmat(kz,1,nFFTx)*zOffset);
            end
        end

        % Inverse transform
        pzxy = ifftn(flip(Pkzkxky,1));

    end

    % Calculate number of image planes in each layer, calc. z axis coordinates
    if ii == 1
        nPlanesZ = ceil((thick(ii)-zOffset)/dzl(ii)) + 1;
        zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zOffset;
    else
        nPlanesZ = ceil(thick(ii)/dzl(ii)) + 1;
        zIm{ii} = (0:(nPlanesZ-1))*dzl(ii) + zIF(ii);
    end

    % Cut out part of pzxy corresponding to current layer
    if dataIs3D
        im{ii} = pzxy(1:nPlanesZ,1:nX,1:nY);
    else
        im{ii} = pzxy(1:nPlanesZ,1:nX);
    end

    % Migrate to next layer (if not last layer)
    if ii < nL
        KK_psm = OMEGA_psm/(cc(ii)/2);
        if dataIs3D
            KZ2_psm = (KK_psm.^2 - KX_psm.^2 - KY_psm.^2);
        else
            KZ2_psm = (KK_psm.^2 - KX_psm.^2);
        end
        KZ_psm = sqrt(KZ2_psm .* (KZ2_psm > 0));
        Pokxky = Pokxky .* exp(-1i*KZ_psm*thick(ii));
    end
end

%% Assign output
varargout{1} = im;                                      % Focused image
varargout{2} = param.xStart + (0:(nX-1))*sss(1);        % X-axis pixel positions
if dataIs3D
    varargout{3} = param.yStart + (0:(nY-1))*sss(2);    % Y-axis pixel positions
    varargout{4} = zIm;                                 % Z-axis pixel positions
else
    varargout{3} = zIm;                                 % Z-axis pixel positions
end
