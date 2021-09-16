function varargout = stolt_mig_2d(ptx,fs,tDelay,cc,fLow,fHigh,xStep,varargin)
% MULOK - Multi-Layer Omega-K algorithm for synthetic aperture ultrasound
%
%   Usage:
%   ...
%
%   Input parameters:
%   ptx         -   2D pulse-echo data
%   tDelay      -   Time delay between pulse transmission and measurement [s]
%   cc          -   Vector of comp. wave velocities for each layer [m/s]
%   thick       -   Vector of thickness for each layer [m/s]
%                   (set last boundary equal to end of ROI)
%   fLow        -   Lower limit of transducer frequency band [Hz]
%   fHigh       -   Higher limit of transducer frequency band [Hz]
%   xStep       -   spatial step size(s) (distance between measurements) [m]
%
%   Optional parameter-value input pairs:
%   tFftMult    -   multiplier for t-axis FFT size. Default: 1
%   xFftMult    -   multiplier for x-axis FFT size. Default: 1
%   yFftMult    -   multiplier for y-axis FFT size. Default: 1
%   upSamp      -   omega upsamp. factor (preproc. before interp.) Default: 4
%                   Used only for linear interpolation
%
%   Output parameters (varargout):
%   im          -   Reconstructed image as cell array - one cell per layer
%   xIm         -   x positions for pixels in im
%
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
%   2021-09-14  MHS
%
%   Copyright (C) 2021  Martin H. Skjelvareid


%% Parse optional input parameters
% p = inputParser;                        % Input parsing object
% p.addParameter('xFftMult',1);           % Multiplier for FFT size in x dir.
% p.addParameter('yFftMult',1);           % Multiplier for FFT size in y dir.
% p.addParameter('tFftMult',1);           % Multiplier for FFT size in t dir.
% p.addParameter('zFftMult',1)            % Multiplier for FFT size in z dir.
% p.addParameter('upSamp',4)              % Multiplier for interpolation FFT size
% p.addParameter('hh',1);                 % Impulse response / waveform
% p.addParameter('xStart',0)              % First x value of scan
% p.addParameter('yStart',0)              % First x value of scan
% validator = @(str) any(strcmp(str,{'linear','chirpz'}));
% p.addParameter('interpol','linear',validator);     % Interpolation method
% p.addParameter('fc',mean([fHigh,fLow]));           % Center frequency
% p.parse(varargin{:});                               % Parse param.-val. pairs
% param = p.Results;                                  % Store results in "param"
% clear p

%% Get dataset dimension
[nT,nX] = size(ptx);

%% Calculate dependent variables
dzl = (cc/2)./(fHigh-fLow);                     % Z resolution for each layer
zOffset = tDelay*(cc(1)/2);                     % Z offset of meas. window

%% Set FFT sizes
% nFFTx = param.xFftMult * 2^nextpow2(nX);        % FFT size in x dimension
% nFFTy = param.yFftMult * 2^nextpow2(nY);        % FFT size in y dimension
% nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.hh)-1);    % FFT size time
nFFTx = 2^nextpow2(nX);        % FFT size in x dimension
nFFTt = 2^nextpow2(nT+length(param.hh)-1);    % FFT size time

%% Create omega and k vectors
dOmega = ((2*pi*fs)/nFFTt);                             % Omega step size
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*dOmega;       % Omega vector
omega = ifftshift(omega);                               % Shift as fft output

kxs = (2*pi)/sss(1);                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);

%% Calculate bandpass mask for omega
omegaBandIndex = (omega >= (2*pi*fLow)) & (omega <= (2*pi*fHigh));
nOmegaBand = nnz(omegaBandIndex);
fBand = fs*(nOmegaBand/nFFTt);

%% Fourier transform along time dimension
Poxy = fft(ptxy,nFFTt,1);                               % P(omega,x,y)
clear ptxy                                              % Free up space

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

