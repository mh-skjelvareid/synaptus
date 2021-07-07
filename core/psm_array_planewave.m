function varargout = psm_array_planewave(ptxa,fs,txDelay,cc,thick,fLow,fHigh,aPitch,varargin)
% psm - Phase Shift Migration for linear arrays with steered plane waves  
%
%   Usage:
%   [im,...] = psm_array_planewave(ptxa,fs,txDelay,cc,thick,fLow,fHigh,aPitch,...)
%
%   Input parameters:
%   ptxa    -    ultrasonic data, dimensions (t,x,angle). It is assumed
%                that there is no "measurement delay", i.e. the data starts
%                at t=0, when the first pulse is fired. 
%   fs      -    Sampling frequency [Hz]
%   txDelay -    2D array of pulse time delays [s], dimensions (x, angle)
%   cc      -    Vector of compressional wave velocities for each layer [m/s]
%   thick   -    Vector of thickness for each layer [m/s]
%                  (set last boundary equal to end of ROI)
%   fLow    -    Lower limit of transducer frequency band [Hz]
%   fHigh   -    Higher limit of transducer frequency band [Hz]
%   aPitch  -    array pitch (distance between array element centers) [m]
%
%   Optional parameter-value input pairs:
%   tFftMult    -    multiplier for t-axis FFT size. Default: 1
%   xFftMult    -    multiplier for x-axis FFT size. Default: 2
%   zResMult    -    multiplier for final image z resolution. Default: 2
%   wf          -    2-way waveform, for matched filt. Default: 1 (no filt.)
%   skipLayers  -    vector w/indices of layers to skip in processing 
%                   (empty output). Default: [] (all layers included).
%
%   Output parameters (varargout):
%   im          -    Reconstructed image as cell array - one cell per layer
%   xIm         -    x positions for pixels in im
%   zIm         -    z positions for pixels in im
%                       (cell array, different z res. for each layer)
%
%   Note: The default z resolution of the output image of layer L is
%       dz = (1/zResMult)*(cc(L)/2)/(fHigh-fLow).
%   This is done to minimize the number of z depths processed for each
%   layer and speed up the imaging. If this resolution is too coarse, 
%   try resampling the image to a higher resolution (see resample, interp)
%   or increase the zResMult parameter.  
%
%   The derivation of the PSM algorithm is described in sections 3.5 and 5.1
%   the PhD thesis "Synthetic aperture ultrasound imaging with application to
%   interior pipe inspection" (2012) by Martin H. Skjelvareid, University of
%   Tromsø, Norway. However, ahe adaptation to array imaging is not 
%   included in the PhD thesis.
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
%   2011-07-05  MHS - Initial version, based on previous PSM code
%
%   2021  Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com


%    A note on code notation:
%    ------------------------
%    In this implementation of PSM, several versions of the discretely
%    sampled wavefield is used. A short explanation of the code notation is
%    included here to make reading of the code easier.
%
%    Code variable          Represents wave field
%    ptxa                   p(t,x,angle)
%    Poxa                   P(omega,x,angle)
%    Pokxa                  P(omega,k_x,angle)
%    Pox                    P(omega,x)


%% Parse optional input arguments
p = inputParser;                        % Create input parser
p.addParameter('xFftMult',2);           % Multiplier for FFT size in x dir.
p.addParameter('tFftMult',1);           % Multiplier for FFT size in t dir.
p.addParameter('zResMult',2);           % Multiplier for z resolution
p.addParameter('wf',1)                  % 2-way waveform
p.addParameter('skipLayers',[])         % Boolean index of layers to skip
p.parse(varargin{:});                   % Parse possible parameter-value pairs
param = p.Results;                      % Transfer results to "param" structure
clear p

%% Calculate dependent variables
[nT,nX,nA] = size(ptxa);        % # time samples, # array elements, # angles    
nL = length(thick);             % Number of layers
zIF = cumsum([0; thick(:)]);    % z coordinates of interfaces

dzl = (1/param.zResMult)*((cc/2)./(fHigh-fLow));    % Z resolution in each layer
nZ = ceil(thick(:)./dzl(:)) + 1;                    % Num. of Z depths in each layer (with some margin)
zz = cell(nL,1);                                    % Z axes (1 cell per layer)
for ii = 1:nL
    zz{ii} = (0:(nZ(ii)-1))*dzl(ii) + zIF(ii);      % Z coordinates in each layer
end

nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.wf)-1);     % FFT size time
nFFTx = param.xFftMult * 2^nextpow2(nX);                        % FFT size x

omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);    % Omega vector, neg. to pos.
omega = ifftshift(omega);                                       % Shift corr. to fft output
omegaBandIndex = (omega >= 2*pi*fLow) & (omega <= 2*pi*fHigh);

kxs = (2*pi)/aPitch;                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);                                 % Shift corr. to fft. output

%% Estimate steering angle in each layer based on time delays
sAngle = zeros(nA,nL);
for ii = 1:nL
    sAngle(:,ii) = asin(mean(diff(txDelay))*cc(ii)/aPitch).';   % theta = asin(dt*c/dx)
end

%% Fourier transform along time dimension for both scan and waveform
Poxa = fft(ptxa,nFFTt,1);
WF = fft(param.wf,nFFTt,1);

%% Matched filtering
matchedFilter = conj(WF);
Poxa = Poxa.*matchedFilter;       

%% Extract frequency band of interest
Poxa = Poxa(omegaBandIndex,:,:);

%% Fourier transform in x direction
Pokxa = fft(Poxa,nFFTx,2);

%% Create grids
[OMEGA,KX] = ndgrid(omega(omegaBandIndex),kx);

%% Phase shift migration through each layer
im = cell(nL,1);            % Preallocate image cell structure

for ll = 1:nL
    % Progress update
    disp('*************************')
    disp(['Processing layer ' num2str(ll) ' of ' num2str(ll)])
    disp('*************************')
    
    % Calculate z-axis wave number KZ for scattered (upgoing) waves
    KZ_up_2 = (1/cc(ll))^2*OMEGA.^2 - KX.^2;    % KZ squared
    realWaveIndex = (KZ_up_2 >= 0);             % Index of real kz
    KZ_up = sqrt(KZ_up_2.*realWaveIndex);       % Calculate KZ
    Pokxa = Pokxa.*realWaveIndex;               % Mask out evanescent waves. 

    % Preallocate image matrix for individual angles
    iml = zeros(nZ(ll),nX,nA);   

    for ii = 1:nA
        % Progress update
        disp(['Processing angle ' num2str(ii) ' of ' num2str(nA)])

        % Create "local" copy of wavefield for current layer and angle
        Pokx = Pokxa(:,:,ii);                 
    
        % Z-axis wave number for transmitted downward plane wave
        KZ_down = (OMEGA/cc(ll)) * cos(sAngle(ii,ll));

        % Pre-calculate phase shift matrices
        PS_Z = exp(1i*(KZ_down + KZ_up)*dzl(ll));        % Z extrapolation, down and up
        PS_T = exp(1i*(OMEGA(:,1:nX).*txDelay(:,ii).'));  % Compensate for tx delay

        % Step through each z coordinate within layer and focus
        if not(ismember(ll,param.skipLayers))           
            for jj = 1:nZ(ll)            
                % Create image at current depth
                Pox = ifft(Pokx,[],2);       % Inverse transform along x axis
                Pox = Pox(:,1:nX).*PS_T;     % Time shift back according to tx delay
                iml(jj,:,ii) = sum(Pox);     % Sum, corr. to inv. Fourier transform at t=0
                
                % Extrapolate to next depth
                Pokx = Pokx.*PS_Z;
            end 
        end
        
        % Extrapolate original wavefield to next interface
        if ll < nL
            Pokxa(:,:,ii) = Pokxa(:,:,ii).*exp(1i*(KZ_down + KZ_up)*thick(ll));
        end
    end
    
    % Sum across all steering angles and copy to cell structure
    im{ll} = sum(iml,3);
end


%% Assign output
varargout{1} = im;                 % Focused image
varargout{2} = (0:(nX-1))*aPitch;  % X-axis pixel positions
varargout{3} = zz;                 % Z-axis pixel positions (1 cell per layer)

