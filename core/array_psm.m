function varargout = array_psm(ptxw,fs,txDelay,cc,thick,fLow,fHigh,aPitch,varargin)
% array_psm - Phase Shift Migration for linear arrays with arbitrary shaped transmitted waves
%
%   Usage:
%   [im,...] = array_psm(ptxw,fs,txDelay,cc,thick,fLow,fHigh,aPitch,...)
%
%   Input parameters:
%   ptxw    -   ultrasonic data, dimensions (t,x,txwave). The x dimension
%               corresponds to spatial position of the array elements. The 
%               "txwave" dimension corresponds to different transmit wave 
%               sequences, e.g. plane waves with different steering angles. 
%               If there is only one transmitted wave, ptxw is a 2D matrix, 
%               otherwise, it is 3D. The signal is assumed to be a "raw"
%               ultrasonic signal, e.g. real-valued and not demodulated.
%               It is also assumed that there is no "measurement delay", 
%               i.e. the data starts at t=0, when the first pulse is fired 
%               (zero-pad if necessary).
%   fs      -   Sampling frequency [Hz]
%   txDelay -   2D array of pulse time delays [s], dimensions (x, txwave)
%   cc      -   Vector of compressional wave velocities for each layer [m/s]
%   thick   -   Vector of thickness for each layer [m/s]
%                  (set last boundary equal to end of ROI)
%   fLow    -   Lower limit of transducer frequency band [Hz]
%   fHigh   -   Higher limit of transducer frequency band [Hz]
%   aPitch  -   array pitch (distance between array element centers) [m]
%
%   Optional parameter-value input pairs:
%   tFftMult    -   multiplier for t-axis FFT size. Default: 1
%   xFftMult    -   multiplier for x-axis FFT size. Default: 2
%   zResMult    -   multiplier for final image z resolution. Default: 2
%   wf          -   2-way waveform, for matched filt [vector]. Default: 1 
%   pulseDelay  -   Used if 2-way waveform is not available. Approximate
%                   delay from transducer excitation to pulse maximum [s],
%                   for transmit and receive processes combined.
%                   Default: 0.
%   skipLayers  -   vector w/indices of layers to skip in processing 
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
%   Tromsø, Norway. However, the adaptation to array imaging is not 
%   included in the PhD thesis.
%
%   In order to focus arbitrarily shaped waves transmitted from an array,
%   the transmitted and received (scattered) wavefields are _both_
%   extrapolated downwards from the array position (z=0) into the imaged
%   region (z>0). Note that this requires careful handling of how
%   positive and negative frequencies are extracted after the initial FFT
%   along the time dimension. In this code, downwards travelling waves
%   correspond to negative omega, and upwards travelling waves correspond
%   to positive omega.
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
%   2011-07-07  MHS - Initial version, based on previous PSM code
%
%   2021 Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com / martin.skjelvareid@uit.no


%    A note on code notation:
%    ------------------------
%    In this implementation of PSM, several versions of the discretely
%    sampled wavefield is used. A short explanation of the code notation is
%    included here to make reading of the code easier.
%
%    Code variable          Represents wave field
%    ptxw                   p(t,x,txwave)          [received wavefield]
%    ptxw_d                 p_d(t,x,txwave)        [transmitted wavefield]
%    Poxw_u                 P_u(omega,x,txwave)    [upwards travelling waves]
%    Pokxw_u                P_u(omega,k_x,txwave)  [upwards travelling waves]
%    Pox_u                  P_u(omega,x)           [upwards travelling waves]
%    Poxw_d                 P_d(omega,x,txwave)    [downwards travelling waves]
%    Pokxw_d                P_d(omega,k_x,txwave)  [downwards travelling waves]
%    Pox_d                  P_d(omega,x)           [downwards travelling waves]



%% Parse optional input arguments
p = inputParser;                        % Create input parser
p.addParameter('xFftMult',2);           % Multiplier for FFT size in x dir.
p.addParameter('tFftMult',1);           % Multiplier for FFT size in t dir.
p.addParameter('zResMult',2);           % Multiplier for z resolution
p.addParameter('wf',1)                  % 2-way waveform
p.addParameter('pulseDelay',0)          % Delay from excitation to pulse 
p.addParameter('skipLayers',[])         % Boolean index of layers to skip
p.parse(varargin{:});                   % Parse possible parameter-value pairs
param = p.Results;                      % Transfer results to "param" structure
clear p

%% Calculate dependent variables
[nT,nX,nW] = size(ptxw);        % # time samples, # array elements, # transmitted waves    
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
omegaBandIndex = (omega >= 2*pi*fLow) & (omega <= 2*pi*fHigh);        % Ind. for upwards travelling waves

kxs = (2*pi)/aPitch;                                % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);                                 % Shift corr. to fft. output

nShiftX = floor((nFFTx-nX)/2);                      % Shift for centering output along x-axis

%% Start building transmitted (downgoing) wavefield p_d(t,x,z=0)
% If no waveform is given, use simple time-shifted pulse
if isequal(param.wf,1) && not(param.pulseDelay == 0)
    param.wf = [zeros(round(param.pulseDelay*fs),1);1]; 
end
% Repeat waveform for all x and txwaves
ptxw_d = repmat(param.wf(:),[1,nX,nW]);     

%% Fourier transform along time dimension for both wavefields
Poxw_d = fft(ptxw_d,nFFTt,1);
Poxw_u = fft(ptxw,nFFTt,1);
       
%% Extract frequency bands of interest (correspond to upgoing/downgoing waves)
Poxw_u = Poxw_u(omegaBandIndex,:,:);                    % Pos. omega band
Poxw_d = Poxw_d(circshift(flip(omegaBandIndex),1),:,:); % Neg. omega band

%% Create omega and kx grids
[OMEGA,KX] = ndgrid(omega(omegaBandIndex),kx);    % Grid of pos. omega

%% Flip spectrum for negative omega so that OMEGA matrix matches
% Note: After reversing the omega axis like this, the inverse Fourier
% transform (along the omega axis) corresponds to fft(), not ifft().
Poxw_d = flip(Poxw_d,1);

%% Time shift downgoing wavefield according to tx delays
dt = reshape(txDelay,[1,nX,nW]);
Poxw_d = Poxw_d.*exp(1i*OMEGA(:,1:nX).*dt);

%% Fourier transform in x direction
Pokxw_d = fft(Poxw_d,nFFTx,2);
Pokxw_u = fft(Poxw_u,nFFTx,2);

%% Phase shift migration through each layer
im = cell(nL,1);            % Preallocate image cell structure

for ll = 1:nL
    % Progress update
    disp('*************************')
    disp(['Processing layer ' num2str(ll) ' of ' num2str(ll)])
    disp('*************************')

    % Calculate z-axis wave number KZ for transmitted (downgoing) waves
    KZ_2 = (1/cc(ll))^2*OMEGA.^2 - KX.^2;   % KZ squared
    realWaveIndex = (KZ_2 >= 0);            % Index of real kz
    KZ = sqrt(KZ_2.*realWaveIndex);       % Calculate KZ
    
    % Pre-calculate phase shifts
    PS = exp(1i*KZ*dzl(ll));            % Phase shift for every dz in layer
    PSL = exp(1i*KZ*thick(ll));        % Phase shift for whole layer, up
    
    % Mask out evanescent waves
    Pokxw_u = Pokxw_u.*realWaveIndex;               
    Pokxw_d = Pokxw_d.*realWaveIndex; 

    % Preallocate image matrix for individual transmitted waves
    iml = zeros(nZ(ll),nFFTx,nW);   

    for ii = 1:nW
        % Progress update
        disp(['Processing transmitted wave ' num2str(ii) ' of ' num2str(nW)])

        % Extract local copies of wavefields, for current layer and txwave
        Pokx_d = Pokxw_d(:,:,ii);   
        Pokx_u = Pokxw_u(:,:,ii);

        % Step through each z coordinate in layer, create focused image lines
        if not(ismember(ll,param.skipLayers))           
            for jj = 1:nZ(ll)            
                % Match wavefields
                Pox = ifft(Pokx_u,[],2).*ifft(Pokx_d,[],2);
                
                % Sum, corresponds to inv. Fourier transform at t=0
                iml(jj,:,ii) = sum(Pox);                     
                
                % Extrapolate to next depth
                Pokx_d = Pokx_d.*PS;
                Pokx_u = Pokx_u.*PS;
            end 
        end
        
        % Extrapolate original wavefield to next interface
        if ll < nL
            Pokxw_d(:,:,ii) = Pokxw_d(:,:,ii).*PSL;
            Pokxw_u(:,:,ii) = Pokxw_u(:,:,ii).*PSL;
        end
    end
    
    % Sum across all transmit sequences, x shift, and copy to cell structure
    iml_sum = sum(iml,3);
    im{ll} = circshift(iml_sum,nShiftX,2);
end


%% Assign output
varargout{1} = im;                                  % Focused image
varargout{2} = ((0:(nFFTx-1)) - nShiftX)*aPitch;    % X-axis pixel positions
varargout{3} = zz;                                  % Z-axis pixel positions (1 cell per layer)

