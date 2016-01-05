function varargout = cpsm(ptpz,fs,tDelay,cc,fLow,fHigh,sss,r0,...
    rStep,rStart,rEnd,varargin)
% CPSM - Phase shift migration for circular geometry
%
%   Usage:
%   [im,phiIm,rIm,zIm] = cpsm3d(ptpz,fs,tDelay,cc,fLow,fHigh,sss,r0,...
%                      rStep,rStart,rEnd)
%
%   Input arguments:
%   ptpz    -   2D or 3D ultrasonic data, dimensions (t,phi) or (t,phi,z)
%   fs      -   sampling frequency [Hz]
%   tDelay  -   measurement window begin [s]
%   cc      -   wave velocity [m/s]
%   fLow    -   lower edge frequency of transducer passband [Hz]
%   fHigh   -   upper edge frequency of transducer passband [Hz]
%   sss     -   spatial step size. 2D: [phiStep], 3D: [phiStep zStep]
%   r0      -   scanning radius [m]
%   rStep   -   imaging range step size [m]
%   rStart  -   start of imaging range [m]
%   rEnd    -   end of imaging range [m]
%
%   Optional parameter-value pairs
%   'transFunc' - Choice of transfer func. for wavefield extrapolation. Options:
%                 'haun'   - Approx. transfer function [fastest, default]
%                 'gardner'- Approx. Hankel trans. func. [slower, more accurate]
%                 'exact'  - Exact Hankel transfer function [slowest]
%
%   Output arguments:
%   im      -   complex-valued focused image
%               Dimensions: [r,phi] (2D), [r,phi,z]
%   phiIm   -   phi coordinates for im
%   rIm     -   range coordinates for im
%   zIm     -   z coordinates for im
%
%   Notes:
%   - The 'haun' transfer function is described in "Efficient three-dimensional
%   imaging from a small cylindrical aperture," by Haun, M.A. et al., in
%   IEEE Trans. Ultrason., Ferroelect., Freq. Control, 49(7), 2002
%   - The 'gardner' approximate Hankel function is described in "An Accurate
%   Closed-Form Approximate Representation for the Hankel Function of the Second
%   Kind", by J. Gardner and R. E. Collin, IEEE Trans. Antennas Propag., 48(10),
%   2000
%
%   The derivation of the CPSM algorithm is described in Chapter 6 in the PhD
%   thesis "Synthetic aperture ultrasound imaging with application to interior
%   pipe inspection" (2012) by Martin H. Skjelvareid, University of Tromsø,
%   Norway.
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
%   2011-06-30  M. H. Skjelvareid - Initial version
%   2015-12-14  M. H. Skjelvareid - Adapt for both 2D and 3D data
%
%   Copyright (C) 2016  Martin H. Skjelvareid
%   martin.hansen.skjelvareid@gmail.com

%% Parse optional input arguments
p = inputParser;                                    % Create input parser
validator = @(str) any(strcmp(str,{'haun','gardner','exact'}));
p.addParamValue('transFunc','haun',validator);      % Migration transfer func.
p.addParamValue('tFftMult',1);                      % Multip. FFT in t dir.
p.parse(varargin{:});                               % Parse param.-value pairs
param = p.Results;                                  % Transfer res. to structure
clear p

%% Dimension testing (2D is a special case)
dataIs3D = (ndims(ptpz) == 3);

% Check that step sizes are correctly given
if dataIs3D && length(sss) == 1
    error('The sss input parameter must have 2 elements for 3D data sets')
end

%% Auxillary variables
% Input size
if dataIs3D
    [nT,nPhi,nZ] = size(ptpz);
else
    [nT,nPhi] = size(ptpz);
    nZ = 1;
end

% Focus ranges
rr = rStart:rStep:rEnd;
nR = length(rr);

nFFTt = param.tFftMult*2^nextpow2(nT);  % FFT size time
nFFTn = 2^nextpow2(nPhi);               % FFT size phi
nFFTz = 2^nextpow2(nZ);                 % FFT size z

%% Fourier transform
if dataIs3D
    Ponkz = fftn(ptpz,[nFFTt,nFFTn,nFFTz]);
else
    Ponkz = fftn(ptpz,[nFFTt,nFFTn]);
end

% Omega vector
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);
omega = ifftshift(omega);           % Shift to correspond to FFT output

% Angular wave number vector
nns = (2*pi)/(sss(1));             % Calculate phi "sampling wave number"
nn = ((0:(nFFTn-1)) - floor(nFFTn/2))*(nns/nFFTn);
nn = ifftshift(nn);                % Shift to correspond to FFT output

% Z wave number vector
if dataIs3D
    kz_s = (2*pi)/sss(2);                                 % Z samp. wave number
    kz = ((0:(nFFTz-1)) - floor(nFFTz/2))*(kz_s/nFFTz);   % Z-dir wave numbers
    kz = ifftshift(kz);                                   % Shift to FFT output
end

% Crop spectrum to transducer frequency band
omegaBandIndex = (omega >= -(2*pi*fHigh)) & (omega <= -(2*pi*fLow));
Ponkz = Ponkz(omegaBandIndex,:,:);

%% Create spectral coordinate grids
if dataIs3D
    [OMEGA,NN,KZ] = ndgrid(omega(omegaBandIndex),nn,kz);% 3D coordinate grids
    KR2 = (OMEGA/(cc/2)).^2 - KZ.^2;                    % k_r squared
    realWaveIndex = (KR2 >= 0);                         % Index of real kr
    KR = sqrt(KR2.*realWaveIndex);                      % Calculate k_r
else
    [OMEGA,NN] = ndgrid(omega(omegaBandIndex),nn);      % 2D coordinate grids
    KR = abs(OMEGA)/(cc/2);                             % Calculate k_r
end

%% Time shift back to t=0
Ponkz = Ponkz.*exp(-1i*OMEGA*tDelay);

%% Preparations before wavefield extrapolation
if strcmp(param.transFunc,'gardner') || strcmp(param.transFunc,'exact')
    if strcmp(param.transFunc,'gardner')    % Use approximate Hankel function
        Hr0 = approxHankel2(KR*r0,NN);      % See subfunction
    else                                    % Exact
        Hr0 = besselh(NN,2,KR*r0);          % Hankel function of the 2. kind
    end
    Ponkz = Ponkz./Hr0;     % Divide wave field by Hankel func. for scan radius
    Ponkz(Hr0 == 0) = 0;    % Set non-valid elements to zero
end

%% Extrapolate wave field to each image range
im = complex(zeros(nR,nPhi,nZ));        % Preallocate matrix for foc. image

for ii = 1:nR
    switch param.transFunc
        case 'haun'
            KR_mod = KR.^2 - (NN.^2)/(r0*rr(ii));   % "Modified" k_r, squared
            KR_mod = sqrt(KR_mod.*(KR_mod>0));      % Sq. root, check real-val.
            GG = sqrt(rr(ii)/r0).*exp(-1i*(rr(ii)-r0)*KR_mod);
            PonkzShifted = Ponkz.*GG;               % Shift wavefield to r

        case 'gardner'
            Hr = approxHankel2(KR*rr(ii),NN);       % Approx. Hankel func.
            PonkzShifted = Ponkz.*Hr;               % Shift wavefield to r

        case 'exact'
            Hr = besselh(NN,2,KR*rr(ii));           % Exact Hankel func.
            PonkzShifted = Ponkz.*Hr;               % Shift wavefield to r

    end % Switch

    Ponkz_t0 = sum(PonkzShifted,1);                 % IFFT at t=0

    if dataIs3D
        imPlane = ifft(ifft(Ponkz_t0,[],2),[],3);   % IFFT along phi and z
        im(ii,:,:) = imPlane(1,1:nPhi,1:nZ);        % Crop to original dim.
    else
        imLine = ifft(Ponkz_t0,[],2);               % IFFT along phi
        im(ii,:) = imLine(1,1:nPhi);                % Crop to original dim.
    end

    if not(mod(ii,ceil(nR/20)))                     % Update progress
        disp([num2str(ii) ' of ' num2str(nR) ' image ranges processed.'])
    end
end

disp('All image ranges processed')

%% Output
varargout{1} = im;
varargout{2} = ((0:(nPhi-1)) - (nPhi-1)/2)*sss(1);
varargout{3} = rr;
if dataIs3D
    varargout{4} = (0:(nZ-1))*sss(2);
end

end % Function cpsm

%% SUBFUNCTION approxHankel2

function H = approxHankel2(z,n)
% approxHankel2 - accurate approximation to Hankel function of second kind
%
%   Usage: H = approxHankel2(z,n)
%
%   approxHankel2 returns an accurate approximation of the Hankel function
%   of the second kind of order n, H2_n(z). z and n must be of the same
%   size. Where abs(n) > abs(z), a zero is returned.
%
%   Reference: Judd Gardner and R. E. Collin, "An Accurate Closed-Form
%   Approximate Representation for the Hankel Function of the Second Kind",
%   IEEE Transactions on Antennas and Propagation, vol. 48, no. 10, october
%   2000
%
%   2011-06-30  Martin H. Skjelvareid

% Find elements to calculate for
ri = abs(z) > 1.1*abs(n);

% Preallocate
sqTerm = zeros(size(z));
amp = zeros(size(z));
phase = zeros(size(z));
H = zeros(size(z));

% Calculate
sqTerm(ri) = sqrt(z(ri).^2 - n(ri).^2);     % Sq. root term used several places
amp(ri) = sqrt(2./(pi*sqTerm(ri)));                     % Amplitude term
phase(ri) = sqTerm(ri) - pi/4 - n(ri)*(pi/2) + ...      % Phase term
    n(ri) .* atan(n(ri)./sqTerm(ri));

H(ri) = amp(ri).*exp(-1i*phase(ri));        % Calculate approx. Hankel func.

end
