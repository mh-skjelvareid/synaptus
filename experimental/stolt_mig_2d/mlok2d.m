function varargout = mlok2d(ptx,fs,mwb,cc,thick,fLow,fHigh,xStep,varargin)
%
%	Usage:
%	[im,...] = mlok2d(ptx,fs,mwb,cc,thick,fLow,fHigh,xStep,...)
% 
%	Input parameters:
%	ptx			-	p(t,x): matrix with ultrasonic data, A-scans along columns
%	fs			-	Sampling frequency [Hz]
%	mwb			-	Measurement Window Begin [s] - time delay between
%					pulse transmission and beginning of data logging
%	cc			-	Vector of compressional wave velocities for each layer [m/s]
%	thick		-	Vector of thickness for each layer [m/s] (set last boundary equal to end of ROI)
%	fLow		-	Lower limit of transducer frequency band [Hz]
%	fHigh		-	Higher limit of transducer frequency band [Hz]
%	xStep		-	x step value (distance between measurements) [m]
%
%	Optional parameter-value input pairs:
%	xFftMult	-	multiplier for FFT size in x dir. Default: 1
%	tFftMult	-	multiplier for FFT size in t dir. Default: 1
%	zFftMult	-	multiplier for z FFT size. Default: 1 
%   upSamp      -   omega upsampling factor (preproc. before interp.) Default: 4
%	hh			-	Transducer impulse response, sampled at fs. Used in matched filter. Default: 1 (no filt.)
%	xStart		-	First x value of scan [m]. Default: 0. Only relevant for "xIm" output
%
%	Output parameters (varargout):
%	im			-	Reconstructed image as cell array - one cell per layer
%	xIm			-	x positions for pixels in im (single vector)
%	zIm			-	z positions for pixels in im (cell array, different z res. for each layer)
%
%   Note: The z resolution of the output image for layer l is set to 
%       dz = (cc(l)/2)/fHigh. 
%   This is done to minimize the required number of z depths for each
%   layer. If this resolution is too coarse, try resampling the image to a
%   higher resolution (see resample, interp)
%
%	2010-06-16	Martin H. Skjelvareid
%
%   2010-10-20 Added a term effectively doubling nFFTz. It is assumed that
%              this will correct aliasing errors seen earlier. 

%	A note on code notation:
%	------------------------
%	What Stolt migration does is essentially to transform the sampled
%	wavefield p(t,z=0,x) to p(t=0,z,x). To understand how this is done,
%	consult f.ex. Margrave (2003) - "Numerical Methods in Exploration
%	Seismology".
%
%	In the process of doing this transformation, many versions of the
%	wavefield are used. A short explanation of the code notation is
%	included here to make reading of the code easier.
%
%	Code variable			Represents wave field
%	ptx						p(t,x)
%	Pox						P(omega,x)
%	Pokx					P(omega,k_x)
%	Pkzkx					P(k_z,k_x)
%	pzx						p(z,x)


%% Parse optional input parameters
p = inputParser;									% Input parsing object
p.addParamValue('xFftMult',1);						% Multiplier for FFT size in x dir.
p.addParamValue('tFftMult',1);						% Multiplier for FFT size in t dir.
p.addParamValue('zFftMult',1)						% Multiplier for FFT size in z dir.
p.addParamValue('upSamp',4)                         % Multiplier for interpolation FFT size
p.addParamValue('hh',1);							% Impulse response
p.addParamValue('xStart',0)							% First x value of scan

p.parse(varargin{:});								% Parse all parameter-value pairs
param = p.Results;									% Store parsed results in structure "param"
clear p

%% Get sizes of input
[nT,nX] = size(ptx);								% Number of time samples / x positions
nL = length(cc);									% Number of layers

%% Derived variables
dzl = (cc/2)./(fHigh);                              % Z resolution for each layer
zOffset = mwb*(cc(1)/2);                            % Z offset due to measurement window begin
zIF = cumsum([0; thick(:)]);						% z position of each interface

%% Set FFT sizes
nFFTx = param.xFftMult * 2^nextpow2(nX);				% FFT size in x dimension
nFFTt = param.tFftMult * 2^nextpow2(nT);                % FFT size in time dimension 

%% Create omega and k vectors
dOmega = ((2*pi*fs)/nFFTt);                             % Omega step size
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*dOmega;       % Angular frequency vector
omega = ifftshift(omega);                               % Shift to correspond to fft output

kxs = (2*pi)/xStep;										% Sampling wavenumber, x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);		% Vector of wave numbers in x dir
kx = ifftshift(kx);                                     % Shift to correspond to fft output

%% Calculate bandpass mask for omega 
omegaBandIndex = (omega >= -(2*pi*fHigh)) & (omega <= -(2*pi*fLow));	% Index vector for transducer freq. band (only neg. freq.)
omegaBand = omega(omegaBandIndex);
nOmega = nnz(omegaBandIndex);

%% Fourier transform along time dimension
Pox = fft(ptx,nFFTt,1);						% P(omega,x)
clear ptx									% Free up space

%% Cut out frequency band
Pox = Pox(omegaBandIndex,:);

%% Optional matched filtering
if param.hh ~= 1
    HH = fft(param.hh,nFFTt,1);                         % Fourier transformed impulse response
    Pox = Pox.*repmat(conj(HH(omegaBandIndex)),1,nX);   % Mathced filter (multiply with conjugate)
end

%% Interpolate band to higher resolution in frequency (useful for Stolt interpol.)
if param.upSamp ~= 1
    Pox = fft(ifft(Pox),param.upSamp*nOmega);
end
omegaIP = (0:(nOmega*param.upSamp-1))*(dOmega/param.upSamp);
omegaIP = -omegaIP(end:-1:1) + omegaBand(end);

%% MWB compensation (MUST be performed AFTER upsamping)
Pox = Pox.*repmat(exp(-1i*omegaIP(:)*mwb),[1 nX]);				% Shift wavefield "upwards" in time to compensate for measurement delay

%% Fourier transform in x direction
Pokx = fft(Pox,nFFTx,2);
clear Pox														

%% Create PSM grids (correspond to high resolution P(omega,kx))
[OMEGA_psm,KX_psm] = ndgrid(omegaIP,kx);

%% Extraplolate wavefield to each interface, and image using Stolt transf.
im = cell(nL,1);
zIm = cell(nL,1);

for ii = 1:nL
    nFFTz = param.zFftMult * 2^nextpow2(nT*(fHigh/(fs/2)));     % Size of kz vector ised in interpolation
    kzs = (2*pi)/dzl(ii);									% The maximum kz in the final image is given by the highest frequency and the lowest sound speed
    kz = (0:(nFFTz-1))*(kzs/nFFTz);							% Create a vector of the kz values of interest
    kz = -kz(end:-1:1);                                     % Reverse to correspond to neg. omega
    
    % "Stolt grids" (correspond to low resolution P(kz,kx) wavefield)
    [KZ_stolt,KX_stolt] = ndgrid(kz,kx);
    KK_stolt = (KZ_stolt.^2 + KX_stolt.^2);
    KK_stolt = -sqrt(KK_stolt .* (KK_stolt > 0));
    
    Akzkx = 1./(1 + (KX_stolt.^2)./(KZ_stolt.^2));          % Amplitude factor 
    Akzkx(isnan(Akzkx)) = 0;                                % Set NaN values to zero
    clear KX_stolt 
    
	% Calc. omega values to be interpolated for
	OMEGA_stolt = (cc(ii)/2)*KK_stolt;	    
    
	% Interpolate for each kx
	Pkzkx = complex(zeros(nFFTz,nFFTx));							% Preallocate for interpolated wavefield
    for jj = 1:nFFTx
% 		Pkzkx(:,jj) = interp1q(omegaIP(:),Pokx(:,jj),OMEGA_stolt(:,jj));
        Pkzkx(:,jj) = nakeinterp1(omegaIP(:),real(Pokx(:,jj)),OMEGA_stolt(:,jj));
        Pkzkx(:,jj) = Pkzkx(:,jj) + 1i*nakeinterp1(omegaIP(:),imag(Pokx(:,jj)),OMEGA_stolt(:,jj));
%         Pkzkx(:,jj) = nakeinterp1(omegaIP(:),Pokx(:,jj),OMEGA_stolt(:,jj));
    end
	
	% Set nonvalid elements to zero
%     Pkzkx(isnan(Pkzkx)) = 0;
    nonValid = (isnan(Pkzkx)) | (OMEGA_stolt < omegaIP(1)) | (OMEGA_stolt > omegaIP(end));
    Pkzkx(nonValid) = 0;
    

	% Amplitude scaling (because of variable change from omega to kz)
	Pkzkx = Pkzkx.*(cc(ii)*Akzkx);		

	% If first layer, phase shift back according to MWB (first z line corresponds to beginning of measurement window)
	if ii == 1
		Pkzkx = Pkzkx.*exp(1i*KZ_stolt*zOffset);
	end
				
	% Inverse transform from (kz,kx) to (z,x)
	pzx = ifft2(Pkzkx);		
	    
	% Calculate number of image lines in current layer
    if ii == 1
		nLinesZ = min([ceil((thick(ii)-zOffset)/dzl(ii)) nFFTz]) + 1;       % Calculate number of lines. Add one extra for slight overlap
        zIm{ii} = (0:(nLinesZ-1))*dzl(ii) + zOffset;                        % Z axis vector
    else
		nLinesZ = ceil(thick(ii)/dzl(ii)) + 1;                              % Calculate number of lines. Add one extra for slight overlap
        zIm{ii} = (0:(nLinesZ-1))*dzl(ii) + zIF(ii);                        % Z axis vector
    end    
    
	% Cut out part of pzx corresponding to current layer
    im{ii} = pzx(1:nLinesZ,1:nX);

    % Migrate to next layer (if not last layer)
	if ii < nL
		KK_psm = OMEGA_psm/(cc(ii)/2); 
		KZ2_psm = (KK_psm.^2 - KX_psm.^2);
		KZ_psm = sqrt(KZ2_psm .* (KZ2_psm > 0));        
		Pokx = Pokx .* exp(-1i*KZ_psm*thick(ii));
	end
end

%% Output
xIm = (0:(nX-1))*xStep + param.xStart;		% X axis vector
varargout{1} = im;		
varargout{2} = xIm;
varargout{3} = zIm;
