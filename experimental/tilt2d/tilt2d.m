function varargout = tilt2d(ptx,fs,mwb,cc,fLow,fHigh,xStep,aa,bb,varargin)
% TILT2D - transform wavefield from line z=0 to line z = ax+b
%
%	Usage:
%	[ptxTilt,...] = tilt2d(ptx,fs,mwb,cc,fLow,fHigh,xStep,aa,bb,...)
%
%	Input parameters:
%	ptx			-	p(t,x): matrix with ultrasonic data, A-scans along columns
%	fs			-	Sampling frequency [Hz]
%	mwb			-	Measurement Window Begin [s] - time delay between
%					pulse transmission and beginning of data logging
%	cc			-	Compressional wave velocity [m/s]
%	thick		-	Vector of thickness for each layer [m/s] (set last boundary equal to end of ROI)
%	fLow		-	Lower limit of transducer frequency band [Hz]
%	fHigh		-	Higher limit of transducer frequency band [Hz]
%	xStep		-	x step value (distance between measurements) [m]
%	aa			-	line slope parameter of surface
%	bb			-	line offset of surface [m] 
%	
%	Optional parameter-value input pairs:
%	ipMethod	-	interpolation method [string] (see INTERP1Q) Default: 'linear'
%	hh			-	Transducer impulse response, sampled at fs. Used in matched filter. Default: 1 (no filt.)
%	tFftMult	-	multiplier for FFT size in t dir. Default: 1
%	xFftMult	-	multiplier for FFT size in x dir. Default: 1
%	ipFftMult	-	multiplier for z FFT size. Default: 4 
%
%	Output parameters (varargout):
%	ptxTilt		-	Reconstructed image [complex]
%	xStepNew	-	new x positions for pixels in im
%
%	This function is used to extrapolate the wavefield measured at z=0 to
%	the line z = aa*x + bb. This is useful if the scan line is not parallel
%	to the object under test. The output wave field is what would have been
%	measured has the transducer been scanned along the surface.
%
%	Note that the returned wavefield is complex - this is because it is
%	nonzero only for negative omega, following the assumption of one-way
%	waves. The wavefield can still be used directly in other migration
%	algorithms
%
%	2010-29-04 Martin H. Skjelvareid

%% Parse optional input parameters
p = inputParser;									% Input parsing object
p.addParamValue('ipMethod','linear');				% Interoplation method
p.addParamValue('xFftMult',1);						% Multiplier for FFT size in x dir.
p.addParamValue('tFftMult',1);						% Multiplier for FFT size in t dir.
p.addParamValue('ipFftMult',4);						% Multiplier for FFT size in z dir.
p.addParamValue('hh',1);							% Impulse response

p.parse(varargin{:});								% Parse all parameter-value pairs
param = p.Results;									% Store parsed results in structure "param"
clear p

%% Get sizes of input
[nT,nX] = size(ptx);								% Number of time samples / x positions

%% Set FFT sizes
nFFTx = param.xFftMult * 2^nextpow2(nX);								% FFT size in x dimension
nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.hh)-1);				% FFT size in time dimension 
nFFTxip = param.ipFftMult*nFFTx;								

%% Create omega and k vectors
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);	% Angular frequency vector

kxs = (2*pi)/xStep;												% Sampling wavenumber, x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);				% Vector of wave numbers in x dir
kxip = ((0:(nFFTxip-1)) - floor(nFFTxip/2))*(kxs/nFFTxip);		% Vector of wave numbers in x dir, interpolated (high res.)

%% Calculate bandpass mask for omega 
omegaBandIndex = (omega >= -(2*pi*fHigh)) & (omega <= -(2*pi*fLow));	% Index vector for transducer freq. band (only neg. freq.)
nOmega = nnz(omegaBandIndex);

%% Fourier transform along time dimension
Pox = fftshift(fft(ptx,nFFTt,1),1);								% P(omega,x)
HH = fftshift(fft(param.hh,nFFTt,1),1);							% Fourier transformed impulse response
clear ptx														% Free up space

%% Cut out frequency band
Pox = Pox(omegaBandIndex,:);

%% Matched filtering
Pox = Pox.*repmat(conj(HH(omegaBandIndex)),1,nX);								

%% MWB compensation 
Pox = Pox.*repmat(exp(-1i*omega(omegaBandIndex)*mwb),[1 nX]);				% Shift wavefield "upwards" in time to compensate for measurement delay

%% Fourier transform in x direction (high resolution for interpolation)
Pokx = fftshift(fft(Pox,nFFTxip,2),2);
clear Pox														

%% Create grids
kk = omega(omegaBandIndex)/cc;					% Wave number vector
[KK,KX] = ndgrid(kk,kx);						% Grid of kk and kk for interpolated values
KZ = (4*KK.^2 - KX.^2);							% Grid for squared kz
KZ = sqrt(KZ .* (KZ > 0));						% kz, mask out values for evanescent waves

%% Calculate new KX and scale factor
TMP = (-KX.^2 + 4*(1+aa^2)*KK.^2);				% Temporary calc. variable
TMP = sqrt(TMP.*(TMP>0));						% Square root, mask out values for evanescent waves
KXTILT = (1/(1+aa^2))*(KX + aa * TMP);			% Kx for tilted wavefield
Akzkx = (1/(1+aa^2))*(1 - (aa*KX)./TMP);		% Scale factor for tilted wavefield

%% Rotate by interpolation in omega-k domain
PokxTilt = complex(zeros(nOmega,nFFTx));						
for ii = 1:nOmega
	PokxTilt(ii,:) = interp1(kxip,Pokx(ii,:),KXTILT(ii,:),'linear',0);
end

%% Scaling due to transformation
PokxTilt = PokxTilt.*Akzkx;

%% Extrapolate wavefield down to interface 
PokxTilt = PokxTilt.*exp(-1i*KZ*bb);

%% Insert tilted wavefield into full frequency range
PokxTiltFull = zeros(nFFTt,nFFTx);
PokxTiltFull(omegaBandIndex,:) = PokxTilt;
ptxTilt = ifftn(ifftshift(PokxTiltFull));

%% Calculate new xStep value
xStepNew = xStep/cos(atan(aa));

%% Output
varargout{1} = ptxTilt(1:nT,1:nX);		% Tilted image, cropped to original dimensions
varargout{2} = xStepNew;				% New xStep value
