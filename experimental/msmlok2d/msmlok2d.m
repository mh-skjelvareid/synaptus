function varargout = msmlok2d(ptx,fs,mwb,cc,thick,fLow,fHigh,xsStep,xrStep,xsrOffset,varargin)
% MSMLOK2D - omega-k migration for 2D multistatic data (full-matrix array)
% 
%	Usage:
%	[im,...] = msmlok2d(ptx,fs,mwb,c,fLow,fHigh,xStep,...)
%
%	Input parameters:
%   ptx			-	p(t,xr,xs): 3D matrix with ultrasonic data
%	fs			-	Sampling frequency [Hz]
%	mwb			-	Measurement Window Begin [s] - time delay between
%					pulse transmission and beginning of data logging
%	cc			-	Wave velocity [m/s]
%	thick		-	Vector of thickness for each layer [m/s] (set last boundary equal to end of ROI)
%	fLow		-	Lower limit of transducer frequency band [Hz]
%	fHigh		-	Higher limit of transducer frequency band [Hz]
%	xsStep		-	x step between source positions [m]
%	xrStep		-	x step between receiver positions [m]
%	xsrOffset	-	offset between first source position and first receiver position
%
%	Optional parameter-value input pairs:
%	xStart		-	Lowest receiver position [m]. Default: 0. Only relevant for "xIm" output
%	ipMethod	-	interpolation method [string] (see INTERP1 ) Default: 'linear'
%	hh			-	Transducer impulse response, sampled at fs. Used in matched filter. Default: 1 (no filt.)
%	tFftMult	-	multiplier for FFT size in t dir. Default: 1
%	xFftMult	-	multiplier for FFT size in both x dir. Default: 2
%	zFftMult	-	multiplier for z FFT size. Default: 1 
%   ipFftMult   -   multiplier for interp. to higher resolution before interpolation from omega to kz
%	zRes		-	Z resolution in the final image. The default resolution
%					is given by (min(cc)/2)/fHigh, but if specified, the
%					final image is interpolated to a resolution of zRes.
%
%	Output parameters (varargout):
%	im			-	Reconstructed image
%	xIm			-	x positions for pixels in im
%	zIm			-	z positions for pixels in im
%
%	2010-01-21	Martin H. Skjelvareid

%	A note on code notation:
%	------------------------
%	In the process of doing this transformation, many versions of the
%	wavefield are used. A short explanation of the code notation is
%	included here to make reading of the code easier.
%
%	Code variable			Represents wave field
%	ptx						p(t,xr,xs)
%	Pox						P(omega,xr,xs)
%	Pokx					P(omega,kxr,kxs)
%	Pkzkx					P(kz,kx,kxs)
%	Pkzx					P(kz,x)
%	pzx						p(z,x) (Image)


%% Parse optional input parameters
p = inputParser;
p.addParamValue('xStart',0)							% Lowest receiver x value
p.addParamValue('ipMethod','linear');				% Interoplation method
p.addParamValue('xFftMult',2);						% Multiplier for FFT size in x dir.
p.addParamValue('tFftMult',1);						% Multiplier for FFT size in x dir.
p.addParamValue('ipFftMult',8)						% Multiplier for interpolation FFT size
p.addParamValue('zFftMult',1)						% Multiplier for FFT size in z dir.
p.addParamValue('hh',1);							% Impulse response
p.addParamValue('zRes',[]);							% Z resolution

p.parse(varargin{:});
param = p.Results;
clear p

%% Get size of input
[nT,nXr,nXs] = size(ptx);
nL = length(cc);

%% Set FFT sizes
nFFTxr = param.xFftMult * 2^nextpow2(nXr);								% FFT size in x dimension
nFFTxs = param.xFftMult * 2^nextpow2(nXs);								% FFT size in x dimension

nFFTt = param.tFftMult * 2^nextpow2(nT+length(param.hh)-1);				% FFT size in time/z dimension
nFFToip = param.ipFftMult * 2^nextpow2(nFFTt*((fHigh-fLow)/fs));		% FFT size for interpolation in omega (prop. to rel. bandwidth)
nFFTz = param.zFftMult * 2^nextpow2((nT+length(param.hh)-1)*(fHigh/fs)*(max(cc)/cc(1)));	% Size of kz vector ised in interpolation

%% Calc. z offset corresponding to measurement window begin
zOffset = mwb*(cc(1)/2);

%% Create omega and k vectors
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);			% Angular frequency vector

kxr_s = (2*pi)/xrStep;													% Sampling wavenumber for xr
kxr = ((0:(nFFTxr-1)) - floor(nFFTxr/2))*(kxr_s/nFFTxr);				% Vector of xr wave numbers

kxs_s = (2*pi)/xsStep;													% Sampling wavenumber for xr
kxs = ((0:(nFFTxs-1)) - floor(nFFTxs/2))*(kxs_s/nFFTxs);				% Vector of wave numbers in x dir

kzs = (2*pi*fHigh)/(min(cc)/2);									% The maximum kz in the final image is given by the highest frequency and the lowest sound speed
kzBand = (-(nFFTz):-1)*(kzs/nFFTz);								% Create a vector of the kz values of interest
zStep = (2*pi)/kzs;												% Step size is given by "sampling wavenumber" kzs. Do not confuse with param.zRes.

%% Calculate bandpass masks for omega and kz (to only include transducer frequency band)
omegaBandIndex = (omega >= -(2*pi*fHigh)) & (omega <= -(2*pi*fLow));	% Index vector for transducer freq. band (only neg. freq.)

%% Fourier transform along time dimension
Pox = fftshift(fft(ptx,nFFTt,1),1);								% P(omega,xr,xs)
clear ptx														% Free up space

%% Matched filtering (optional)
if param.hh ~= 1
	H = fftshift(fft(param.hh,nFFTt,1),1);						% Fourier transformed impulse response
	Pox = Pox.*repmat(conj(H),[1 nXr,nXt]);						% Multiply with data column by column
end

%% Cut out frequency band of interest
Pox = Pox(omegaBandIndex,:,:);

%% Fourier transform in x directions
Pokx = fftshift(fft(Pox,nFFTxr,2),2);
Pokx = fftshift(fft(Pokx,nFFTxs,3),3);
clear Pox														% Free up space

%% Interpolate band to higher resolution in frequency (useful for Stolt interpol.)
Pokx = fft(ifft(Pokx),nFFToip);
omegaBand = omega(omegaBandIndex);
omegaBandWidth = abs(omegaBand(end) - omegaBand(1));
omegaIP = (0:(nFFToip-1))*(omegaBandWidth/nFFToip) + omegaBand(1);

%% Calculate omega and k matrices
[KZ,KXR,KXS] = ndgrid(kzBand,kxr,kxs);										% Wave number matrices

OMEGAIP2 = (KZ.^2 + 2*(KXS.^2 + KXR.^2) - ((KXS.^2 - KXR.^2).^2)./(KZ.^2));		% Note: Sound speed is not included here, but is rather multiplied for each layer
realWaveIndex = OMEGAIP2 >= 0;											% Identify positive elements
OMEGAIP = -sqrt(OMEGAIP2.*realWaveIndex);

clear OMEGAIP2 KXS KXR													% Free up space

% "PSM grids" (correspond to high resolution P(omega,kx))
[OMEGA_psm,KXR_psm,KXS_psm] = ndgrid(omegaIP,kxr,kxs);

%% Phase shift to compensate for measurement time delay and source offset
Pokx = Pokx.*exp(-1i*(OMEGA_psm*mwb + KXS_psm*xsrOffset));

%% Extrapolate wavefield to each interface, and image
im = complex(zeros(nFFTz,nFFTxr));                                         % Preallocate for final image
Pkzkx = complex(zeros(nFFTz,nFFTxr,nFFTxs));							% Preallocate for interpolated wavefield
layerBeginIndex = 0;

for ii = 1:nL
	% Calculate number of image lines in each layer
	if ii == 1
		nLinesZ = max([0 round((thick(ii)-zOffset)/zStep)]);
	else
		nLinesZ = round(thick(ii)/zStep);
	end

	% Calc. omega values to be interpolated for - multiply with local sound speed
	OMEGAIP_loc = (cc(ii)/2)*OMEGAIP;
	
	% Interpolate for each kxr and kxs
    for jj = 1:nFFTxs
        for ll = 1:nFFTxr
			kxrIndex = mod((ll+jj-2),nFFTxr)+1;						% Shift due to substitution of kxr with kx
			Pkzkx(:,kxrIndex,jj) = interp1q(omegaIP(:),Pokx(:,ll,jj),OMEGAIP_loc(:,ll,jj));
        end
    end
    
    % Remove NaN elements due to interpolation out of range
    Pkzkx(isnan(Pkzkx)) = 0;
            
	% Scaling
	% Not yet implemented - not essential to focusing
	
	% Sum across xs dimension
	Pkzkx = sum(Pkzkx,3);
    
	% If first layer, phase shift back according to MWB (first z line corresponds to beginning of measurement window)	
	if ii == 1
		Pkzkx = Pkzkx.*exp(1i*KZ(:,:,1)*(mwb*(cc(1)/2)));
	end
	
	% Inverse transform
    pzx = ifft2(ifftshift(Pkzkx));
    
	% Cut out part of pzx corresponding to current layer
	im((1:nLinesZ) + layerBeginIndex,:) = pzx(1:nLinesZ,:);
	layerBeginIndex = layerBeginIndex + nLinesZ;	% Update index to beginning of next layer
	
	if ii < nL
		% Migrate to next layer
		KK_psm = OMEGA_psm/(cc(ii));
		KZR2_psm = (KK_psm.^2 - KXR_psm.^2);
		KZS2_psm = (KK_psm.^2 - KXS_psm.^2);
		rwi_psm = (KZR2_psm > 0) & (KZS2_psm > 0);
		KZ_psm = sqrt(KZR2_psm .* rwi_psm) + sqrt(KZS2_psm .* rwi_psm);
		Pokx = Pokx .* exp(-1i*KZ_psm*thick(ii)) .* rwi_psm;
	end	
end

nZ = layerBeginIndex;		% Total number of image planes 
im = im(1:nZ,:);			% Crop image

%% Optional interpolation to another z resolution
if not(isempty(param.zRes))
	ipFactor = zStep/param.zRes;                                                    % Calc. interpolation factor
	zStep = zStep/ipFactor;                                                         % Adjust zStep
	IM = fftshift(fft(im,[],1),1);													% Fourier transform image in z dir.
	if ipFactor < 1                                                                 % If coarser resolution,
		im = ifft(ifftshift(IM(1:round(nZ*ipFactor),:,:),1),[],1);                  % Truncate and inverse transform 
    else                                                                            % If finer resolution
		im = ifft(ifftshift([IM; zeros(round(nZ*(ipFactor-1)),nFFTxr)],1),[],1);	% Zeropad and inverse transform
	end
end

%% Output
xShift = round((nFFTxr-nXr)/2);								% x shift to center image
xIm = ((0:(nFFTxr-1))-xShift)*xrStep + param.xStart;		% X axis vector
zIm = (0:(size(im,1)-1))*zStep + zOffset;                   % Z axis vector

varargout{1} = circshift(im,[0 xShift]);
varargout{2} = xIm(:);
varargout{3} = zIm(:);
