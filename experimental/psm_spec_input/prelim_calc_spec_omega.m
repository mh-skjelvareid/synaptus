function [Pox,omega] = calc_spec_omega(ptx,fs,fLow,fHigh,varargin)
% calc_spec_omega - calculate (bandpass) omega spectrum for 2D dataset
%
%   Usage:
%   [Pox,omega] = calc_spec_omega(ptx,fs,fLow,fHigh,...)
%
%   Input parameters:
%   ptx     -    MxN matrix, dimensions (t,x)
%   fs      -    Sampling frequency [Hz]
%   fLow    -    Lower limit of transducer frequency band [Hz]
%   fHigh   -    Higher limit of transducer frequency band [Hz]
%
%   Optional parameter-value input pairs:
%   tFftMult    -    multiplier for t-axis FFT size. Default: 1
%
%   Output parameters:
%   Pox     -   Spectrum matrix, dimensions (omega,x)
%   omega   -   vector of omega values in output
%
%   Notes:
%   The function calculates the Fourier spectrum along the first dimension
%   of ptx, corresponding to the time axis. The spectrum is calculated at
%   nextpow2(M) points, where M is the number of time samples (the fft()
%   function is most efficient for input sizes that are powers of 2).
%   Optionally, this number can be modified using the 'tFftMult' parameter.
%   The FFT size is then given by tFftMult * nextpow2(M).
%
%   After calculating the "full" FFT, the passband (given by fLow and
%   fHigh) is extracted. Note that only positive omega (one-sided spectrum)
%   are used, i.e. the spectrum is extracted from the "first half" of the
%   output from fft().
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


%% Parse optional input arguments
p = inputParser;                        % Create input parser
p.addParameter('tFftMult',1);           % Multiplier for FFT size in t dir.
p.parse(varargin{:});                   % Parse possible parameter-value pairs
param = p.Results;                      % Transfer results to "param" structure
clear p

%% 
nFFTt = param.tFftMult * 2^nextpow2(size(ptx,1));       % FFT size, time axis

omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);    % Freq. vactor
omega = ifftshift(omega);                                       % Shift as fft
omegaBandIndex = (omega >= 2*pi*fLow) & (omega <= 2*pi*fHigh);

%% Fourier transform along time dimension
Pox = fft(ptx,nFFTt,1);

%% Extract passband
Pox = Pox(omegaBandIndex,:);
omega = omega(omegaBandIndex);
