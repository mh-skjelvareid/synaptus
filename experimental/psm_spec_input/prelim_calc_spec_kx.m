function [Pokx,kx] = calc_spec_kx(Pox,xStep,varargin)
% calc_spec_kx - calculate kx spectrum for 2D dataset
%
%   Usage:
%   [Pokx,kx] = calc_spec_kx(Pox,xStep,...)
%
%   Input parameters:
%   Pox     -    MxN matrix, dimensions (omega,x)
%   xStep   -    spatial sampling step [m]. Assumed to be constant.
%
%   Optional parameter-value input pairs:
%   xFftMult    -    multiplier for t-axis FFT size. Default: 1
%
%   Output parameters:
%   Pokx    -   Spectrum matrix, dimensions (omega,x)
%   kx      -   vector of kx values in output. Contains both positive
%               values (first half) and negative values (second half).
%
%   Notes:
%   The function calculates the Fourier spectrum along the second dimension
%   of Pox, corresponding to the x (spatial) axis. The spectrum is calculated at
%   nextpow2(N) points, where N is the number of x samples (the fft()
%   function is most efficient for input sizes that are powers of 2).
%   Optionally, this number can be modified using the 'xFftMult' parameter.
%   The FFT size is then given by xFftMult * nextpow2(N).
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
p.addParameter('xFftMult',1);           % Multiplier for FFT size in t dir.
p.parse(varargin{:});                   % Parse possible parameter-value pairs
param = p.Results;                      % Transfer results to "param" structure
clear p

%% 
nFFTx = param.xFftMult * 2^nextpow2(nX);            % FFT size x

kxs = (2*pi)/xStep;                                 % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);                                 % Shift to correspond to fft() output

%% Fourier transform along time dimension
Pokx = fft(Pox,nFFTx,2);

