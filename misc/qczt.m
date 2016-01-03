function g = qczt(x, k, w, a)
%QCZT  Chirp z-transform.
%   G = QCZT(X, M, W, A) is a specialized version of the M-element
%   z-transform of sequence X, where M, W and A are scalars which specify
%   the points along the unity circle on which the z-transform is computed.
%   M is the length of the transform, W is the phase difference between
%   points on the contour, and A is the phase starting point.  More
%   explicitly, the points on the unit circle are described by
%       z = exp(1i*(A + W.^(-(0:M-1)))
%
%   The parameters M, W, and A are optional; their default values are M =
%   length(X), W = -j*2*pi/M, and A = 0.  These defaults cause CZT to
%   return the z-transform of X at equally spaced points around the unit
%   circle, equivalent to FFT(X).
%
%   If X is a matrix, the chirp z-transform operation is applied to each
%   column.
%
%   See also CZT, FFT, FREQZ.

%   References:
%     [1] Oppenheim, A.V. & R.W. Schafer, Discrete-Time Signal
%         Processing,  Prentice-Hall, pp. 623-628, 1989.
%     [2] Rabiner, L.R. and B. Gold, Theory and Application of
%         Digital Signal Processing, Prentice-Hall, Englewood Cliffs, New
%         Jersey, pp. 393-399, 1975.

%% Dimension checking
[m, n] = size(x); oldm = m;
if m == 1, x = x(:); [m, n] = size(x); end

if nargin < 2, k = length(x); end
if nargin < 3, w = -1i .* 2 .* pi ./ k; end
if nargin < 4, a = 0; end

if any([size(k) size(w) size(a)]~=1),
    error(generatemsgid('InvalidDimensions'),'Inputs M, W, A must be scalars.')
end

%% Compute FFT length
nfft = 2^nextpow2(m+k-1);

%% Premultiply data
kk = ( (-m+1):max(k-1,m-1) ).';
kk2 = (kk .^ 2) ./ 2;

ww = exp(1i*w.*(kk2));                  % Chirp filter is 1./ww

nn = (0:(m-1))';
aa = exp(1i*a .* ( -nn ));
aa = aa.*ww(m+nn);
y = x .* aa(:,ones(1,n));

%% Fast convolution via FFT.
fy = fft(  y, nfft );
fv = fft(conj(ww(1:(k-1+m))), nfft );   % Chirp filter
fy = fy .* fv(:,ones(1, n));
g  = ifft( fy );

%% Final multiply
g = g( m:(m+k-1), : ) .* ww( m:(m+k-1),ones(1, n) );

if oldm == 1, g = g.'; end
