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

%% Find elements to calculate for
ri = abs(z) > 1.1*abs(n);

%% Preallocate
sqTerm = zeros(size(z));
amp = zeros(size(z));
phase = zeros(size(z));
H = zeros(size(z));

%% Calculate
sqTerm(ri) = sqrt(z(ri).^2 - n(ri).^2);                                 % Square root term used in several places
amp(ri) = sqrt(2./(pi*sqTerm(ri)));                                  % Amplitude term
phase(ri) = sqTerm(ri) - pi/4 - n(ri)*(pi/2) + n(ri) .* atan(n(ri)./sqTerm(ri));    % Phase term

H(ri) = amp(ri).*exp(-1i*phase(ri));

