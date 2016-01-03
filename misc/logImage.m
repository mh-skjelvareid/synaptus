function imOut = logImage(imIn)
% logImage - return normalized logarithimc image (dB)
%
%	Martin H. Skjelvareid, 2010-04-27

imOut = 20*log10(abs(imIn));
imOut = imOut - max(max(imOut));