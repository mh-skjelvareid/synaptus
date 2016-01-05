close all
clear all

%% Load data
disp('Loading data');
load('../datasets/PlaneScan3D_PlexiAluFBH.mat',...
    'fs','ptxy','tDelay','fLow','fHigh','fc','thick','xStep','yStep','cc')

interpol = 'chirpz';    % Stolt interp. method, options 'linear', 'chirpz'
% interpol = 'linear';    % Stolt interp. method, options 'linear', 'chirpz'

%% Focus
disp('Processing data');
tic
[im,xIm,yIm,zIm] = mulok(ptxy,fs,tDelay,cc,thick,fLow,fHigh,[xStep yStep],...
    'fc',fc,'interpol',interpol);
toc;

%% Plot raw and focused image
test3D_plot
