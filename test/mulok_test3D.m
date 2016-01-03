close all
clear all

%% Load data
disp('Loading data'); fflush(stdout);
load('../datasets/PlaneScan3D_PlexiAluFBH.mat',...
    'fs','ptxy','tDelay','fLow','fHigh','fc','thick','xStep','yStep','cc')

interpol = 'chirpz';    % Stolt interpolation method, options 'linear', 'chirpz'

%% Focus
disp('Processing data'); fflush(stdout);
tic
[im,xIm,yIm,zIm] = mulok(ptxy,fs,tDelay,cc,thick,fLow,fHigh,[xStep yStep],...
    'fc',fc,'interpol',interpol);
toc; fflush(stdout);

%% Plot raw and focused image
test3D_plot
