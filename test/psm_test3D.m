close all
clear all
clc

%% Load data
disp('Loading data'); fflush(stdout);
load('../datasets/PlaneScan3D_PlexiAluFBH.mat',...
    'fs','ptxy','tDelay','fLow','fHigh','thick','xStep','yStep','cc')

%% Focus using PSM
disp('Processing data'); fflush(stdout);
tic
[im,xIm,yIm,zIm] = psm(ptxy,fs,tDelay,cc,thick,fLow,fHigh,[xStep yStep]);
toc; fflush(stdout);

%% Plot raw and focused image
test3D_plot
