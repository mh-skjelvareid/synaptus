close all
clear all
clc

%% Load data
disp('Loading data'); 
load('../datasets/PlaneScan3D_PlexiAluFBH.mat',...
    'fs','ptxy','tDelay','fLow','fHigh','thick','xStep','yStep','cc')

%% Focus using PSM
disp('Processing data'); 
tic
[im,xIm,yIm,zIm] = psm(ptxy,fs,tDelay,cc,thick,fLow,fHigh,[xStep yStep]);
toc; 

%% Plot raw and focused image
test3D_plot
