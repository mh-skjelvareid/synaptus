close all
clearvars
clc

%% Show test suite header
disp('###############################################################################')
disp('#####                                                                     #####')
disp('#####                     SYNAPTUS TOOLBOX TEST SUITE                     #####')
disp('#####                                                                     #####')
disp('#####         No errors indicates that everything is working OK.          #####')
disp('#####          To view data plots: Run each test individually.            #####')
disp('#####                                                                     #####')
disp('###############################################################################')
disp(' ')

%% Run PSM tests
psm_test2D
psm_test3D

%% Run MULOK tests
mulok_test2D
mulok_test3D

%% Run CPSM tests
cpsm_test2D
cpsm_test3D
cpsm_test3D_ObjectsInPipe

%% Run "full matrix" PSM tests
array_psm_testFullMatrix_1layer
array_psm_testFullMatrix_2layer

