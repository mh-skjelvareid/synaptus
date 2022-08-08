%   psm_2d_batch - batch migration using PSM and preprocessed data
%
%   Note that this script uses two local functions,
%   loadAndConcatenateBlocks() and calc_spec_kx(), defined at the end of
%   the file.
%
%   2021-09-16  MHS - Initial version
%

%% Reset
close all
clearvars

%% Define data paths, get all file names
% Files are assumed to be named so that they are sorted in the same order
% that they were originally generated (i.e. the first file in a
% alphabetically sorted list contains the first part of the scan, etc.).

% Update this path if re-processing on your local computer
dataPath = 'E:\Synaptus\RadarAntarcticaData\Preprocessed\';

%% Create folder for saving output (if not already in place)
imagePath = [dataPath 'Images\'];
if not(exist(imagePath,'dir'))
    mkdir(imagePath)
end

%% Constants
c_medium_orig = 168913914;          % Original wave speed (from "settings")
c_medium = c_medium_orig*1.45;      % Modified wave speed, used in processing, seems better(?)
tDelay = 0;
tMax = 4.6016e-05;
thick = (c_medium/2)*tMax;

%% Get file names
matFiles = dir([dataPath '*.mat']);
matFiles = {matFiles.name};
nFiles = length(matFiles);

%% Process data - concatenate neigboring blocks, use image from middle block

for ii = 1:nFiles
    % Set file indices for neighboring blocks, handle edge cases
    fileInd = [-1 0 1] + ii;
    fileInd = fileInd((fileInd>=1) & (fileInd<=nFiles));
    
    % Concatenate blocks and Fourier transform along x axis
    [Pox,index,omega,xStep] = loadAndConcatenateBlocks(dataPath, matFiles(fileInd));
    [Pokx,kx] = calc_spec_kx(Pox,xStep);
    
    % Create focused image (migration)
    disp(['Creating focused image based on files numbered ' num2str(fileInd)])
    [im,xIm,zIm] = psm_spec_input(Pokx,omega,kx,xStep,tDelay,c_medium,thick);
    
    % Crop image to second block (except for left edge)
    for ll = 1:length(im)
        if ii == 1
            im{ll} = im{ll}(:,index(1,1):index(1,2));   % Only two blocks, keep left
        else
            im{ll} = im{ll}(:,index(2,1):index(2,2));   % Crop to second block
        end
    end

    % Save focused image
    [~,fileBase,~] = fileparts(matFiles{ii});
    imageFileName = [imagePath fileBase '_focused.mat'];
    save(imageFileName,'im','xIm','zIm','c_medium')

end

function [Pox_concat,index,omega,xStep] = loadAndConcatenateBlocks(matFileFolder,matFiles)
% loadAndConcatenateBlocks - load preprocessed files, tile horizontally
%
% Usage:
% [Pox_concat,index,omega,xStep] = loadAndConcatenateBlocks(matFileFolder,matFiles)
%
% Input parameters:
% matFileFolder -   path to folder containing matfiles
% matFiles      -   cell array with names to matfiles containing 'Pox',
%                   'omega', 'xStep' variables
%
% Output parameters:
% Pox_concat    - blocks concatenated horizontally
% index         - indices to start and end of each block (horizontally),
%                 size [nFiles,2]
% omega         - vector of omega values, corresponding to vertical axis of
%                 Pox_concat
% xStep         - spatial sampling size for Pox_concat
% 
% NOTE:
% 'omega' and 'xStep' are loaded from each file, but are overwritten when
% the next file is loaded. Thus, only the variables from the last file are
% kept. These variables are assumed to be identical for each file.
%

%% Status update
disp('Loading preprocessed data blocks.')

%% Get number of files
nFiles = length(matFiles);

%% Load and concatenate
Pox_concat = [];
nX = zeros(nFiles,1);
for ii = 1:nFiles
    load([matFileFolder matFiles{ii}],'Pox','omega','xStep');
    nX(ii) = size(Pox,2);
    Pox_concat = [Pox_concat Pox]; %#ok<AGROW>
end

%% Build index
index = zeros(nFiles,2);
index(:,1) = cumsum([1; nX(1:(end-1))]);
index(:,2) = cumsum(nX);
end


function [Pokx,kx] = calc_spec_kx(Pox,xStep)
% calc_spec_kx - calculate kx spectrum for 2D dataset
%
%   Usage:
%   [Pokx,kx] = calc_spec_kx(Pox,xStep)
%
%   Input parameters:
%   Pox     -    MxN matrix, dimensions (omega,x)
%   xStep   -    spatial sampling step [m]. Assumed to be constant.
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

%% Status update
disp('Fourier transform along x axis.')

%% Create kx vector
nFFTx = 2^nextpow2(size(Pox,2));                    % FFT size x

kxs = (2*pi)/xStep;                                 % Sampling wavenum., x dir.
kx = ((0:(nFFTx-1)) - floor(nFFTx/2))*(kxs/nFFTx);  % X-axis wave number vector
kx = ifftshift(kx);                                 % Shift to correspond to fft() output

%% Fourier transform along time dimension
Pokx = fft(Pox,nFFTx,2);
end

