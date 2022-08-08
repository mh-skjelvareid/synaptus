%   fix_bad_samples - replace bad samples with the mean of their neighbors
%
%   2021-09-16  MHS - Initial version
%

%% Reset
clearvars

%% Path
% Update this path if re-processing on your local computer
preprocessed_path = 'E:\Synaptus\RadarAntarcticaData\Preprocessed\';

%% Fix bad sample in file 2
load([preprocessed_path '20191204_030602_VHF_Radar_Antarctica_2019_241-260_combined_preprocessed.mat'],...
    'Pox','xStep','omega');

% Column 1439 of file 2 (downsampled) is corrupted. Set to mean of neighbors.
Pox(:,1439) = mean(Pox(:,[1438 1440]),2);

save([preprocessed_path '20191204_030602_VHF_Radar_Antarctica_2019_241-260_combined_preprocessed.mat'],...
    'Pox','xStep','omega');

%% Fix bad sample in file 3
load([preprocessed_path '20191204_030602_VHF_Radar_Antarctica_2019_261-280_combined_preprocessed.mat'],...
    'Pox','xStep','omega');

% Column 1371 of file 3 (downsampled) is corrupted. Set to mean of neighbors.
Pox(:,1371) = mean(Pox(:,[1370 1372]),2);

save([preprocessed_path '20191204_030602_VHF_Radar_Antarctica_2019_261-280_combined_preprocessed.mat'],...
    'Pox','xStep','omega');
