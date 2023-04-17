% Calculate MAE for all signals
% Also saves reconstructions in csv

clear ; close all ;
addpath('./tool') ;
addpath('./Morse') ;

files = dir('../../data/0*.mat');


Hz = 300;
% number of chosen orthonormal windows for ConceFT
NoWindowsInConceFT = 1 ;
% number of random linear combinations of chosen windows
NoConceFT = 1 ;
% the window length. Ideally, it should be chosen so that
% roughly 7-10 oscillations (ignore the multiples) are
% included in the window.
WindowLength = Hz * 10 + 1 ; %377 ;
% this is the bandwith of the chosen window. See hermf.m
% in the attached code for details.
WindowBandwidth = 14;
SamplingRate = Hz ;
% Setup the frequency range for the analysis
% The allowed range is 0-0.5
% This part might be tricky. This 0-0.5 limitation is
% setup under the assumption that the sampling rate is 1Hz
% After the analysis, the sampling rate should be adjusted
% so that the result is associated with the original
% sampling rate.
% In this example, the true range is [0, 0.5]*SamplingRate
HighFrequencyLimit = 0.1 ;
LowFrequencyLimit = 0 ;
% the frequency axis resolution in the final time-frequency representation
FrequencyAxisResolution = 1e-4 ;

HOP = 10;

Band = 0.02;

%% MAE 2

mae_vec = zeros(1, length(files));

recon_mat = zeros(14401, length(files));
path = "../../fig/";

for i = 1:1
    load(strcat(files(i).folder, '/', files(i).name))
    strcat('Reconstructing :', files(i).name)

    [recon, mae] = get_mae_2(signal.pleth.y, signal.co2.y, Hz, NoWindowsInConceFT, NoConceFT, ...
WindowLength, WindowBandwidth, HighFrequencyLimit, LowFrequencyLimit, ...
FrequencyAxisResolution, HOP, Band, 0, 10);
    %mae_vec(i) = mae;
    %recon_mat(:, i) = real(recon)';
end

%co2 = signal.co2.y - mean(signal.co2.y);
%mean(abs(co2(1:HOP:end) - recon))

%SamplingRate = Hz;
xm = signal.pleth.y - mean(signal.pleth.y);
time = (1:length(xm))' / Hz;

figure;
plot(time, signal.co2.y - mean(signal.co2.y));
hold on;
plot(time(1:HOP:end), recon - mean(recon),'r');


% names = {files.name};
% mae_tab = array2table(mae_vec);
% mae_tab.Properties.VariableNames(:) = names;
% writetable(mae_tab,'../../data/mae_reconstruction.csv');
% 
% signal_tab = array2table(recon_mat);
% signal_tab.Properties.VariableNames(:) = names;
% writetable(signal_tab,'../../data/signal_reconstruction.csv');