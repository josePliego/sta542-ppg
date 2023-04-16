% 2023 03 09
%
% in this simulation example, how to run spectrogram
% SST and ConceFT is demonstrated
% The same flow could be appied to conceFT for the synchrosqueezed CWT or others.

%%
pwd
clear ; close all ;
addpath('./tool') ;
addpath('./Morse') ;
%cd 'Time-Frequency Analysis Matlab code'/;

if 1
    load ../../data/0009_8min.mat
    Hz = 300;
    xm = signal.pleth.y - mean(signal.pleth.y);
    time = (1:length(xm))' / Hz;

end

%% setup parameters for the SST or ConceFT

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
WindowBandwidth = 14 ;
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

% call the main code, which is the ConceFT based on
% synchrosqueezed short time Fourier transform (STFT)
% Output:
% tfr: STFT result
% tfrtic: frequency axis tic for the STFT
% tfrsq: synchrosqueezed STFT (it is equivalent to running ConceFT only one time)
% ConceFT: ConceFT of synchrosqueezed STFT.
% tfrsqtic: frequency axis tic for the tfrsq and ConceFT
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(xm, LowFrequencyLimit, ...
    HighFrequencyLimit, FrequencyAxisResolution, HOP, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, 0, 0) ;

% plot the time frequency representation determined by
% ConceFT. .995 is the quantile truncation set to avoid
% possible outliers in the final analysis result.
% see the code for details.
figure ;
imageSQ(time(1:HOP:end), tfrtic*SamplingRate, abs(tfr), .995) ; colormap(1-gray) ; title('spectrogram') ;
axis([-inf inf 0 10])

figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
axis([-inf inf 0 10])






%% curve extraction
[c] = CurveExt_M(abs(tfrsq(1:30, :))', 0.5);
%[c] = CurveExt_M(abs(tfrsq)', .5);
figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
hold on;
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate, 'r', 'linewidth', 3);
axis([-inf inf 0 5])





%0.3 Hz for resp rate
%% reconstruction
[h, Dh, t] = hermf(WindowLength, 1, WindowBandwidth) ;
Band = 0.02 ;

figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
hold on;
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate, 'r', 'linewidth',3);
axis([-inf inf 0 20])
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate+Band, 'b', 'linewidth',3);
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate-Band, 'b', 'linewidth',3);


[recon] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, SamplingRate, c, Band, h((WindowLength+1)/2)) ;

figure;
plot(time, xm)
hold on;
plot(time(1:HOP:end), real(recon),'r')

%% Compare to true CO2
co2 = signal.co2.y;
co2 = co2 - mean(co2);

figure;
plot(time, co2 - mean(co2));
hold on;
plot(time(1:HOP:end), real(recon),'r');

% MAE
mean(abs(co2(1:HOP:end) - real(recon)'));

%%
%recon2 = resample(recon, 100, 100/HOP) ;
%plot(time, real(recon2),'b', 'linewidth', 2)



% if more harmonics
if 1
    % in this example, if the input is BP signal, consider the first 5 harmonics
    % for a comparison with the upcoming SAMD algorithm
    for jj = 2:2
        [tmp] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, SamplingRate, c*jj, Band, h((WindowLength+1)/2)) ;
        recon = recon + tmp ;
    end

    figure;
    plot(time, xm)
    hold on;
    plot(time(1:HOP:end), real(recon),'r')
    axis([110 120 -inf inf])

    figure ;
    imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
    hold on;
    plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate, 'r', 'linewidth',3);
    axis([110 120 0 10])
end




% get CWT.
[tfr, tfrtic] = CWT2023(xm) ;

subplot(121)
imagesc(time, 1:length(tfrtic), abs(tfr')); axis xy;
colormap(1-gray) ;
% number of ticks in the y-axis
yN = 15 ;

% it is log-linear scale in the y-axis
set(gca,'YTick',floor(length(tfrtic)/yN):floor(length(tfrtic)/yN):length(tfrtic));
scaletic = (time(end)-time(1))./tfrtic ;
set(gca,'YTickLabel',round(100*scaletic(get(gca,'ytick')))/100);

xlabel('Time (sec)') ; title('scalogram')
ylabel('Scale (sec)') ;
set(gca,'fontsize', 24) ;


% rescale the y-axis to the linear tick
freqmin = 1./scaletic(end) ;
freqmax = 1./scaletic(1) ;
freqtic = linspace(freqmin, freqmax, 1000) ;
tfr2 = zeros(size(tfr,1), length(freqtic)) ;
for jj = 1: size(tfr2,1)
    tmp = interp1(1./scaletic, abs(tfr(jj,:)), freqtic, 'linear') ;
    tfr2(jj, :) = tmp ;
end

subplot(122)
imageSQ(time, freqtic, abs(tfr2'), .999) ;
colormap(1-gray) ;
xlabel('Time (sec)') ; title('TFR of CWT')
ylabel('Frequency (Hz)') ;
set(gca,'fontsize', 24) ;
