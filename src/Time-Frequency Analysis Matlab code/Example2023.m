% 2023 03 09
%
% in this simulation example, how to run spectrogram
% SST and ConceFT is demonstrated
% The same flow could be appied to conceFT for the synchrosqueezed CWT or others.

%%
pwd
%%
clear ; close all ;
addpath('./tool') ;
addpath('./Morse') ;
%cd 'Time-Frequency Analysis Matlab code'/;

if 0
    cflow = csvread('CFlow.csv');
    flow = highpass(x(:,2)-mean(x(:,2)), 0.1, 100) ;
    flow = resample(flow, 64, 100) ;
    xm = flow(33001: 33000+4096) ;
    Hz = 100 ;
    time = [1:length(xm)]' / Hz ;
end


if 0
    load BP.mat ;
    xm = resample(bp, 100, 1000) ;
    xm = xm' - mean(xm) ;
    Hz = 100 ;
    time = [1:length(xm)]' / Hz ;

    % try ARMA(1,1)
    snrdb = 5 ; 
    sigma = exp(-snrdb/20) * std(xm) ;
    Mdl = arima('Constant',0.5,'AR',{0.7 0.25},'Variance',.1);
    noise = simulate(Mdl, length(xm)/2);
    %noise = armaxfilter_simulate(length(xm)/2, 0, 1, .95, 1, -.5) ;
    xm(1:end/2) = xm(1:end/2) + sigma * noise ./ std(noise) ;
end

if 1
    load 0009_8min.mat
    Hz = 300;
    xm = signal.pleth.y - mean(signal.pleth.y);
    time = [1:length(xm)]' / Hz;

end
%%

if 0
    load ta-mECG.mat
    mecg = resample(mecg-mean(mecg), 100, 1000) ;
    trend = zeros(size(mecg)) ;
    for jj = 1:length(trend)
        idx = max(1, jj-25) : min(jj+25, length(trend)) ;
        trend(jj) = median(mecg(idx)) ;
    end
    trend = smooth(trend, 100, 'loess') ;
    xm = mecg - trend ;
    Hz = 100 ;
    time = [1:length(xm)]' / Hz ;
end

% simulated signal
if 0
    % the sampling rate for the simulated signal
    Hz = 100 ;
    % the sampling time of the simulated signal
    time = [1/Hz:1/Hz:17]' ;
    % the number of sampling points of the simulated signal
    N = length(time) ;

    % fix the random seed for the reproducibility issue
    initstate(2) ;


    % the amplitude modulation of the simulated signal
    % simulate 1 oscillatory components with dynamics
    am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
    am1 = 2 + am1 ./ max(abs(am1)) ;

    % the instantaneous frequency of the simulated signal
    if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
    if1 = 4 + 1 * if1 ./ max(abs(if1)) ;
    phi1 = cumsum(if1) / Hz + time.^2 / 10 ;
    if1 = if1 + time / 5 ;

    % the simulated signal.
    % this one is very interesting! get very clear double bumps
    %s1 = am1 .* cos(2*pi*phi1) + max(2, am1.^2) .* cos(2*pi*2*phi1+2) ...
    %   + min(3, am1.^2) .* cos(2*pi*3*phi1+1) ;

    s1 = am1 .* cos(2*pi*phi1) ; %+ max(2, am1.^2) .* cos(2*pi*2*phi1+1.1) ...
    %+ min(3, am1.^2) .* cos(2*pi*3*phi1+1) ;
    s2 = am1 .* exp(i*2*pi*phi1) + max(2, am1.^2) .* exp(i*2*pi*2*phi1+1.1) ...
        + min(3, am1.^2) .* exp(i*2*pi*3*phi1+1) ;


    %s1 = mecg ;
    clean = s1 ;
    %clean = cos(2*pi*44*(time + 0.0 * time.^2)) + 0*cos(2*pi*2*(time + 0.025 * time.^2)) ...
    %    + 0*cos(2*pi*3*(time + 0.025 * time.^2)) ;

    % add noise
    sigma = 0.5 ;
    noise = random('T',4,N,1) ;
    noise = sigma * noise ;
    var(noise)
    snrdb = 20 * log10(std(clean)./std(noise)) ;
    fprintf(['snrdb = ',num2str(snrdb),'\n']) ;

    % simulated observed time series
    xm = clean + noise ;

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
WindowBandwidth = 10 ;
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

HOP = 10 ;

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




keyboard

%% curve extraction
[c] = CurveExt_M(abs(tfrsq)', .5);

figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
hold on;
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate, 'r', 'linewidth',3);
axis([-inf inf 0 20])





%0.3 Hz for resp rate
% reconstruction
[h, Dh, t] = hermf(WindowLength, 1, WindowBandwidth) ;
Band = 0.2 ;

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
%%
recon2 = resample(recon, 100, 100/HOP) ;
plot(time, real(recon2),'b', 'linewidth', 2)



% if more harmonics
if 0
    % in this example, if the input is BP signal, consider the first 5 harmonics
    % for a comparison with the upcoming SAMD algorithm
    for jj = 2:5
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
