function [recon, mae] = get_mae_2(signal, co2, Hz, NoWindowsInConceFT, NoConceFT, ...
WindowLength, WindowBandwidth, HighFrequencyLimit, LowFrequencyLimit, ...
FrequencyAxisResolution, HOP, Band, showplot, harmonics)

SamplingRate = Hz;
xm = signal - mean(signal);
time = (1:length(xm))' / Hz;

[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(xm, LowFrequencyLimit, ...
HighFrequencyLimit, FrequencyAxisResolution, HOP, WindowLength, NoWindowsInConceFT, WindowBandwidth, NoConceFT, 0, 0) ;

if showplot
figure ;
imageSQ(time(1:HOP:end), tfrtic*SamplingRate, abs(tfr), .995) ; colormap(1-gray) ; title('spectrogram') ;
axis([-inf inf 0 10]);

figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
axis([-inf inf 0 10])
end

% curve extraction
[c] = CurveExt_M(abs(tfrsq)', .5);


if showplot
figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
hold on;
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate, 'r', 'linewidth', 3);
axis([-inf inf 0 5]);
end

% reconstruction
[h, Dh, t] = hermf(WindowLength, 1, WindowBandwidth) ;

if showplot
figure ;
imageSQ(time(1:HOP:end), tfrsqtic*SamplingRate, abs(tfrsq), .995) ; colormap(1-gray) ; title('SST') ;
hold on;
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate, 'r', 'linewidth',3);
%axis([-inf inf 0 20])
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate+Band, 'b', 'linewidth',3);
plot(time(1:HOP:end), tfrsqtic(c)*SamplingRate-Band, 'b', 'linewidth',3);
end

[recon] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, SamplingRate, c, Band, h((WindowLength+1)/2)) ;


if harmonics > 1
for jj = 2:harmonics
    [tmp] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, SamplingRate, c*jj, Band, h((WindowLength+1)/2)) ;
    recon = recon + tmp ;
end
end

%subtract out heart rate and above
recon = xm - recon;

if showplot
figure;
plot(time, xm);
hold on;
plot(time(1:HOP:end), real(recon),'r');
end


% Compare to true CO2
co2 = co2 - mean(co2);

if showplot
figure;
plot(time, co2 - mean(co2));
hold on;
plot(time(1:HOP:end), real(recon),'r');
end


% MAE
mae = mean(abs(co2(1:HOP:end) - real(recon)'));


end