function [tfr, scaletic] = CWT2023(x, opts)
%
% Continuous wavelet transform ver 0.1
% Modified from wavelab 85
%
% by Hau-tieng Wu 2023-02-14
%
%	To plot the scalogram, remember that scaletic is the inverse of
%	the scale a, and tfr is evaluated in the log scale.
%	be careful of the imagesc usage of the log scale tick
%	imagesc
%	where
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% prepare for the input

if nargin < 3
    opts = struct();
    opts.motherwavelet = 'Cinfc' ;
    opts.CENTER = 1 ;
    opts.FWHM = 0.3 ;
end

nvoice = 32;
scale = 1; % the starting scale
oct = 1;




%% start to do CWT

n = length(x);
% assume the original signal is on [0,L].
% scale the signal to be on [0,1].
% Frequencies are rescaled to xi*L
xi = [(0:(n/2)) (((-n/2)+1):-1)]; xi = xi(:);
xhat = fft(x);

noctave = floor(log2(n)) - oct;
tfr = zeros(n, nvoice .* noctave);

kscale = 1;
%scaletic = scale .* 2.^((1 : nvoice .* noctave)/nvoice) ;
scaletic = 2.^((1 : nvoice .* noctave)/nvoice) ;


% scaletic = zeros(1, nvoice .* noctave) ;
% for jj = 1 : nvoice .* noctave
%     scaletic(jj) = scale .* (2^(jj/nvoice));
% end



for jo = 1:noctave	% # of scales
    for jv = 1:nvoice

        qscale = scale .* (2^(jv/nvoice));
        omega =  xi ./ qscale ;


        if strcmp(opts.motherwavelet, 'Cinfc');

            tmp0 = (omega-opts.CENTER)./opts.FWHM;
            tmp1 = (tmp0).^2-1;

            windowq = exp( 1./tmp1 );
            windowq( find( omega >= (opts.CENTER+opts.FWHM) ) ) = 0;
            windowq( find( omega <= (opts.CENTER-opts.FWHM) ) ) = 0;

        elseif strcmp(opts.motherwavelet, 'morlet')

            windowq = 4*sqrt(pi)*exp(-4*(omega-0.69*pi).^2)-4.89098d-4*4*sqrt(pi)*exp(-4*omega.^2);

        elseif strcmp(opts.motherwavelet, 'gaussian');

            psihat = @(f) exp( -log(2)*( 2*(f-opts.CENTER)./opts.FWHM ).^2 );
            windowq = psihat(omega);


        elseif strcmp(opts.motherwavelet, 'meyer');  %% Meyer

            windowq = zeros(size(omega));
            int1 = find((omega>=5./8*0.69*pi)&(omega<0.69*pi));
            int2 = find((omega>=0.69*pi)&(omega<7./4*0.69*pi));
            windowq(int1) = sin(pi/2*meyeraux((omega(int1)-5./8*0.69*pi)/(3./8*0.69*pi)));
            windowq(int2) = cos(pi/2*meyeraux((omega(int2)-0.69*pi)/(3./4*0.69*pi)));

        elseif strcmp(opts.motherwavelet, 'BL3');    %% B-L 3

            phihat = (2*pi)^(-0.5)*(sin(omega/4)./(omega/4)).^4; phihat(1) = (2*pi)^(-0.5);
            aux1 = 151./315 + 397./840*cos(omega/2) + 1./21*cos(omega) + 1./2520*cos(3*omega/2);
            phisharphat = phihat.*(aux1.^(-0.5));

            aux2 = 151./315 - 397./840*cos(omega/2) + 1./21*cos(omega) - 1./2520*cos(3*omega/2);
            aux3 = 151./315 + 397./840*cos(omega) + 1./21*cos(2*omega) + 1./2520*cos(3*omega);
            msharphat = sin(omega/4).^4.*(aux2.^(0.5)).*(aux3.^(-0.5));
            windowq = phisharphat.*msharphat.*exp(i*omega/2).*(omega>=0);

    	end

        windowq = windowq ./ sqrt(qscale);

        what = windowq .* xhat;

        w = ifft(what);
        tfr(:,kscale) = transpose(w);
        kscale = kscale+1;

    end

    scale = scale .* 2;

end


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% calculate the constant for reconstruction (very rough estimate. Needs to be improved)
xi = [0.05:1/10000:10];

if strcmp(opts.motherwavelet, 'gaussian');   %% Gaussian (not really wavelet)

    psihat = @(f) exp( -log(2)*( 2*(f-opts.CENTER)./opts.FWHM ).^2 );
    windowq = psihat(xi);
    Rpsi = sum(windowq./xi)/10000;

elseif strcmp(opts.motherwavelet, 'morlet');

    windowq = 4*sqrt(pi)*exp(-4*(xi-0.69*pi).^2)-4.89098d-4*4*sqrt(pi)*exp(-4*xi.^2);
    Rpsi = sum(windowq./xi)/10000;

elseif strcmp(opts.motherwavelet, 'Cinfc');

    tmp0 = (xi - opts.CENTER)./opts.FWHM;
    tmp1 = (tmp0).^2-1;

    windowq = exp( 1./tmp1 );
    windowq( find( xi >= (opts.CENTER+opts.FWHM) ) ) = 0;
    windowq( find( xi <= (opts.CENTER-opts.FWHM) ) ) = 0;
    Rpsi = sum(windowq./xi)/10000;

else

    % Normalization is not implemented for Other mother wavelets
    Rpsi = 1 ;

end


tfr = tfr ./ Rpsi;


end
