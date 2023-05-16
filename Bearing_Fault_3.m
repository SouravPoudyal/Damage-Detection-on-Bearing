%Initialization and loading files
clc
clear
close all

fontname='times new roman';
fontsize=12;
linewidth=1.5;


akt_path = pwd;
savepath = fullfile(akt_path,"WS2223_4");
x = what(savepath);
%%
%Loading signal
newload = x.mat(3);
signalname = newload;
load(fullfile(savepath, signalname));
D_3 = datensatz;
[n,m] = size(D_3);

figure(1)
plot(D_3,'-k','linewidth',linewidth,'displayname',['|{\itx}_3|'])
xlabel('Time in s')
ylabel('{\itx}_3 (m/s^2)')
set(gca,'fontname',fontname,'fontsize',fontsize)
title('Signal_3')
legend toggle
%%
%calculation for Damaged frequencies for a Bearing
nwx = 13;
Dw = 3.7;
Dt = 26.15;
a = 0;
fk = fn/2*(1-(Dw/Dt)*cos(a));
fk_ = fn/2*(1+(Dw/Dt)*cos(a));
fa = fn*nwx/2*(1-(Dw/Dt)*cos(a));
fi = fn*nwx/2*(1+(Dw/Dt)*cos(a));
fwa = fn*nwx/2*(Dt/Dw)*(1-(Dw/Dt*cos(a))^2);
fw = 2*fwa;

%%
%using fast fourier transform to plot the signal in frequency domain

X=fft(D_3);
df=fs/m; %resolution in Hz
f_achs=[0:df:fs-df]';

figure(2)

plot(f_achs,sqrt(conj(X).*X)*df/(fs/2),'-k','linewidth',linewidth,'markersize',10)
xlabel('Frequency in Hz')
title('Accleration Spectrum')
ylabel('|{\itA-x3}({\itf}) | in m/s^2')
set(gca,'xlim',[0 fs/2])
set(gca,'ylim',[0 4])
set(gca,'fontname',fontname,'fontsize',fontsize)
legend('A-x3')
%%
%Plotting the Envelop specturm of the signal

ht_D=hcurve_fun(D_3, 0, m); %Find envelop of the signal using hilbert transfrom
figure(3)

X=fft(ht_D); % finding fft
df=fs/m; %resolution in Hz
f_achs=[0:df:fs-df]';
plot(f_achs,sqrt(conj(X).*X)*df/(fs/2),'-k','linewidth',linewidth,'markersize',5)
xlabel('Frequency in Hz')
title('Accleration Envelop Spectrum')
ylabel('|{\itA-x3}({\itf}) | in m/s^2')
set(gca,'xlim',[0 fs/20])
set(gca,'ylim',[0 4])
set(gca,'fontname',fontname,'fontsize',fontsize)
legend('A-x3')
%%
%Ploting the roll over frequency harmonics
harmImpact = (0:1:20)*fa;
[X,Y] = meshgrid(harmImpact,ylim);

hold on
plot(X,Y,'g', 'DisplayName', 'fa harmonics')
legend('A-x3','fa-harmonics')
hold off
%%
%Finding the Short Time Fourier Transform of the Signal
figure(4)
nfft = 4463;
ht_D=hcurve_fun(D_3, 0, m);
w = hamming(nfft);
[spec_,f_,t_] = stft(ht_D,fs, Window=w,OverlapLength=2231,FFTLength=nfft,FrequencyRange="onesided");

imagesc(t_, f_, 10*log10(abs(spec_)));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
a = colorbar;
a.Label.String = 'Power (dB)';
%%
%Ploting the Waterfall plot of the signal
figure(5)
waterplot(spec_,f_,t_)
title("stft")
a = colorbar;
a.Label.String = 'Power (dB)';

%%
%Custome code to generate Envelop Power Spectral density
% Reference: https://www.dsprelated.com/showarticle/1221.php

% Given signal for example x(1:8192),
%
%            |x(1)    x(1025) . . . . x(7169)|
%            |x(2)    x(1026) . . . . x(7170)|
% xMatrix=   | .         .               .   |
%            | .         .               .   |
%            |x(1024) x(2048) . . . . x(8192)|
%
%

%with 50% overlap
D_ = D_3(1:1620060);
D_o = D_(806:1619253);
nfft= 1612;                      % number of frequency samples
Navg=2009;



% Createing a matrix xMatrix with Navg columns,
% each column a segment of x of length nfft


xMatrix= reshape(D_,[],1005);
xMatrix_o= reshape(D_o, 1612, []);
xMatrix= [xMatrix xMatrix_o];
window= hamming(1612); 
magsq_sum= zeros(nfft/2);
for i= 1:Navg
    x_column= xMatrix(:,i);
    xw= 2.*x_column.*window;
    ht_D=hcurve_fun(xw, 0, 1612);        % applying window of length nfft
    X= fft(ht_D);                        % DFT
    X= X(1:nfft/2);                      % retaining samples from 0 to fs/2
    magsq= real(X).^2 + imag(X).^2;      % DFT magnitude squared
    magsq_sum= magsq_sum + magsq;        % sum of DFT mag squared
end
mag_sq_avg= magsq_sum/Navg;              % average of DFT mag squared
P_bin= 2/nfft.^2 *mag_sq_avg;            % W/bin power spectrum
P_Hz= P_bin*nfft/fs;                     % W/Hz power spectrum
PdB_bin= 10*log10(P_bin);                % dBW/bin
PdB_Hz= 10*log10(P_Hz);                  % dBW/Hz
k= 0:nfft/2 -1;                          % frequency index
f= k*fs/nfft;                            % Hz frequency vector

figure(6)
h = plot(f,PdB_bin),grid
axis([0 fs/2 -60 10])
xlabel('Hz'),ylabel('dBW/bin')
title('Envelop PSD in dBW/bin')
legend([h(1)],'PSD-x3')

figure(7)
h1 = plot(f,PdB_Hz);
axis([0 fs/2 -60 10]);

xlabel('Hz'),ylabel('dBW/Hz')
title('Envelop PSD in dBW/Hz')
hold on
xlim([0 20*fa])
%ylim([-40 10])
harmImpact = (0:1:50)*fa;
[X,Y] = meshgrid(harmImpact,ylim);
h2 = plot(X,Y,'g');
legend([h1(1),h2(1)],'PSD-x3','fa-harmonics')
hold off

%%
%Calculating parameters for a kutogram
[kgram, f, w, fc, wc, bw] = kurtogram(D_3, fs);
%%
%Plotting the Kutrogram of the Signal
figure(7)
kurtogram(D_3,fs)
%%
%Filtering the Signal using Spectral Kutrosis
fc = fc/fs;
b = hanning(wc)';
b = b/sum(bw);
b = b.*exp(2i*pi*(0:wc-1)*fc);
xfilt = fftfilt(b,D_3);

%%
%Plotting the filtered Signal
figure(8)
plot(abs(xfilt),'-k')

%%
function waterplot(s,f,t)
% Waterfall plot of spectrogram
    waterfall(f,t,10*log10(abs(s)'.^2))
    set(gca,View=[30 50])
    xlabel("Frequency in Hz")
    ylabel("Time in S")
end
