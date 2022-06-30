signal = 'blow';
fs = 9600;
freqLimits = [1 4800];
bpo = 3;

% Load Signal
s = load(signal);
x = s.(signal);
t = (0:length(x)-1)*1/fs;

% Bandpass Filtering
FftFilterBank = fftBankFilterDesign(fs,bpo,freqLimits);
[xb_fft,fc] = fftBankFilter(FftFilterBank,x,false);
[xb_fft_zp,~] = fftBankFilter(FftFilterBank,x,true);
Xb = 20*log10(xb_fft(1,:));
Xb_zp = 20*log10(xb_fft_zp(1,:));

% Plotting
figure
hold on
plot(fc,Xb,'b')
plot(fc,Xb_zp,'m')
xlabel('Central Frequency [Hz]')
ylabel('Band Level [dBV]')
title(sprintf('FFT Bandpass Filter Bank \\rm(%s)',signal))
set(gca,'XScale','log')
box on
axis([fc(1) fc(end) -10 170])
legend({'Without Zero Padding','With Zero Padding'},'Location','SouthEast')
fname = sprintf ('FFT BPF (%s)',signal);
savefig(fname)
print(fname,'-dpng','-r200')
