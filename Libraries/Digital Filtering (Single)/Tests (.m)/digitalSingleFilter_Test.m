tau = 1;
fs = 9600;
freqCutoff = [50 200];
filtMode = 'filter';
filtType = 'bandpass';
t = 0:1/fs:tau-1/fs;
nSamples = length(t);
x = randn(nSamples,1);

DigitalFilter = digitalSingleFilterDesign(fs,freqCutoff,'FilterOrder',10,...
    'FilterType',filtType);
xf = digitalSingleFilter(DigitalFilter,x,'MetricsOutput',0,...
    'FilterMode',filtMode,'ZeroPadding',true,'DataWrap',true); 
            
[X,f] = periodogram(x,hamming(length(x)),8192,fs,'psd','onesided');
[Xf,~] = periodogram(xf,hamming(length(xf)),8192,fs,'psd','onesided');

% Plotting
figure
hold on
plot(f,10*log10(X),'b')
plot(f,10*log10(Xf),'m')
xlabel('Frequency [Hz]')
ylabel('PSD [dBV]')
title(sprintf('Bandpass Filter \\rm(wnoise, %0.0f-%0.0f Hz)',freqCutoff(1),freqCutoff(2)))
legend({'Original','Filtered'},'Location','SouthEast')
set(gca,'XScale','log')
box on
axis([f(1) f(end) -80 -20])
fname = 'BPF (wnoise)';
savefig(fname)
print(fname,'-dpng','-r200')
