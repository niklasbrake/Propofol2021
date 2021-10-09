function [freq,time,psd,FFT] = eegfft(t,X,lWindow,overlap);

	Fs = 1/mean(diff(t)); % sampling frequency
	T = lWindow; % window length in seconds
	iwindow = -floor(Fs*T/2)+1:floor(T/2*Fs);
	L = length(iwindow);
	f = Fs * linspace(0,1/2,L/2+1);
	f = f(2:end);
	minT = floor(Fs*T/2);
	maxT = size(X,1)-floor(Fs*T/2);
	difT = ceil(Fs*(T-overlap));
	rachTs = minT:difT:maxT;
	hmw = hamming(L);
	K = zeros(L/2,length(rachTs));
	thet = K;
	for i = 1:length(rachTs);
		temp = X(rachTs(i)+iwindow);
		if(any(isnan(temp)))
			K(:,i) = nan;
			continue;
		end
		temp = temp.*hmw;
		temp2 = fft(temp);
		thet(:,i) = angle(temp2(2:L/2+1));		
		P = abs(temp2).^2*(1/Fs).^2/T;
		K(:,i) = 2*P(2:L/2+1);
		FFT(:,i) = temp2;
	end
	freq = f(:);
	time = t(rachTs(:));
	psd = K;