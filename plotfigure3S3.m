function plotfigure3

myPath = fileparts(mfilename('fullpath'));
addpath(fullfile(myPath,'functions'));
addpath(fullfile(myPath,'data'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gi = @(k,a1,a2) 1e3/(a1-a2)*((1e3/a1+sqrt(-1)*2*pi*k).^(-1)-(1e3/a2+sqrt(-1)*2*pi*k).^(-1));
lgNorm = @(k,m1,s1) exp(m1-s1^2/2)*exp(-(log(k)-m1).^2/(2*s1^2))./k;
gausses = @(m1,s1,b1,m2,s2,b2,m4,s4,b4,k) b1*lgNorm(k,m1,s1)+b2*lgNorm(k,m2,s2)+b4*lgNorm(k,m4,s4);
gz = @(k,a1,n) 1./(1+2*pi*sqrt(-1)*(k/a1).^n);
ftfun = @(m1,s1,m2,s2,m4,s4,b1,b2,b4,Z1,gamI,tauI,tauI2,wf,theta,f) 2* log(abs(gz(f,wf,theta))) +  ...
	log( Z1 + gamI*(1+gamI/2*gausses(m1,s1,b1,m2,s2,b2,m4,s4,b4,f)).* abs(gi(f,tauI,tauI2)).^2);
F = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),f);
FBL = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),0,0,0,x(10),x(11),x(12),x(13),x(14),x(15),f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('sampleTimeSeries2.mat')
load('ModelParameters.mat');
load('pt12SegFit.mat');
c = pt12Params; c(14) = Inf;
c1 = c; c1(10) = 0;
c2 = c; c2(11) = 0;

dt = 1./nanmean(diff(time));
t = 0:1/dt:20;

gam = @(t1,t2) (1-t1/t2)*(t1/t2)^(-t1/(t1-t2));
TI1 = c(13)*1e-3;
TI2 = c(12)*1e-3;
aI = @(t) (exp(-t/TI1)-exp(-t/TI2)).*(t>=0);
icov  = aI(0:1/dt:1);
icov = icov/sum(icov);

Wi = conv(randn(length(t)+length(icov)-1,1)*sqrt(c(11)*dt),icov,'valid');
Wa = randn(length(t),1)*sqrt(c(10)*dt);

[b,a] = butter(2,0.1/dt*2,'high');
[b2,a2] = butter(2,300/dt*2,'low');

ln = 0.8*sin(2*pi*t*60)'+0.3*sin(2*pi*t*180)'+0.3*sin(2*pi*t*300)'+0.3*sin(2*pi*t*420)';
Xsim = filtfilt(b2,a2,filtfilt(b,a,Wi+Wa+ln));
Xsim = Xsim(and(t>5,t<=15));

[freq,~,psd] = eegfft(time,pt12EEG,2,1.9);
[freq,~,psd2] = eegfft(t,Xsim,2,1.9);


m = 3;
fig = figure('color','w','units','centimeters');
fig.Position = [0,0,18.3,8];
subplot(2,4,1:3);
	plot(downsample(time,m),downsample(pt12EEG,m),'color','k');
	xlim([-170,-160]);
	ylim([-50,50]);
	axis off;
	line([-170,-169],[-50,-50],'color','k','linewidth',1.5)
	line([-170,-170],[-50,-30],'color','k','linewidth',1.5)
	text(-170.1,-40,['20 ' char(956) 'V'],'HorizontalAlignment','center','VerticalAlignment','bottom','rotation',90,'FontSize',7)
	text(-169.5,-52,'1 s','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7)

ax=subplot(2,4,4);
[xData, yData] = prepareCurveData( freq, log(nanmean(psd,2)) );
	plot(xData,exp(yData),'k','LineWidth',1); hold on;
	plot(freq,exp(FBL(c,freq)),'r');
	plot(freq,exp(FBL(c1,freq)),'b','LineWidth',1);
	plot(freq,exp(FBL(c2,freq)),'g','LineWidth',1);
	set(gca,'yscale','log')
	set(gca,'xscale','log')
	xlim([0.5,400]);
	xticks([0.01,0.1,1,10,100]);
	xticklabels([0.01,0.1,1,10,100]);
	hold on;
	ylim([1e-4,1e4]);
	yticks([1e-4,1e-2,1,1e2,1e4]);
	ylabel(['PSD (' char(956) 'V^2/Hz)'])
	xlabel('Frequency (Hz)')
	gcaformat;
	grid on;
	ax.YMinorGrid = 'off';
	ax.MinorGridLineStyle = '-';
	ax.MinorGridColor = [0.8,0.8,0.8];
	ax.GridColor = [0,0,0]; ax.GridAlpha = 0.3;
	ax.XAxis.Color = [0.5,0.5,0.5]; ax.XAxis.Label.Color = 'k';
	ax.YAxis.Color = [0.5,0.5,0.5]; ax.YAxis.Label.Color = 'k';
	box on
	set(gca,'Tickdir','in')



subplot(2,4,5:7);
	plot(downsample(time,m),downsample(Xsim,m),'color','k');
	xlim([-170,-160]);
	ylim([-50,50]);
	axis off;
ax=subplot(2,4,8);
[xData, yData] = prepareCurveData( freq, log(nanmean(psd2,2)) );
	plot(xData,exp(yData),'k','LineWidth',1); hold on;
	plot(freq,exp(FBL(c,freq)),'r');
	plot(freq,exp(FBL(c1,freq)),'b','LineWidth',1);
	plot(freq,exp(FBL(c2,freq)),'g','LineWidth',1);
	set(gca,'yscale','log')
	set(gca,'xscale','log')
	xlim([0.5,400]);
	xticks([0.01,0.1,1,10,100]);
	xticklabels([0.01,0.1,1,10,100]);
	hold on;
	ylim([1e-4,1e4]);
	yticks([1e-4,1e-2,1,1e2,1e4]);
	ylabel(['PSD (' char(956) 'V^2/Hz)'])
	xlabel('Frequency (Hz)')
	gcaformat;
	grid on;
	ax.YMinorGrid = 'off';
	ax.MinorGridLineStyle = '-';
	ax.MinorGridColor = [0.8,0.8,0.8];
	ax.GridColor = [0,0,0]; ax.GridAlpha = 0.3;
	ax.XAxis.Color = [0.5,0.5,0.5]; ax.XAxis.Label.Color = 'k';
	ax.YAxis.Color = [0.5,0.5,0.5]; ax.YAxis.Label.Color = 'k';
	box on
	set(gca,'Tickdir','in')
