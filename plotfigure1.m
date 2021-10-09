function plotfigure1

myPath = fileparts(mfilename('fullpath'));
addpath(fullfile(myPath,'functions'));
addpath(fullfile(myPath,'data'));


	
fig = figure('color','w','units','centimeters');
fig.Position = [0,0,12.5,11];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% \begin{Example Time Series} %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('sampleTimeSeries1.mat');
axes('Position',[0.076,0.74,0.85,0.22]);
	plot(downsample(time,5),downsample(timedomain,5),'LineWidth',0.2); hold on;
	xlim([time(1)-10,time(end)]);
	ylim([-100,100]);
	line([time(1)-10,time(1)-10],[-25,25],'color','k','linewidth',1);
	line([time(1),time(1)+60],[-60,-60],'color','k','linewidth',1);
	text(time(1)-15,0,['50 ' char(956) 'V'],'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',7,'color','k','Rotation',90);
	text(time(1)+30,-65,'60 s','VerticalAlignment','top','HorizontalAlignment','center','fontsize',7,'color','k');
	scatter(0,90,10,'vk','filled'); text(-5,95,sprintf('Object dropped (LOC)'),'FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','left','color','k');
	line([-200,-1],[75,75],'color','r','linewidth',2);
	text(-100,80,'Propofol Infusion','color','r','FontSize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
	axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% \end{Example Time Series} %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% \begin{Baseline-normalized analysis} %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('timeInformation.mat','timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;

% Compute spectrogram
%{
% load('FlaviePtsDropAlignedTimeDomain20191012.mat');
for i = 1:14
	[freq,time,psd(:,:,i)] = eegfft(Time,squeeze(TimeDomainAligned(:,2,i)),2,1.9);
end
%}
% Load presaved results
% load('PSD_2s_1.9overlap.mat),'freq','time','psd');
load(fullfile('psd_channel_Cz.mat'),'freq','time','psd');

% Compute baseline spectrum
for i = 1:14
	BL(:,i) = nanmedian(psd(:,time<infusionTime(i),i),2);
end

% Convert to baseline-normalized decibels
dB = 10*log(psd./permute(BL,[1,3,2]))/log(10);
medDB = smoothdata(dB,1,'movmean',3);

% Plot spectrogram heatmap
axes('Position',[0.1,0.44,0.36,0.28]);
	imagesc(time,freq,nanmedian(medDB,3));
	axis xy
	ylim([0.5,40])
	CM = colorbar;
	CM.Position = [0.4687 0.4846 0.0282 0.1867];
	% CM.Title.String = '(dB)';
	set(gca,'CLim',[-1,15])
	% set(gca,'CLim',[0,1])
	colormap(jet(1e3))
	set(gca,'fontsize',7);
	xlabel('LOC-aligned time (s)');
	ylabel('Frequency (Hz)');
	line([0,0],get(gca,'ylim'),'color','k','linestyle','-','linewidth',1);
	line([0,0],get(gca,'ylim'),'color','w','linestyle','--','linewidth',1);
	xlim([-400,60]);

% Compute narrowband power 
dRange = [0.5,3];
A = squeeze(nanmedian(psd(and(freq>=8,freq<=15),:,:),1)); 	A = smoothdata(A,1,'movmean',40);
B = squeeze(nanmedian(psd(and(freq>=15,freq<=30),:,:),1)); 	B = smoothdata(B,1,'movmean',40);
D = squeeze(nanmedian(psd(and(freq>=dRange(1),freq<=dRange(2)),:,:),1)); 	D = smoothdata(D,1,'movmean',40);
for i = 1:14
	ABL(i) = nanmedian(A(time<infusionTime(i),i));
	BBL(i) = nanmedian(B(time<infusionTime(i),i));
	DBL(i) = nanmedian(D(time<infusionTime(i),i));
end
AdB = 10*log(A./ABL)/log(10);
BdB = 10*log(B./BBL)/log(10);
DdB = 10*log(D./DBL)/log(10);

% Plot change in alpha, beta, and delta
axes('Position',[0.6,0.44,0.36,0.28]);
	plotwitherror(time,AdB,'CI');
	plotwitherror(time,BdB,'CI');
	plotwitherror(time,DdB,'CI');
	xlabel('LOC-aligned time (s)');
	xlim([-400,60]);
	yax = get(gca,'yaxis');
	yax.TickValues = [0,5,10,15];
	ylim([-2,15]);
	line(get(gca,'xlim'),[0,0],'color','k','linewidth',1);
	line([0,0],get(gca,'ylim'),'color','k','linewidth',1,'linestyle','--');
	gcaformat
	ylabel('Normalized power (dB)');
	clrs = lines(3);
	text(-375,14,'Delta (0.5-3 Hz)','color',clrs(3,:),'fontsize',7);
	text(-375,12.5,'Alpha (8-15 Hz)','color',clrs(1,:),'fontsize',7);
	text(-375,11,'Beta (15-30 Hz)','color',clrs(2,:),'fontsize',7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% \end{Baseline-normalized analysis} %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMax = 55;

for i = 1:13
preX(:,i) = nanmedian(psd(:,time<infusionTime(i),i),2);
postX(:,i) = nanmedian(psd(:,and(time>0,time<60),i),2);
end

[oofPre,funPre,rsquaredPre] = getFOOOF(freq(freq<=fMax),preX(freq<=fMax,:),false);
[oofPost,funPost,rsquaredPost] = getFOOOF(freq(freq<=fMax),postX(freq<=fMax,:),false);

% Compute 1/f over time
%{
m=20; N = floor((size(psd,2)-m/2)/m);
X = zeros(size(psd,1),N,size(psd,3));
it = zeros(m,N);
for i = 1:N
it(:,i) = [i*m-m/2+1:i*m+m/2];
X(:,i,:) = nanmean(psd(:,it,:),2);
end
fitTime = interp1(1:length(time),time,mean(it));
for j = 1:13
goodIdcs = find(any(~isnan(X(:,:,j)))); % Times that were not NaNs 
temp = getFOOOF(freq(freq<fMax),X(freq<fMax,goodIdcs,j),false);
alpha(goodIdcs,j) = temp(:,2);
end
%}
% Load presaved results
load('fooof_params_Cz_overtime.mat');

alpha = [];
for i = 1:14
alpha(:,i) = OOF(:,2,i);
end

idcs = [];
for i = 1:floor(max(freq)/60)
idcs(:,i) = and(freq>60*i-7,freq<60*i+5);
end
idcs = find(sum(idcs,2));
preX(idcs,:) = nan; postX(idcs,:) = nan;

axes('Position',[0.1 0.09 0.16 0.22]);
	plotwitherror(freq,log(preX)/log(10),'CI');
	plotwitherror(freq,log(postX)/log(10),'CI');
	xlim([0.5,300])
	set(gca,'xscale','log')
	xticks([1,10,100]);
	xticklabels([1,10,100]);
	xlim([0.5,100]);
	gcaformat;
	xlabel('Frequency (Hz)');
	ylabel(['Power (' char(956) 'V^2/Hz)'])
	yticks([-4:4:4]);
	yticklabels({'10^{-4}','10^{0}','10^{4}'})
	text(0.75,-2,'Baseline','color',clrs(1,:),'fontsize',7)
	text(0.75,-3,'post-LOC (0-60 s)','color',clrs(2,:),'fontsize',7)
pt = 1;
axes('Position',[0.32 0.09 0.16 0.22]);
	plot(freq,preX(:,pt),'color','k'); hold on;
	plot(freq,funPre(freq,oofPre(pt,:)),'b','LineWidth',1)
	set(gca,'xscale','log');
	set(gca,'yscale','log');
	xlim([0.5,freq(end)]);
	xlabel('Frequency (Hz)');
	xticks([0.5,5,50]);
	xticklabels([0.5,5,50]);
	xlim([0.5,fMax]);
	ylim([1e-4,1e4])
	yticks(10.^[-4:4:4])
	set(gca,'yminortick','off');
	gcaformat;
axes('Position',[0.54 0.09 0.16 0.22]);
	plot(freq,postX(:,pt),'color','k'); hold on;
	plot(freq,funPost(freq,oofPost(pt,:)),'b','LineWidth',1)
	plot(freq,funPre(freq,oofPre(pt,:)),'color',[0.6,0.6,1],'LineWidth',0.5,'linestyle','--')
	set(gca,'xscale','log');
	set(gca,'yscale','log');
	xlim([0.5,freq(end)]);
	xlabel('Frequency (Hz)');
	xticks([0.5,5,50]);
	xticklabels([0.5,5,50]);
	xlim([0.5,fMax]);
	gcaformat;
	ylim([1e-4,1e4])
	yticks(10.^[-4:4:4])
	set(gca,'yminortick','off');
axes('Position',[0.8 0.09 0.16 0.22]);
	plotwitherror(fitTime,-alpha,'CI');
	gcaformat;
	ylabel('Spectral Exponent');
	xlabel('LOC-aligned time (s)');
	xlim([-300,60]);
	ylim([-2.5,-0.5]);


labelpanel(0.02,0.93,'A');
labelpanel(0.02,0.69,'B');
labelpanel(0.52,0.69,'C');
labelpanel(0.02,0.29,'D');
labelpanel(0.25,0.29,'E');
labelpanel(0.72,0.29,'F');
