% function plotfigure2

myPath = fileparts(mfilename('fullpath'));
addpath(myPath,'functions');
addpath(myPath,'data');

	
load('timeInformation.mat','timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For access to spectrogram and time-domain data, 
% please contact Gilles Plourde at gilles.plourde@mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% load('psd_channel_Cz.mat','freq','time','psd');
	freq = 0.5:0.5:200;
	psd = nan*zeros(400,7491,14);
	time = linspace(-525.0010,228.3877,7491);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=10; N = size(psd,2)-m;
X = zeros(size(psd,1),N,size(psd,3));
it = zeros(m,N);
for i = 1:N
	it(:,i) = [i+1:i+m];
	X(:,i,:) = nanmean(psd(:,it(:,i),:),2);
end
fitTime = interp1(1:length(time),time,mean(it));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load pre-saved results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fMax = 55;
	% OOF = nan*zeros(size(X,2),2,size(X,3));
	% for j = 1:13
		% goodIdcs = find(any(~isnan(X(:,:,j)))); % Times that were not NaNs 
		% [OOF(goodIdcs,:,j),fun] = getFOOOF(freq(freq<fMax),X(freq<fMax,goodIdcs,j),false);
	% end
	load('fooof_params_Cz_overtime.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(OOF,3)
	for j = 1:size(OOF,1)
		oof(:,j,i) = fun(freq(freq<fMax),OOF(j,:,i));
	end
end
detrendedX = 10*log(X(freq<fMax,:,:)./oof)/log(10);


load('spectra_pre_post.mat');
[oofPre,funPre,rsquaredPre] = getFOOOF(freq(freq<=fMax),preX(freq<=fMax,:),false);
[oofPost,funPost,rsquaredPost] = getFOOOF(freq(freq<=fMax),postX(freq<=fMax,:),false);
dRange = [1,4];
for i = 1:13
	detrendPreX(:,i) = 10*log(preX(:,i)./funPre(freq,oofPre(i,:)))/log(10);
	detrendPostX(:,i) = 10*log(postX(:,i)./funPost(freq,oofPost(i,:)))/log(10);
	preD(i) = nanmean(detrendPreX(and(freq>=dRange(1),freq<dRange(2)),i));
	postD(i) = nanmean(detrendPostX(and(freq>=dRange(1),freq<dRange(2)),i));
end



fig = figure('color','w','units','centimeters');
% fig.Position(3) = 13.5;
fig.Position(3) = 18.3;
fig.Position(4) = 8.9;
axes('Position',[0.06,0.6,0.19,0.33]);
	imagesc(time(:),freq,log(nanmean(psd,3)));
	ylim([0.5,40]);
	axis xy;
	set(gca,'clim',[-4,7])
	xlim([-300,60]);
	xlabel('Drop-Aligned Time (s)')
	ylabel('Frequency (Hz)')
	TTL = title('Group-level spectrogram','FontWeight','bold','FontSize',7,...
				'HorizontalAlignment','center'); 
	TTL.Position(2) = 41;

	set(gca,'FontSize',7);
axes('Position',[0.315,0.6,0.19,0.33]);
	imagesc(fitTime,freq(freq<fMax),log(nanmedian(oof,3)));
	ylim([0.5,40])
	% set(gca,'CLim',[-5,1])
	xlim([-300,60]);
	xlabel('Drop-Aligned Time (s)')
	ylabel('Frequency (Hz)')
	TTL = title('1/f fit','FontWeight','bold','FontSize',7,...
				'HorizontalAlignment','center'); 
	TTL.Position(2) = 41;
	set(gca,'FontSize',7);
	axis xy;

axes('Position',[0.57,0.6,0.19,0.33]);
	imagesc(fitTime,freq(freq<fMax),nanmedian(detrendedX,3));
	set(gca,'CLim',[3,12]);
	ylim([0.5,40]);
	axis xy;
	xlim([-300,60]);
	xlabel('Drop-Aligned Time (s)')
	ylabel('Frequency (Hz)')
	TTL = title('1/f corrected','FontWeight','bold','FontSize',7,...
				'HorizontalAlignment','center'); 
	TTL.Position(2) = 41;
	set(gca,'FontSize',7);
axes('Position',[0.82,0.6,0.15,0.33]);
	plot([preD(:),postD(:)]','.-','color',[0.6,0.6,0.6],'markerSize',5); hold on;
	errorbar([0.9,2.1],mean([preD(:),postD(:)]),stderror([preD(:),postD(:)]),'color','k','LineWidth',1);
	xticks([1,2]);
	xticklabels({'BL','post-LOC'});
	xlim([0.5,2.5]);
	gcaformat;
	ylabel('Delta power (dB)')
	xlabel(sprintf('n.s. (p = %0.2f)',ranksum(preD(:),postD(:))));
	ylim([-1,5])

colormap jet

example = load('sampleTimeSeries1.mat');
idcs = find(and(time>35,time<60));
idcs2 = find(and(example.time>35,example.time<60));
T = example.time(idcs2);
X = example.timedomain(idcs2);

fMax = 55;
[~,~,psd1] = eegfft(T,X,2,1.9);
psd1 = nanmedian(psd1,2);
[OOF,oofF,pks] = getFOOOF(freq(freq<fMax), ...
	psd1(freq<fMax),false);

axes('Position',[0.06,0.12,0.15,0.33]);
	plot(freq,log(psd1)/log(10),'k','LineWidth',1);
	xlim([0.5,fMax]);
	ylabel('log power');
	xlim([0.5,55]);
	xticks([0.5,5,50])
	xlabel('Frequency (Hz)');
	set(gca,'xscale','log')
	ylabel(['Power (' char(956) 'V^2/Hz)'])
	ylim([-2,4])
	yticks([-2:2:4]);
	yticklabels({'10^{-2}','10^{0}','10^{2}','10^{4}'})
	hold on;
	plot(freq,log(oofF(freq,OOF))/log(10),'color','b','linestyle','--');
	gcaformat;

	for i = 1:length(pks.freq)
		scatter(pks.freq(i),1.25/log(10)+interp1(freq,log(psd1)/log(10),pks.freq(i)),20,[1,0,0],'filled','v')
		tx(i) = text(pks.freq(i),1.5/log(10)+interp1(freq,log(psd1)/log(10),pks.freq(i)),[num2str(pks.freq(i),2) ' Hz'], ...
				'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',7);
	end
	tx(1).Position = [2.6,1.5];
	tx(2).Position = [10.5,1.2];
	tx(3).Position = [18.9,0];


axes('Position',[0.28,0.11,0.465,0.34]);
	plot(T,X,'k','LineWidth',0.2);
	hold on;

	[b,a] = butter(4,5/512);
	plot(T,filtfilt(b,a,X),'r','linewidth',1);
	plot(T,smoothdata(X),'r','LineWidth',1);
	xlim([T(1),T(end)]);
	axis off;
	gcaformat;
	line([35,37],[-75,-75],'color','k','linewidth',1)
	text(36,-80,'2 s','VerticalAlignment','top','HorizontalAlignment','center','fontsize',7);

	[y,t] =	findpeaks(smoothdata(X),T,'MinPeakDistance',0.25,'MinPeakProminence',2);
	hold on;
	scatter(t,interp1(T,movmax(X,500),t)+10,10,'*r');

axes('Position',[0.82,0.12,0.15,0.33]);
	BE = exp(linspace(log(0.25),log(4),20));
	h=histogram(1./diff(t),'BinEdges',BE);
	h.EdgeColor = [1,0.5,0.5]; h.LineWidth = 0.75;
	h.FaceColor = 'r';
	xlim([0.5,4])
	set(gca,'xscale','log');
	xlabel('1/inter-peak interval (Hz)')
	ylabel('Count')
	xax = get(gca,'xaxis');
	xax.MinorTickValues = [0.5:0.5:4];
	gcaformat





labelpanel(0.01,0.93,'A');
labelpanel(0.265,0.93,'B');
labelpanel(0.52,0.93,'C');
labelpanel(0.77,0.93,'D');
labelpanel(0.01,0.43,'E');
labelpanel(0.28,0.43,'F');
labelpanel(0.77,0.43,'G');
