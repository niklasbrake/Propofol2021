function plotfigure2S1

myPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(myPath,'data');
	
load(fullfile(dataPath,'timeInformation.mat'),'timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For access to time-domain data, please contact 
% Gilles Plourde at gilles.plourde@mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('FlaviePtsDropAlignedTimeDomain20191012.mat');
error('For access to time-domain data, please contact Gilles Plourde at gilles.plourde@mcgill.ca');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fMax = 55;

figure('color','w');
% P = [1,3,4,5,6,9,10,12,13];
P = 1:14;
for ipts = 1:length(P);
	pt = P(ipts);
	ax1 = axes('Position',[0.01+mod((15-ipts)-1,3)/3,0.02+0.2*floor((15-ipts)/3),0.24,0.18]);
	idcs = find(and(Time>35,Time<60));
	if(pt==9)
		idcs = find(and(Time>15,Time<40));
	end
	if(pt==11)
		idcs = find(and(Time>45,Time<60));
	end
	if(pt == 14)
		idcs = find(and(Time>22,Time<37));
	end
	X = TimeDomainAligned(idcs,2,pt);
	I = find(~isnan(X));
	T = linspace(0,25,length(idcs));
	plot(T,X-170,'r','LineWidth',0.2);
	hold on;
	[b,a] = butter(4,5/512);
	plot(T,filtfilt(b,a,X)-170,'k','linewidth',0.75);
	axis off;
	gcaformat;

	T = Time(idcs);
	temp = interp1(T(I),movmean(X(I),100),T,'linear');
	X(isnan(X)) = temp(isnan(X));
	[freq,~,psd] = eegfft(Time(idcs),X,2,1.9);
	postX = nanmean(psd,2);
	

	% PRE
	idcs = find(Time<infusionTime(pt));
	X = TimeDomainAligned(idcs,2,pt);
	n = ~isnan(X)';
	f = find(diff([false,n==1,false])~=0);
	[m,ix] = max(f(2:2:end)-f(1:2:end-1));

	idcs = find(and(Time>Time(f(2*ix-1)),Time<Time(f(2*ix-1))+25));
	if(ipts==6)
		idcs = idcs+5000;
	end

	T = linspace(0,25,length(idcs));
	X = TimeDomainAligned(idcs,2,pt);
	plot(T,X,'b','LineWidth',0.2);
	hold on;

	[b,a] = butter(4,5/512);
	plot(T,filtfilt(b,a,X),'k','linewidth',0.75);
	axis off;
	gcaformat;
	xlim([0,25]);


	T = Time(idcs);
	temp = interp1(T(I),movmean(X(I),100),T,'linear');
	X(isnan(X)) = temp(isnan(X));
	[freq,~,psd] = eegfft(Time(idcs),X,2,1.9);
	preX = nanmean(psd,2);


	[oofPre,funPre,rsquaredPre] = getFOOOF(freq(freq<=fMax),preX(freq<=fMax,:),false);
	[oofPost,funPost,rsquaredPost] = getFOOOF(freq(freq<=fMax),postX(freq<=fMax,:),false);


	ax2 = axes('Position',[0.265+mod((15-ipts)-1,3)/3,0.04+0.2*floor((15-ipts)/3),0.06,0.16]);
	plot(freq,postX./funPost(freq,oofPost),'LineWidth',0.75,'color','r'); hold on;
	plot(freq,preX./funPre(freq,oofPre),'LineWidth',0.75,'color','b');
	xlim([0,55]);
	set(gca,'xscale','log');
	gcaformat
	yticks([]);
	xticks([1,10]);
	xticklabels({});
	if(sum(ipts == 12:14))
		xticklabels({'1','10 Hz'});
	end
	drawnow;
end


