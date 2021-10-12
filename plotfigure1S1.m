function plotfigure1

myPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(myPath,'data');
	
LUT = [{'Fz'}, {'Pz'}, {'C3'}, {'C4'}, {'CP3'}, {'CP4'}];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For access to time-domain data, please contact 
% Gilles Plourde at gilles.plourde@mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('FlaviePtsDropAlignedTimeDomain20191012.mat');
error('For access to time-domain data, please contact Gilles Plourde at gilles.plourde@mcgill.ca');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% \begin{Baseline-normalized analysis} %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(dataPath,'timeInformation.mat'),'timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;
	

fig = figure('color','w','units','centimeters');
for k = 1:6
	load(fullfile(dataPath,['psd_channel_' LUT{k}]));

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
	subplot(2,3,k);
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

		if(k==1)
			clrs = lines(3);
			text(-375,14,'Delta (0.5-3 Hz)','color',clrs(3,:),'fontsize',7);
			text(-375,12.5,'Alpha (8-15 Hz)','color',clrs(1,:),'fontsize',7);
			text(-375,11,'Beta (15-30 Hz)','color',clrs(2,:),'fontsize',7);
			drawnow;
		end

		title(['Channel ' LUT{k}])
end