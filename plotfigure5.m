function plotfigure5

myPath = fileparts(mfilename('fullpath'));
addpath(fullfile(myPath,'functions'));
addpath(fullfile(myPath,'data'));

load('GA_FitOverTime202105051147.mat')
fitTime = time;

load('timeInformation.mat','timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For access to spectrogram and time-domain data, 
% please contact Gilles Plourde at gilles.plourde@mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% load('psd_channel_Cz.mat','time','freq','psd')
	% for pts = 1:13;
	% 	idcs = [];
	% 	for i = 1:floor(max(freq)/60)
	% 		idcs = [idcs,find(and(freq>60*i-5,freq<60*i+5))];
	% 	end
	% 	idcs = setdiff(1:length(freq),idcs);
	% 	for i = 1:size(psd,2)
	% 		psd(:,i,pts) = interp1(freq(idcs),psd(idcs,i,pts),freq,'linear');
	% 	end
	% 	logPSD(:,:,pts) = log(psd(freq<=200,:,pts));
	% end
	freq = 0.5:0.5:200;
	logPSD = nan*zeros(400,302,13);
	time = linspace(-525.0010,228.3877,302);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gi = @(k,a1,a2) (1e3/a1+sqrt(-1)*2*pi*k).^(-1)-(1e3/a2+sqrt(-1)*2*pi*k).^(-1);
gz = @(k,a1,n) 1./(1+2*pi*sqrt(-1)*(k/a1).^n);
ftgauss1 = @(k,m1,s1) exp(-(k-m1).^2/(2*s1^2));
ftgauss2 = @(k,m1,s1) exp(m1-s1^2/2)*exp(-(log(k)-m1).^2/(2*s1^2))./k;
gausses = @(m1,s1,b1,m2,s2,b2,m4,s4,b4,k) b1*ftgauss2(k,m1,s1)+b2*ftgauss2(k,m2,s2)+b4*ftgauss2(k,m4,s4);
ftfun = @(m1,s1,m2,s2,m4,s4,b1,b2,b4,Z1,gamI,tauI,tauI2,w,n,f) 2*log(abs(Z1*gz(f,w,n)+gamI*(1+gausses(m1,s1,b1,m2,s2,b2,m4,s4,b4,f)).*gi(f,tauI,tauI2)));
F = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),f);

for k = 1:size(params,1)
	for pts = 1:13
		params2(k,:,pts) = interp1(fitTime,params(k,:,pts),time,'linear');
	end
end
params = movmedian(params2,1,2);
params = params2;

for i = 1:15
	params(i,:,1) = fillmissing(params(i,:,1),'linear');
end

M=1; % Sample patient
for pts = 1:13
	for i = 1:size(logPSD,2)
		parBestBL = params(:,i,pts);
		parBestBL(7:9) = 0;
		blPW(:,i,pts) = F(parBestBL,freq);
		peakPW(:,i,pts) = logPSD(:,i,pts) - blPW(:,i,pts);
		if(pts==M)
			P1FT(:,i) = F(params(:,i,pts),freq);
		end
	end
end


for i = 1:13
	logPSDbl(:,:,i) = repmat(nanmedian(logPSD(:,time<infusionTime(i),i),2),[1,size(logPSD,2)]);
end
normPW = 10*(logPSD -  logPSDbl)/log(10);
normPW = smoothdata(normPW,1,'movmean',2);

for i = 1:13
	peakPW2(:,:,i) = peakPW(:,:,i) - nanmedian(peakPW(:,time<infusionTime(i),i),2);
end
peakPW2 = smoothdata(peakPW,1,'movmean',2);


for i = 1:100
	P1(i,:) = fillgaps(logPSD(i,:,M));
end
P1 = imgaussfilt(P1,0.5);
P1 = logPSD(i,:,M);


dRange = [0.5,3];
A = squeeze(nanmean(normPW(and(freq>=8,freq<=15),:,:),1));
B = squeeze(nanmean(normPW(and(freq>=15,freq<=30),:,:),1));
D = squeeze(nanmean(normPW(and(freq>=dRange(1),freq<=dRange(2)),:,:),1));
A_Adj = squeeze(nanmean(peakPW2(and(freq>=8,freq<=15),:,:),1));
B_Adj = squeeze(nanmean(peakPW2(and(freq>=15,freq<=30),:,:),1));
D_Adj = squeeze(nanmean(peakPW2(and(freq>=dRange(1),freq<=dRange(2)),:,:),1));
pltParams = permute(params,[2,3,1]);
tau = squeeze(params(12,:,:)); tau = smoothdata(tau,1,'gauss',5,'includenan');
gam = squeeze(params(11,:,:)); gam = smoothdata(gam,1,'gauss',5,'includenan');
Z = squeeze(params(10,:,:)); Z = smoothdata(Z,1,'gauss',5,'includenan');
for i = 1:13
	Z(:,i) = Z(:,i)./nanmedian(Z(time(:)<infusionTime(i),i));
	gam(:,i) = gam(:,i)./nanmedian(gam(time(:)<infusionTime(i),i));
	t(:,i) = -time(:)/infusionTime(i);
end
tRescaled = linspace(-2,1,300);
for i = 1:13
	tauRescaled(:,i) = interp1(t(:,i),movmean(tau(:,i),3),tRescaled);
	gamRescaled(:,i) = interp1(t(:,i),movmean(gam(:,i),3),tRescaled);
	zRescaled(:,i) = interp1(t(:,i),Z(:,i),tRescaled);
	aRescaled(:,i) = interp1(t(:,i),movmean(A(:,i),3),tRescaled);
	bRescaled(:,i) = interp1(t(:,i),movmean(B(:,i),3),tRescaled);
	dRescaled(:,i) = interp1(t(:,i),movmean(D(:,i),3),tRescaled);
	aAdjRescaled(:,i) = interp1(t(:,i),movmean(A_Adj(:,i),3),tRescaled);
	bAdjRescaled(:,i) = interp1(t(:,i),movmean(B_Adj(:,i),3),tRescaled);
	dAdjRescaled(:,i) = interp1(t(:,i),movmean(D_Adj(:,i),3),tRescaled);
	for k = 1:size(pltParams,3)
		pltParamsRescaled(:,i,k) = interp1(t(:,i),pltParams(:,i,k),tRescaled);
	end
end


A_AdjBL = []; ABL = [];
B_AdjBL = []; BBL = [];
D_AdjBL = []; DBL = [];
A_Adjpost = []; Apost = [];
B_Adjpost = []; Bpost = [];
D_Adjpost = []; Dpost = [];


A2 = squeeze(nanmean(logPSD(and(freq>=8,freq<=15),:,:),1));
B2 = squeeze(nanmean(logPSD(and(freq>=15,freq<=30),:,:),1));
D2 = squeeze(nanmean(logPSD(and(freq>=dRange(1),freq<=dRange(2)),:,:),1));
for i = 1:size(tau,2)
	A_AdjBL(i) = log(nanmedian(exp(A_Adj(time<infusionTime(i),i))));
	B_AdjBL(i) = log(nanmedian(exp(B_Adj(time<infusionTime(i),i))));
	D_AdjBL(i) = log(nanmedian(exp(D_Adj(time<infusionTime(i),i))));

	A_Adjpre(i) = log(nanmedian(exp(A_Adj(and(time>=-30,time<-5),i))));
	B_Adjpre(i) = log(nanmedian(exp(B_Adj(and(time>=-30,time<-5),i))));
	D_Adjpre(i) = log(nanmedian(exp(D_Adj(and(time>=-30,time<-5),i))));

	A_Adjpost(i) = log(nanmedian(exp(A_Adj(and(time>=10,time<60),i))));
	B_Adjpost(i) = log(nanmedian(exp(B_Adj(and(time>=10,time<60),i))));
	D_Adjpost(i) = log(nanmedian(exp(D_Adj(and(time>=10,time<60),i))));


	ABL(i) = log(nanmedian(exp(A2(time<infusionTime(i),i))));
	BBL(i) = log(nanmedian(exp(B2(time<infusionTime(i),i))));
	DBL(i) = log(nanmedian(exp(D2(time<infusionTime(i),i))));

	Apre(i) = log(nanmedian(exp(A2(and(time>=-30,time<-5),i))));
	Bpre(i) = log(nanmedian(exp(B2(and(time>=-30,time<-5),i))));
	Dpre(i) = log(nanmedian(exp(D2(and(time>=-30,time<-5),i))));

	Apost(i) = log(nanmedian(exp(A2(and(time>=10,time<60),i))));
	Bpost(i) = log(nanmedian(exp(B2(and(time>=10,time<60),i))));
	Dpost(i) = log(nanmedian(exp(D2(and(time>=10,time<60),i))));
end



fig = figure('color','w','units','centimeters');
fig.Position = [0,0,18.3,16];

ax = axes('Position',[0.09 0.74 0.85 0.24]); axis off;
yl = ax.Position(4)/ax.Position(3);
xlim([0,1]); ylim([0,1]*yl);
rectangle('Position',[0 0 1 yl],'Curvature',0.2,'EdgeColor','none','FaceColor',[222,235,247]/255);

axes('Position',[0.15 0.81 0.14 0.14]);
	% imagesc(time,freq(1:100),P1)
	ylabel('Frequency (Hz)');
	xlabel('Drop-aligned time (s)');
	axis xy
	xlim([-360,60])
	ylim([0,40]);
	set(gca,'clim',[-5,5])
	gcaformat;
	title('Pt. #1 spectrogram')
axes('Position',[0.35 0.81 0.14 0.14]);
	imagesc(time,freq,P1FT)
	ylabel('Frequency (Hz)');
	xlabel('Drop-aligned time (s)');
	xlim([-360,60])
	axis xy;
	ylim([0,40]);
	set(gca,'clim',[-5,5])
	gcaformat;
	title('Fitted model');
	set(get(gca,'yaxis'),'visible','off');
axes('Position',[0.75 0.81 0.14 0.14]);
	% CS=P1-blPW(1:100,:,M);
	% imagesc(time,freq(1:100),CS)
	ylabel('Frequency (Hz)');
	xlabel('Drop-aligned time (s)');
	xlim([-360,60])
	axis xy;
	ylim([0,40]);
	gcaformat;
	set(gca,'clim',[2,5])
	% set(gca,'clim',[-7,3])
	title('Corrected spectrogram')
	set(get(gca,'yaxis'),'visible','off');

axes('Position',[0.55 0.81 0.14 0.14]);
	% imagesc(time,freq,blPW(:,:,M))
	ylabel('Frequency (Hz)');
	xlabel('Drop-aligned time (s)');
	xlim([-360,60])
	axis xy;
	ylim([0,40]);
	gcaformat;
	% set(gca,'clim',[-7,3])
	set(gca,'clim',[-5,5])
	title('Estimated aperiodic power')
	set(get(gca,'yaxis'),'visible','off');

axes('Position',[0.022+0.13 0.54 0.2 0.16]);
	imagesc(time,freq,nanmean(blPW,3));
	ylabel('Frequency (Hz)');
	xlabel('Drop-aligned time (s)');
	title('Avg. aperiodic power')
	xlim([-300,60])
	ylim([0,40]); axis xy;
	gcaformat;
	set(gca,'clim',[-5,5])
axes('Position',[0.022+0.41 0.54 0.2 0.16]);
	plot(time,tau,'color',[0.6,0.6,0.6]); hold on;
	plot(time,nanmedian(tau,2),'k','LineWidth',1);
	xlim([-300,60])
	ylim([0,65]);
	gcaformat;
	xlabel('Drop-aligned time (s)')
	ylabel('\tau_I (ms)');
axes('Position',[0.022+0.69 0.54 0.2 0.16]);
	plot(time,Z,'color',[0.6,0.6,0.6]); hold on;
	plot(time,nanmedian(Z,2),'k','LineWidth',1);
	xlim([-300,60])
	ylim([0,1.2]);
	gcaformat;
	xlabel('Drop-aligned time (s)')
	ylabel('\Lambda_{AP} (norm.)');
colormap jet

axes('Position',[0.022+0.13 0.3 0.2 0.16]);
	imagesc(time,freq,nanmean(logPSD,3));
	xlim([-300,60])
	ylim([0,40]); axis xy;
	gcaformat;
	set(gca,'clim',[-2.5,5])
	ylabel('Frequency (Hz)');
	xlabel('Drop-aligned time (s)');
	title('Avg. raw spectrogram')
axes('Position',[0.022+0.13 0.06 0.2 0.16]);
	imagesc(time,freq,nanmean(peakPW2,3));
	xlim([-300,60])
	ylim([0,40]); axis xy;
	gcaformat;
	set(gca,'clim',[0,7])
	xlabel('Drop-aligned time (s)');
	ylabel('Frequency (Hz)');
	title('Avg. corrected spectrogram')
colormap jet

Dtemp = smoothdata(D,1,'movmean',3);
Atemp = smoothdata(A,1,'movmean',3);
Btemp = smoothdata(B,1,'movmean',3);
axes('Position',[0.022+0.41 0.3 0.2 0.16]);
	plotwitherror(time,Atemp,true,'linewidth',0.75);
	plotwitherror(time,Btemp,true,'linewidth',0.75);
	plotwitherror(time,Dtemp,true,'linewidth',0.75);
	xlim([-360,60]); ylim([-3,15]);
	line([0,0],get(gca,'ylim'),'color','k','LineStyle','--');
	clrs = lines(3);
	text(-325,14,'Delta (0.5-3 Hz)','color',clrs(3,:),'fontsize',7);
	text(-325,12.5,'Alpha (8-15 Hz)','color',clrs(1,:),'fontsize',7);
	text(-325,11,'Beta (15-30 Hz)','color',clrs(2,:),'fontsize',7);
	gcaformat;
	xlabel('Drop-aligned time (s)');
	ylabel('Raw power (dB)');
Dtemp = smoothdata((D_Adj-D_AdjBL)*10/log(10),1,'movmean',5);
Atemp = smoothdata((A_Adj-A_AdjBL)*10/log(10),1,'movmean',5);
Btemp = smoothdata((B_Adj-B_AdjBL)*10/log(10),1,'movmean',5);
axes('Position',[0.022+0.41 0.06 0.2 0.16]);
	plotwitherror(time,Atemp,true,'linewidth',0.75);
	plotwitherror(time,Btemp,true,'linewidth',0.75);
	plotwitherror(time,Dtemp,true,'linewidth',0.75);
	xlim([-360,60]); ylim([-3,8])
	xlabel('Drop-aligned time (s)');
	line([0,0],get(gca,'ylim'),'color','k','LineStyle','--');
	% xticks(0); xticklabels({'LOC'});
	gcaformat;
	ylabel('Model-corrected power (dB)');
Dtemp = smoothdata(dRescaled,1,'movmean',1);
Atemp = smoothdata(aRescaled,1,'movmean',1);
Btemp = smoothdata(bRescaled,1,'movmean',1);
axes('Position',[0.022+0.69 0.3 0.2 0.16]);
	plotwitherror(tRescaled,Atemp,true,'linewidth',0.75);
	plotwitherror(tRescaled,Btemp,true,'linewidth',0.75);
	plotwitherror(tRescaled,Dtemp,true,'linewidth',0.75);
	xlim([-1.5,0.35]); ylim([-3,15]);
	xticks([]);
	line([0,0],get(gca,'ylim'),'color','k','LineStyle','--');
	line([-1,-1],get(gca,'ylim'),'color','k','LineStyle','--');
	gcaformat;
	xticks([-1,0]); xticklabels({'Infusion','LOC'});
	ylabel('Raw power (dB)');

Dtemp = smoothdata(dAdjRescaled,1,'movmean',1);
Atemp = smoothdata(aAdjRescaled,1,'movmean',1);
Btemp = smoothdata(bAdjRescaled,1,'movmean',1);
axes('Position',[0.022+0.69 0.06 0.2 0.16]);
	plotwitherror(tRescaled,(Atemp-A_AdjBL)*10/log(10),true,'linewidth',0.75);
	plotwitherror(tRescaled,(Btemp-B_AdjBL)*10/log(10),true,'linewidth',0.75);
	plotwitherror(tRescaled,(Dtemp-D_AdjBL)*10/log(10),true,'linewidth',0.75);
	xlim([-1.5,0.35]); ylim([-3,8])
	line([0,0],get(gca,'ylim'),'color','k','LineStyle','--');
	line([-1,-1],get(gca,'ylim'),'color','k','LineStyle','--');
	xlabel('Rescaled time');
	xticks([-1,0]); xticklabels({'Infusion','LOC'});
	gcaformat;
	ylabel('Model-corrected power (dB)');



annotation('arrow','headstyle','plain','linewidth',4,'position',[0.3 0.88 0.035 0]);
annotation('arrow','headstyle','plain','linewidth',4,'position',[0.5 0.88 0.035 0]);
annotation('arrow','headstyle','plain','linewidth',4,'position',[0.7 0.88 0.035 0]);


labelpanel(0.1,0.92,'A');
labelpanel(0.1,0.67,'B');
labelpanel(0.38,0.67,'C');
labelpanel(0.65,0.67,'D');
labelpanel(0.1,0.44,'E');
labelpanel(0.38,0.44,'F');
labelpanel(0.65,0.44,'G');
labelpanel(0.1,0.19,'H');
labelpanel(0.38,0.19,'I');
labelpanel(0.65,0.19,'J');