function plotfigure4

	myPath = fileparts(mfilename('fullpath'));
	addpath(fullfile(myPath,'functions'));
	addpath(fullfile(myPath,'data'));

	load('timeInformation.mat','timeInfo');
	infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;
	[F,FBL] = fittingmodel;

	load('PrePostGAFit.mat')
	M = 1;
	preFull = pre(:,M);
	postFull = post(:,M);

	load('GA_FitOverTime202105051147.mat')
	fitTime = time;

	load('spectra_pre_post.mat');
	for i = 1:13
		peaksPre(:,i) = log(preX(:,i)) - FBL(pre(:,i),freq);
		peaksPost(:,i) = log(postX(:,i)) - FBL(post(:,i),freq);
	end
	preX(and(freq>55,freq<64),:) = nan;
	postX(and(freq>55,freq<64),:) = nan;
	
	pre = zeros(15,13);
	post = pre;
	for i = 1:13
		pre(:,i) = nanmedian(params(:,find(fitTime<infusionTime(i)),i),2);
		post(:,i) = nanmedian(params(:,find(and(fitTime<60,fitTime>=0)),i),2);
	end

	delta_b = pre(7,:)./pre(11,:).*pre(2,:).*sqrt(2*pi)./exp(pre(1,:)-pre(2,:).^2/2);
	alpha_b = pre(8,:)./pre(11,:).*pre(4,:).*sqrt(2*pi)./exp(pre(3,:)-pre(4,:).^2/2);
	beta_b = pre(9,:)./pre(11,:).*pre(6,:).*sqrt(2*pi)./exp(pre(5,:)-pre(6,:).^2/2);

	AP = [pre(10,:);post(10,:)]./pre(10,:);
	gam = [pre(11,:);post(11,:)];
	tau = [pre(12,:);post(12,:)];
	tau2 = [pre(13,:);post(13,:)];
	gam = gam./gam(1,:);

	fprintf('L_AP = %.3f +/- %.3f\n',mean(AP(2,:)),stderror(AP(2,:)'))
	fprintf('L_I = %.3f +/- %.3f\n',median(gam(2,:)),stderror(gam(2,:)'))
	fprintf('TauIpre = %.3f +/- %.3f\n',mean(tau(1,:)),stderror(tau(1,:)'))
	fprintf('TauIpost = %.3f +/- %.3f\n',mean(tau(2,:)),stderror(tau(2,:)'))

	fig = figure('color','w','units','centimeters');
	fig.Position = [0,0,18.3,7.6];
	subplot(2,5,1);
		plot(freq,preX(:,M),'k','linewidth',1); hold on;
		plot(freq,exp(F(preFull(:,M),freq)),'color','b','linewidth',0.5); hold on;
		plot(freq,exp(FBL(preFull(:,M),freq)),'color','r','linewidth',1)
		set(gca,'yscale','log');
		gcaformat
		xlim([0.5,100]);
		ylim([1e-4,1e4])
		set(gca,'xscale','log');
		xticks([1,10,100]);
		set(gca,'xscale','log');
		ylabel(['PSD (' char(956) 'V^2/Hz)'])
		xlabel('Frequency (Hz)')
		yticks([1e-4,1,1e4])
		text(0.6,1e-3,'Baseline','fontsize',7)
	subplot(2,5,6);
		plot(freq,postX(:,M),'k','linewidth',1); hold on;
		plot(freq,exp(F(postFull(:,M),freq)),'color','b','linewidth',0.5); hold on;
		plot(freq,exp(FBL(postFull(:,M),freq)),'color','r','linewidth',1)
		set(gca,'yscale','log');
		gcaformat
		xlim([0.5,100]);
		ylim([1e-4,1e4])
		yticks([1e-4,1,1e4])
		set(gca,'xscale','log');
		xticks([1,10,100]);
		set(gca,'xscale','log');
		ylabel(['PSD (' char(956) 'V^2/Hz)'])
		xlabel('Frequency (Hz)')
		yticks([1e-4,1,1e4])
		text(0.6,1e-3,'post-LOC','fontsize',7)
		drawnow;
	subplot(2,5,2);
		plot(tau,'.-','color',[0.6,0.6,0.6],'linewidth',0.25,'MarkerSize',7);
		hold on;
		errorbar([0.8,2.2],mean(tau'),stderror(tau'),'color','k','linewidth',1);
		xlim([0.5,2.5]);
		xticks([1,2]); xticklabels({'BL','post-LOC'})
		gcaformat;
		ylabel('\tau_I (ms)');
	subplot(2,5,3);
		plot(gam,'.-','color',[0.6,0.6,0.6],'linewidth',0.25,'MarkerSize',7);
		hold on;
		errorbar([0.8,2.2],mean(gam'),stderror(gam'),'color','k','linewidth',1);
		xlim([0.5,2.5]);
		xticks([1,2]); xticklabels({'BL','post-LOC'})
		gcaformat;
		ylabel('\Lambda_I (norm.)');
	subplot(2,5,4);
		plot(AP,'.-','color',[0.6,0.6,0.6],'linewidth',0.25,'MarkerSize',7);
		hold on;
		errorbar([0.8,2.2],mean(AP'),stderror(AP'),'color','k','linewidth',1);
		xlim([0.5,2.5]);
		xticks([1,2]); xticklabels({'BL','post-LOC'})
		gcaformat;
		ylabel('\Lambda_{AP} (norm.)');
	subplot(2,5,7);
		h=plotmodel(1,1,1); h.Color = [0.6,0.6,0.6]; hold on;
		h=plotmodel(2.5,1,1);
		set(gca,'yscale','log'); xlim([0,100]);
		ylim([1e-4,1e4])
		yticks([1e-4,1,1e4])
		xticks([1,10,100]);
		set(gca,'xscale','log');
		xlabel('Frequency (Hz)')
		ylabel(['PSD (' char(956) 'V^2/Hz)'])
		text(1,0.8,'Baseline','color',[0.6,0.6,0.6],'fontsize',7)
		text(20,3,'2.5\tau_I','color','k','fontsize',7)
		gcaformat;
	subplot(2,5,8);
		h=plotmodel(1,1,1); h.Color = [0.6,0.6,0.6]; hold on;
		h=plotmodel(1,1.4,1);
		set(gca,'yscale','log'); xlim([0,100]);
		ylim([1e-4,1e4])
		yticks([1e-4,1,1e4])
		xticks([1,10,100]);
		set(gca,'xscale','log');
		xlabel('Frequency (Hz)')
		ylabel(['PSD (' char(956) 'V^2/Hz)'])
		text(20,5,'1.4\Lambda_I','color','k','fontsize',7)
		gcaformat;
	subplot(2,5,9);
		h=plotmodel(1,1,1); h.Color = [0.6,0.6,0.6]; hold on;
		h=plotmodel(1,1,0.01);
		set(gca,'yscale','log'); xlim([0,100]);
		ylim([1e-4,1e4])
		yticks([1e-4,1,1e4])
		xticks([1,10,100]);
		set(gca,'xscale','log');
		xlabel('Frequency (Hz)')
		ylabel(['PSD (' char(956) 'V^2/Hz)'])
		text(20,3,'0.3\Lambda_{AP}','color','k','fontsize',7)
		gcaformat;
	subplot(2,5,10);
		h=plotmodel(1,1,1); h.Color = [0.6,0.6,0.6];hold on;
		h=plotmodel(2.5,1.4,0.01);
		set(gca,'yscale','log'); xlim([0,100]);
		ylim([1e-4,1e4])
		yticks([1e-4,1,1e4])
		xticks([1,10,100]);
		set(gca,'xscale','log');
		xlabel('Frequency (Hz)')
		ylabel(['PSD (' char(956) 'V^2/Hz)'])
		gcaformat;
	subplot(2,5,5);
		plotwitherror(freq,log(preX)/log(10),'CI');
		plotwitherror(freq,log(postX)/log(10),'CI');
		gcaformat;
		xlim([0,100]);
		ylim([-4,4])
		set(gca,'xscale','log')
		yticks([-4,0,4])
		xticks([1,10,100]);
		yticklabels({'10^{-4}','10^0','10^4'})
		xlabel('Frequency (Hz)')
		ylabel(['PSD (' char(956) 'V^2/Hz)'])


	labelpanel(0.1,0.92,'A');
	% labelpanel(0.1,0.44,'B');
	labelpanel(0.265,0.92,'B');
	% labelpanel(0.26,0.44,'D');
	labelpanel(0.43,0.92,'C');
	% labelpanel(0.42,0.44,'F');
	labelpanel(0.595,0.92,'D');
	% labelpanel(0.58,0.44,'H');
	labelpanel(0.76,0.92,'E');
	% labelpanel(0.74,0.44,'J');

function h=plotmodel(tau,G,AP)

	myPath = fileparts(mfilename('fullpath'));
	dataPath = fullfile(myPath,'data');	
	load('ModelParameters.mat');
	lambda_E = lambda_E*AP;
	lambda_AP = lambda_AP*AP;
	gamma_I = gamma_I*G;
	tauI2 = tauI2*tau;

	LE = (a_0^2+a_t^2)^2*rho_N*N_E*lambda_E/(2*pi*96*pi*R_a^2*sig^2*r_0^2);
	LI = (a_0^2+a_t^2)^2*rho_N*N_I*lambda_I/(2*pi*96*pi*R_a^2*sig^2*r_0^2);
	LAP = d_AP^4*rho_N*rho_AP*lambda_AP/(2*pi*768*pi*R_a^2*sig^2*r_0^2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gam = @(t1,t2) (1-t1/t2)*(t1/t2)^(-t1/(t1-t2));
	gam2 = @(t1,t2,B) B*(1-t1/t2)*(B*t1/t2)^(-t1/(t1-t2));

	dt = 1e-6;
	t = -0.01:dt:3;


	tau_ax = Cm_ax*Rm_ax;
	lambda_ax = sqrt(d_AP*Rm_ax/(4*R_a));
	B = tauAP2/tauAP1;

	a = @(t) gamma_AP*sqrt(tau_ax)/(lambda_ax*Cm_ax*gam2(tauAP1,tauAP2,B))*[B * exp(-t/tauAP1)/sqrt(1/tauAP1-1/tau_ax) .* erf(sqrt(abs(t)*(1/tauAP1-1/tau_ax))) - exp(-t/tauAP2)/sqrt(1/tauAP2-1/tau_ax) .* erf(sqrt(abs(t)*(1/tauAP2-1/tau_ax)))].*(t>0);
	L = length(t);
	T = range(t);
	Y = fft(a(t))*dt/sqrt(T);
	pAP = 2*Y(2:ceil(L/2)+1);
	f = linspace(0,1/2,ceil(L/2)+1)/dt; f = f(2:end);

	pI = gamma_I/gam(tauI1,tauI2)*[(2*pi*sqrt(-1)*f+1/tauI1).^(-1) - (2*pi*sqrt(-1)*f+1/tauI2).^(-1)];
	pE = gamma_E/gam(tauE1,tauE2)*[(2*pi*sqrt(-1)*f+1/tauE1).^(-1) - (2*pi*sqrt(-1)*f+1/tauE2).^(-1)];



	pwE = LE.*abs(pE).^2;
	pwI = LI.*abs(pI).^2;
	pwA = LAP.*abs(pAP).^2;

	h=plot(f,pwE+pwI+pwA,'k','LineWidth',1.5); 


function [F,FBL] = fittingmodel
	gi = @(k,a1,a2) 1e3/(a1-a2)*((1e3/a1+sqrt(-1)*2*pi*k).^(-1)-(1e3/a2+sqrt(-1)*2*pi*k).^(-1));
	lgNorm = @(k,m1,s1) exp(m1-s1^2/2)*exp(-(log(k)-m1).^2/(2*s1^2))./k;
	gausses = @(m1,s1,b1,m2,s2,b2,m4,s4,b4,k) b1*lgNorm(k,m1,s1)+b2*lgNorm(k,m2,s2)+b4*lgNorm(k,m4,s4);
	gz = @(k,a1,n) 1./(1+2*pi*sqrt(-1)*(k/a1).^n);

	ftfun = @(m1,s1,m2,s2,m4,s4,b1,b2,b4,Z1,gamI,tauI,tauI2,wf,theta,f) 2* log(abs(gz(f,wf,theta))) +  ...
		log( Z1 + gamI*(1+gamI/2*gausses(m1,s1,b1,m2,s2,b2,m4,s4,b4,f)).* abs(gi(f,tauI,tauI2)).^2);
	F = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),f);



	FBL = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),0,0,0,x(10),x(11),x(12),x(13),x(14),x(15),f);
