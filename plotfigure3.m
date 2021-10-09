function plotfigure3

myPath = fileparts(mfilename('fullpath'));
addpath(myPath,'functions');
addpath(myPath,'data');

load('ModelParameters.mat');


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

fig = figure('color','w','units','centimeters');
fig.Position = [0,0,9,4];
ax = subplot(1,2,1);
	h(1)=plot(f,pwE,'-r','LineWidth',1); hold on;
	h(2)=plot(f,pwI,'-b','LineWidth',1);
	h(3)=plot(f,pwA,'-g','LineWidth',1);
	plot(f,pwE+pwI+pwA,'k','LineWidth',1.5); 
	set(gca,'yscale','log')
	set(gca,'xscale','log')
	xlim([0.5,400]);
	xticks([0.1,1,10,100,1e3]);
	xticklabels([0.1,1,10,100,1e3]);
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
	ax.GridColor = [0,0,0]; ax.GridAlpha = 0.2;
	ax.XAxis.Color = [0.5,0.5,0.5]; ax.XAxis.Label.Color = 'k';
	ax.YAxis.Color = [0.5,0.5,0.5]; ax.YAxis.Label.Color = 'k';
	box on; set(gca,'tickdir','in');


load('PrePostGAFit.mat')
load('timeInformation.mat','timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;

load('psd_channel_Cz.mat','freq','time','psd');
for i = 1:13
	preX(:,i) = nanmedian(psd(:,time<infusionTime(i),i),2);
	postX(:,i) = nanmedian(psd(:,and(time>0,time<60),i),2);
end

[F,FBL] = fittingmodel2;


pts = 12;
c = pre(:,pts);
c(14) = Inf;
cBL = c; cBL(7:9) = 0;
c1 = cBL; c1(10) = 0;
c2 = cBL; c2(11) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = subplot(1,2,2);
	plot(freq,preX(:,pts),'k','LineWidth',1); hold on;
	plot(freq,exp(FBL(c,freq)),'r'); 
	plot(freq,exp(F(c1,freq)),'b','LineWidth',1);
	plot(freq,exp(F(c2,freq)),'g','LineWidth',1);
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
