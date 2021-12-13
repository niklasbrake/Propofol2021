% function plotfigure3S2

myPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(myPath,'data');	
load(fullfile(dataPath,'ModelParameters.mat'));

sc1 = (a_0^2+a_t^2)^2/(192*pi^2*R_a^2*sig^2*r_0^2);
sc2 = d_AP^4/(1536*pi^2*R_a^2*sig^2*r_0^2);

LE = rho_N*N_E*lambda_E;
LI = rho_N*N_I*lambda_I;
LAP = rho_N*rho_AP*lambda_AP;

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


g = @(s,m,x) 1/sqrt(2*pi*s^2) * exp(-(x-m).^2/(2*s^2));
lgNorm = @(k,m1,s1) exp(-(log(k)-m1).^2/(2*s1^2))./(k*s1*sqrt(2*pi));

fig = figure('color','w','units','centimeters');
b=0;
a=1e-4;
tauI2 = 20e-3;
gamma_I = 750;
pI = gamma_I/gam(tauI1,tauI2)*[(2*pi*sqrt(-1)*f+1/tauI1).^(-1) - (2*pi*sqrt(-1)*f+1/tauI2).^(-1)];
pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*lgNorm(f,0.5,0.4)).*abs(pE).^2;
pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*lgNorm(f,0.5,0.4)).*abs(pI).^2;
pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*lgNorm(f,0.5,0.4)).*abs(pAP).^2;
Y = pwE+pwI+pwA;
plot(f,Y,'k','LineWidth',1); hold on;

b=0.6;
a=1e-4;
tauI2 = 60e-3;
gamma_I = gamma_I*2;
pI = gamma_I/gam(tauI1,tauI2)*[(2*pi*sqrt(-1)*f+1/tauI1).^(-1) - (2*pi*sqrt(-1)*f+1/tauI2).^(-1)];
pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*lgNorm(f,0.5,0.4)).*abs(pE).^2;
pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*lgNorm(f,0.5,0.4)).*abs(pI).^2;
pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*lgNorm(f,0.5,0.4)).*abs(pAP).^2;
Y = pwE+pwI+pwA;
plot(f,Y,'r','LineWidth',1); hold on;
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([0.5,500]);

freq = linspace(0,400,1e3);
pEK = interp1(f,abs(pE).^2,freq);
pIK = interp1(f,abs(pI).^2,freq);
pAPK = interp1(f,abs(pE).^2,freq);

B = flip(10.^linspace(-3,-0.1,20));

a=1e-3;
for i = 1:length(B)
	b=B(i);

	pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*g(2,15,freq)).*pEK;
	pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*g(2,15,freq)).*pIK;
	pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*g(2,15,freq)).*pAPK;
	Ypre = pwE+pwI+pwA;

	r(:,i) = Ypre;
end
figure;
CM = parula(length(B));
for i=  1:length(B)
	plot(freq,r(:,i),'color',CM(i,:)); hold on;
end
set(gca,'xscale','log');
set(gca,'yscale','log');


B = flip(10.^linspace(-4,0,20));
for i = 1:length(B)
	b=B(i);
	a = 1e-4/b;

	pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*g(2,15,f)).*abs(pE).^2;
	pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*g(2,15,f)).*abs(pI).^2;
	pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*g(2,15,f)).*abs(pAP).^2;
	Ypre = pwE+pwI+pwA;

	pwE = sc1*(LE*(1-b+b*a)+0*(b*a*LE)^2*g(2,15,f)).*abs(pE).^2;
	pwI = sc1*(LI*(1-b+b*a)+0*(b*a*LI)^2*g(2,15,f)).*abs(pI).^2;
	pwA = sc2*(LAP*(1-b+b*a)+0*(b*a*LAP)^2*g(2,15,f)).*abs(pAP).^2;
	Ypost = pwE+pwI+pwA;

	[m,I] = max(Ypre);
	R(i) = log(m/Ypost(I))/log(10);
end
plot(B,R,'.-')
set(gca,'xscale','log')

for i = 1:5
	for j = 1:5
		b = 10^(1-i);
		a = 10^(1-j);
		pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*g(2,15,f)).*abs(pE).^2;
		pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*g(2,15,f)).*abs(pI).^2;
		pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*g(2,15,f)).*abs(pAP).^2;
		Y = pwE+pwI+pwA;
		ax = subplot(5,5,5*(i-1)+j);
		plot(f,Y,'k','LineWidth',1.5); hold on;
		[y,x] = findpeaks(Y,f,'NPeaks',1);
		scatter(x,y*4,10,[1,0,0],'v','filled');


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
		title(['a = ' num2str(a) ', c = ' num2str(b)]);
	end
end


load('psd_channel_cz');
cla;
hold off;
Y = nanmean(psd(:,5e3:5.1e3,5),2);
plot(freq,Y,'linewidth',1); hold on;
[y,x] = findpeaks(movmean(Y,10),freq,'NPeaks',1,'MinPeakProminence',1);
scatter(x,y*4,10,[1,0,0],'v','filled');
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
title('Data sample')