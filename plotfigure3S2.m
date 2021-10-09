function plotfigure3S2

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

b = 0.01;
a= 0.01;
pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*g(1,10,f)).*abs(pE).^2;
pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*g(1,10,f)).*abs(pI).^2;
pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*g(1,10,f)).*abs(pAP).^2;

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
	title('a = 0.01, c = 0.01');

b = 0.01;
a= 0.1;
pwE = sc1*(LE*(1-b+b*a)+(b*a*LE)^2*g(1,10,f)).*abs(pE).^2;
pwI = sc1*(LI*(1-b+b*a)+(b*a*LI)^2*g(1,10,f)).*abs(pI).^2;
pwA = sc2*(LAP*(1-b+b*a)+(b*a*LAP)^2*g(1,10,f)).*abs(pAP).^2;
ax = subplot(1,2,2);
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
	title('a = 0.1, c = 0.01');