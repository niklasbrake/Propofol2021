% function plotfigure3S1

myPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(myPath,'data');	
load(fullfile(dataPath,'ModelParameters.mat'));

tau_ax= Cm_ax*Rm_ax*1e3;
lambda_ax = sqrt(d_AP*Rm_ax/(4*R_a));
B = tauAP2/tauAP1;

gam2 = @(t1,t2,B) B*(1-t1/t2)*(B*t1/t2)^(-t1/(t1-t2));

dt = 1e-6;
t2 = -0.01:dt:3;
aAP = @(t) gamma_AP*sqrt(tau_ax)/(lambda_ax*Cm_ax*gam2(tauAP1,tauAP2,B))*[B * exp(-t/tauAP1)/sqrt(1/tauAP1-1/tau_ax) .* erf(sqrt(abs(t)*(1/tauAP1-1/tau_ax))) - exp(-t/tauAP2)/sqrt(1/tauAP2-1/tau_ax) .* erf(sqrt(abs(t)*(1/tauAP2-1/tau_ax)))].*(t>0);
L = length(t2);
T = range(t2);
Y = fft(aAP(t2))*dt/sqrt(T);
pAP = 2*Y(2:ceil(L/2)+1);
f = linspace(0,1/2,ceil(L/2)+1)/dt; f = f(2:end);

dx = 1e-3;
dt = 1e-3;
x = -0.2:dx:0.2;
t = 0:dt:4;

V = zeros(length(x),length(t));
IAP = @(t,x) gamma_AP/gam2((tauAP1*1e3),(tauAP2*1e3),B)*(B*exp(-(t-1)/(tauAP1*1e3))-exp(-(t-1)/(tauAP2*1e3))).*(t>=1).*(abs(x)<=1.1e-3);

for j = 1:length(t)
	for i = 1:length(x)
		s = t(1:j)';
		V(i,j) = -60e3 + Rm_ax/tau_ax * nansum(dt*exp(-s/tau_ax)./sqrt(4*pi*lambda_ax^2*s/tau_ax).*IAP(t(j)-s,0).*exp(-x(i)^2./(4*lambda_ax^2*s/tau_ax)));
	end
end

fig = figure('color','w','units','centimeters');
fig.Position = [0,0,12,8.5];
subplot(2,2,1);
	plot(t-1,IAP(t,0),'k','LineWidth',1);
	gcaformat;
	ylabel(['I_{AP} (' char(956) 'A/mm)']);
	xlabel('Time (ms)')
	xlim([-1,3]);
subplot(2,2,2);
	imagesc(t-1,x,V*1e-3);
	gcaformat;
	box on; set(gca,'tickdir','in');
	xlabel('Time (ms)');
	ylabel('Distance (mm)');
	colormap(parula(255));
	C = colorbar('Position',[0.9182 0.3088 0.0266 0.4854]);
	C.Limits = [-9e4,6e4]*1e-3;
	set(gca,'CLim',[-90,60]);
	C.Label.String = ['Membrane Potential (mV)'];
subplot(2,2,3);
	plot(t2*1e3,aAP(t2),'k','LineWidth',1);
	gcaformat;
	xlabel('Time (ms)')
	ylabel('\alpha_{AP}')
	xlim([-1,3])
	yticks([]);
subplot(2,2,4);
	plot(f,10*log(pAP)/log(10),'k','LineWidth',1);
	gcaformat;
	xlabel('Frequency (Hz)')
	xticks([1,10,100,1e3]);
	yticks([]);
	xlim([1,1e3]);
	set(gca,'xscale','log');
	ylabel('log power')


	X = V(:,1.5e3);

	X1d1 = diff(X);
	% X2d2=X(3:end)+X(1:end-2)-X(2:end-1)*2;