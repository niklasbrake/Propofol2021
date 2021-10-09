function slope = plotfigureS4(tau)

Fs = 15e3;
tmax = 20;
t = [0:1/Fs:0.1];

gam = @(t1,t2) (1-t1/t2)*(t1/t2)^(-t1/(t1-t2));

ecov = exp(-1e3*t/2)-exp(-1e3*t/0.1);
mE = gam(0.1,2);
ecov = ecov/mE; sE = sum(ecov);
icov = exp(-1e3*t/tau)-exp(-1e3*t/0.5);
mI = gam(0.5,tau);
icov = icov/mI; sI = sum(icov);

EIR = 4*(2*8e3/Fs*sE)/(5*2000/Fs*sI); % ratio of expected values of E and I

L = tmax*Fs;
hmw = hamming(L);

f = Fs * linspace(0,1/2,L/2+1);
f = f(2:end);

N=40;
pEx = zeros(N,sum(f<150));
pIn = zeros(N,sum(f<150));

for j = 1:N
	E = poissrnd(2/Fs*8e3,floor(tmax*Fs+length(t)-1),1);
	I = poissrnd(5/Fs*2e3,floor(tmax*Fs+length(t)-1),1);

	gE = conv(E,ecov,'valid');
	gI = conv(I,icov,'valid');

	gI = EIR*gI;

	IPSC = gI*15;
	EPSC = gE*65;

	LFP = IPSC + EPSC;

	temp = EPSC;
	L = length(temp);
	temp = temp.*hmw;
	P = sqrt(2)*fft(temp)/sqrt(L*Fs);
	P = P(1:L/2+1);
	pEx(j,:) = P(f<150);

	temp = IPSC;
	L = length(temp);
	temp = temp.*hmw;
	P = sqrt(2)*fft(temp)/sqrt(L*Fs);
	P = P(1:L/2+1);
	pIn(j,:) = P(f<150);

	temp = LFP;
	L = length(temp);
	temp = temp.*hmw;
	P = sqrt(2)*fft(temp)/sqrt(L*Fs);
	P = P(1:L/2+1);
	pLFP(j,:) = P(f<150);
end

F = exp(linspace(log(1),log(150),1e3));

PE = interp1(f(f<150),smoothdata(abs(pEx').^2,'movmean',5),F);
PI = interp1(f(f<150),smoothdata(abs(pIn').^2,'movmean',5),F);
PL = interp1(f(f<150),smoothdata(abs(pLFP').^2,'movmean',5),F);

ft = @(k,a1,a2) (1e3/a1+sqrt(-1)*2*pi*k).^(-1) - (1e3/a2+sqrt(-1)*2*pi*k).^(-1);
A = 65*sqrt(8e3*2/mE)*ft(F,0.1,2);%.*exp(sqrt(-1)*2*pi*rand(size(F)));
B = EIR*15*sqrt(2e3*5/mI)*ft(F,0.5,tau);%.*exp(sqrt(-1)*2*pi*rand(size(F)));

	it = find(and(F>30,F<=50));
	for i = 1:N
		FT = fitlm(log(F(it))/log(10),log(PL(it,i))/log(10));
		slope(i) = FT.Coefficients{2,1};
	end

	% return;
	
fig = figure('color','w','units','centimeters');
fig.Position(3:4) = [8.7,7];
	plotwitherror(F,PE,'CI');
	plotwitherror(F,PI,'CI');
	plotwitherror(F,PL,'CI');
	plot(F,abs(A).^2+abs(B).^2,'-k','LineWidth',1);
	plot(F,abs(A).^2,'-k','LineWidth',1);
	plot(F,abs(B).^2,'-k','LineWidth',1);
	set(gca,'yscale','log')
	set(gca,'xscale','log')
	xlim([1,150]);
	% ylim([6,1.3e3]);
	xlabel('Frequency (Hz)')
	ylabel('log power')
	xticks([1,10,100]);
	xticklabels({'1','10','100'});
	box off;
	set(gca,'fontsize',9); set(gca,'LineWidth',0.75);



	title(sprintf('slope = %.3f',mean(slope)));