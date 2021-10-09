function c0 = fitGAsingle(freq,psd,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myPath = fileparts(mfilename('fullpath'));
dataPath = fullfile(fileparts(myPath),'data');	
load(fullfile(dataPath,'timeInformation.mat'),'timeInfo');
infusionTime = timeInfo.infusion_onset-timeInfo.object_drop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ip = inputParser;

paramName = 'StartPoint';
	defaultval = [0.1,0.4,2.5464,0.3037,3.1887,0.1172,...
				4,4,4, ...
				0.0679,2.2464, ...
				15.2168,3.1644,382.8012,2.4416];
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'N';
	defaultval = 200;
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'M';
	defaultval = 30;
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);


paramName = 'fig';
	defaultval = [];
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) true;
	addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'coolingOnset';
	defaultval = 100;
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

parse(ip,varargin{:});
	x0 = ip.Results.StartPoint;
	N = ip.Results.N;
	M = ip.Results.M;
	fig = ip.Results.fig;
	coolingOnset = ip.Results.coolingOnset;
	if(isempty(fig))
		figure;
		fig = gca;
	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F,FBL] = fittingmodel;

[xData, yData] = prepareCurveData( freq, log(psd) );
yData(xData>250) = [];
xData(xData>250) = [];
w = ones(size(xData));
for i = 1:floor(max(xData)/60)
	idcs = find(and(xData>60*i-5,xData<60*i+5));
	w(idcs) = 0;
end
idcs = find(w==1);
yData = interp1(xData(idcs),yData(idcs),xData,'linear');
f = @(x) sum(abs(yData-F(x,xData))./xData) + sum((yData<FBL(x,xData))./xData);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Proposal distributions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mu of delta peak
	a = -1; b = log(7); sig = (b-a)/50;
	A1 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{1} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A1(x0,T),x0,T*sig);
% Sigma of delta peak
	a = 0.01; b = 1.5; sig = (b-a)/50;
	A2 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{2} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A2(x0,T),x0,T*sig);
% Mu of alpha peak
	a = log(8); b = log(15); sig = (b-a)/10;
	A3 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{3} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A3(x0,T),x0,T*sig);
% Sigma of alpha peak
	a = 0.01; b = 2; sig = (b-a)/50;
	A4 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{4} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A4(x0,T),x0,T*sig);
% Mu of beta peak
	a = log(15); b = log(30); sig = (b-a)/10;
	A5 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{5} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A5(x0,T),x0,T*sig);
% Sigma of beta peak
	a = 0.01; b = 2; sig = (b-a)/50;
	A6 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{6} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A6(x0,T),x0,T*sig);
% Height of delta peak
	a = 0; b = 30; sig = (b-a)/50;
	A7 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{7} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A7(x0,T),x0,T*sig);
% Height of alpha peak
	a = 0; b = 30; sig = (b-a)/50;
	A8 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{8} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A8(x0,T),x0,T*sig);
% Height of beta peak
	a = 0; b = 30; sig = (b-a)/50;
	A9 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{9} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A9(x0,T),x0,T*sig);
% AP firing rate
	a = 0; b = 1; sig =(b-a)/50;
	A10 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{10} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A10(x0,T),x0,T*sig);
% Scaling factor of IPSPs
	a = 0; b = 50; sig = (b-a)/50;
	A11 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{11} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A11(x0,T),x0,T*sig);
% Decay time of IPSPs
	a = 10; b = 80; sig = (b-a)/50;
	A12 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{12} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A12(x0,T),x0,T*sig);
% Rise time of IPSPs
	a = 0.5; b = 7; sig = (b-a)/50;
	A13 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{13} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A13(x0,T),x0,T*sig);
% Z2
	a = 200; b = 500; sig = 5;
	A14 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{14} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A14(x0,T),x0,T*sig);
% Z2
	a = 0; b = 10; sig = 0.1;
	A15 = @(x0,T) cdf('norm',b,x0,T*sig)-cdf('norm',a,x0,T*sig);
	g{15} = @(x0,T) icdf('norm',cdf('norm',a,x0,T*sig)+rand*A15(x0,T),x0,T*sig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gene{1} = [1,2,7];
gene{2} = [3,4,8];
gene{3} = [5,6,9];
gene{4} = 10;
gene{5} = 11;
gene{6} = [12,13];
gene{7} = [14,15];
nGenes = length(gene);
nPars = length(g);

x0 = repmat(x0(:)',[M,1]);
x1 = x0;
L1 = zeros(M,1);
gof = zeros(N,1);
best = zeros(nPars,N);
fValues = linspace(0.5,max(xData),1e3);

% axes(fig);
% 	cla; hold off;
% 	plot(xData,yData,'.-');
% 	hold on;
% 	hPlot1 = plot(fValues,F(x0(1,:),fValues));
% 	c0BL = x0(1,:);
% 	c0BL(7:9) = 0;
% 	hPlot2 = plot(fValues,F(c0BL,fValues));
% 	xlim([min(xData),max(xData)]);
% 	set(gca,'xscale','log');
% 	gcaformat;
% 	xlabel('Frequency (Hz)');
% 	yl = get(gca,'ylim');
% 	txt = text(1,yl(1)+0.1*range(yl),['i = ' int2str(i)]);


for i = 1:N

	for j = 1:M
		Ltemp = f(x0(j,:)) * (1-0.9*exp(-i/10));
		osc = nansum(abs(FBL(x0(j,:),xData)-F(x0(j,:),xData)));
		L1(j) = max(exp(-(Ltemp))*(1+(i>100)*0.5/(1+exp(osc/100))),1e-11);
	end
	[gof(i),I] = max(L1);
	best(:,i) = x0(I,:);

	if(i==N)
		break;
	end

	LS = cumsum(L1)/sum(L1);
	iParents = interp1(LS,1:length(LS),rand(M,2),'next','extrap');
	iGene = randi(nGenes,M,1);

	T = ones(nPars,1);
	if(i>coolingOnset)
		T = T*0.97;
	end
	for j = 1:M
		iGenes1 = randperm(nGenes,iGene(j));
		iGenes2 = setdiff(1:nGenes,iGenes1);
		idcs1 = [gene{iGenes1}];
		idcs2 = [gene{iGenes2}];
		x1(j,idcs1) = x0(iParents(j,1),idcs1);
		x1(j,idcs2) = x0(iParents(j,2),idcs2);

		for k = 1:nPars
			x1(j,k) = g{k}(x1(j,k),T(k));
		end
	end
	x0 = x1;
	x0(1,:) = best(:,i);

	% osc = nansum(abs(FBL(best(:,i),xData)-F(best(:,i),xData)));
	% hPlot1.YData = F(best(:,i),fValues);
	% hPlot2.YData = FBL(best(:,i),fValues);
	% yl = get(gca,'ylim');
	% txt.Position(2) = yl(1)+0.1*range(yl);
	% txt.String = ['i = ' int2str(i)];
	% drawnow;

end

c0 = best(:,end);
