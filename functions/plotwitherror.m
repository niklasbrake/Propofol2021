function h =plotwitherror(x,y,isstdrr,varargin)

	if(strcmp(isstdrr,'CI') || isstdrr == true)
		ymu = nanmean(y,2);
		% ymu = nanmedian(y,2);
		z = icdf('norm',0.975,0,1);
		ylo = ymu-z*stderror(y')';
		yhi = ymu+z*stderror(y')';
	elseif(strcmp(isstdrr,'SE') || isstdrr == false)
		ymu = nanmean(y,2);
		ylo = ymu-stderror(y')';
		yhi = ymu+stderror(y')';
	elseif(strcmp(isstdrr,'Q'))
		ymu = nanmedian(y,2);
		ylo = quantile(y,0.025,2);
		yhi = quantile(y,0.925,2);
	end

	h = plot(x,ymu,varargin{:}); hold on;

	idcs = ~isnan(ymu);
	ylo = ylo(idcs);
	yhi = yhi(idcs);
	ymu = ymu(idcs);
	x = x(idcs);

	xfvec = [x,flip(x)];
	yfvec = [ylo;flip(yhi)];
	F = fill(xfvec(:),yfvec(:),'k');

	if(isempty(F))
		return;
	end
	F.FaceAlpha = 0.2;
	F.FaceColor = h.Color;
	F.EdgeColor = 'none';
