function [F,FBL] = fittingmodel

gi = @(k,a1,a2) 1e3/(a1-a2)*((1e3/a1+sqrt(-1)*2*pi*k).^(-1)-(1e3/a2+sqrt(-1)*2*pi*k).^(-1));
lgNorm = @(k,m1,s1) exp(m1-s1^2/2)*exp(-(log(k)-m1).^2/(2*s1^2))./k;
gausses = @(m1,s1,b1,m2,s2,b2,m4,s4,b4,k) b1*lgNorm(k,m1,s1)+b2*lgNorm(k,m2,s2)+b4*lgNorm(k,m4,s4);
gz = @(k,a1,n) 1./(1+2*pi*sqrt(-1)*(k/a1).^n);

ftfun = @(m1,s1,m2,s2,m4,s4,b1,b2,b4,Z1,gamI,tauI,tauI2,wf,theta,f) 2* log(abs(gz(f,wf,theta))) +  ...
	log( Z1 + gamI*(1+gamI/2*gausses(m1,s1,b1,m2,s2,b2,m4,s4,b4,f)).* abs(gi(f,tauI,tauI2)).^2);
F = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),f);



FBL = @(x,f) ftfun(x(1),x(2),x(3),x(4),x(5),x(6),0,0,0,x(10),x(11),x(12),x(13),x(14),x(15),f);
