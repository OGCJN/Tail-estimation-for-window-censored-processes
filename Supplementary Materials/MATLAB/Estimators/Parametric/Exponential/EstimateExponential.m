%% Estimation of the mean of the Exponential distribution for glance durations under full-parametric censoring model
%   INPUT:  S - glance start, relative to the observation window.
%               For left-censored glances put S = 0.
%           L - glance duration
%           w - length of the observation window 
%
%   OUTPUT: sigma - mean parameter of the Exponential distribution
%           CI    - 1x2 vector that constitutes 95% confidence interval for sigma 
%
%   Ver 1.0, 23 Jan 2014
%   Supplement to "Tail Estimation for Window Censored Processes"

function [sigma, CI]=EstimateExponential(S,L,w)
nlc=sum((S==0)&(S+L<w));
nnc=sum((S>0)&(S+L<w));
sigma=sum(L)/(nlc+nnc);  
C=norminv(0.975,0,1)*sigma/sqrt(nlc+nnc);  
CI=[sigma-C sigma+C];
end