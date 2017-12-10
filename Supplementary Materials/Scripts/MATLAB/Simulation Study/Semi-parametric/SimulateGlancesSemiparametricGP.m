%% Simulate censored glances under semi-parametric Generalized Pareto model
%   INPUT:  mu1 - mean of the Exponential distribution, F1
%           u - threshold that spearates the uniform component from the GP
%               component.
%           gamma - shape parameter of the GP distribution of exceedances of the threshold u.
%                   Note that scale parameter is computed automatically from the expressions for the mixture components and the requirement 
%                   that the resulting distribution should be continuous.
%           p - mixture probability
%           N - desired number of observed censored glances
%           w - length of the observation window 
%
%   OUTPUT: S - 1xN vector of glance starting positions relative to the beginning of the observation window. 
%           For left-censored glances we put S = 0.
%           L - 1xN vector of glance durations. 
%           Each glance k is a pair (S(k),L(k)), k=1,2,...,N.
%
%   Ver 1.0, 23 Jan 2014
%   Supplement to "Tail Estimation for Window Censored Processes"
function [S, L]=SimulateGlancesSemiparametricGP(N,u,p,gamma,mu1,w)
% Sub-routine used in simulation from mixture distribution algorithm 
Choose = @(x0,p0,A0,B0) (x0<=p0).*A0+(x0>p0).*B0;

% Simulating random variables from F0 distribution 
SemiGPrnd = @(u0,p0,gamma0,N0) Choose(rand(1,N0),p0,...
                               u0*rand(1,N0),...
                               gprnd(gamma0,u0*(1-p0)/p0,u0,1,N0));
                           
% Simulating random variables from the corresponding overshoot distribution                            
SemiGPOvershootrnd = @(u0,p0,gamma0,N0) Choose(rand(1,N0),...
    u0*(1-p0/2)/(u0*(1-p0/2)+(1-p0)^2*u0/(p0*(1-gamma0))),...
    u0/p0-sqrt((u0/p0)^2-2*u0/p0*(u0*(1-p0/2))*rand(1,N0)),...
    gprnd(gamma0/(1-gamma0),u0*(1-p0)/(p0*(1-gamma0)),u0,1,N0));

mu0=u*(1-p/2)+(1-p)^2*u/(p*(1-gamma));
p0=mu0/(mu0+mu1);
S=nan(1,N);
L=nan(1,N);
for i=1:N
   while true
      pp=rand;
      if pp<=p0 % Start with the zero glance 
         S(i)=0;
         L(i)=min(w,SemiGPOvershootrnd(u,p,gamma,1)); % it will be lc or dc
      else % Start with the one glance
         S(i)=exprnd(mu1);
         if S(i)>w, continue; end; % if it exceeds w then there are no zero glances and we continue
         L(i)=min(w-S(i),SemiGPrnd(u,p,gamma,1)); % else we get the zero glance - it can be nc or rc
      end
      break; % if we get here then we got a zero glance
   end
end
