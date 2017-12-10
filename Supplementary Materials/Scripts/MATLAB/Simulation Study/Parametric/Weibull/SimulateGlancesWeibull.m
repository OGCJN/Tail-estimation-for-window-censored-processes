%% Simulate censored glances under full-parametric Weibull model
%   INPUT:  mu1 - mean of the Exponential distribution, F1
%           k0 - shape parameter of the Weibull distribution, F0 
%           sigma0 - scale parameter of the Weibull distribution, F0
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
function [S, L]=SimulateGlancesWeibull(N,k0,sigma0,mu1,w)
mu0=sigma0*gamma(1+1/k0);
p0=mu0/(mu0+mu1);
S=nan(1,N);
L=nan(1,N);
for i=1:N
   while true
      p=rand;
      if p<=p0 % Start with the zero glance 
         S(i)=0;
         L(i)=min(w,SimFromWeibullOvershoot(k0,sigma0)); % it will be lc or dc
      else % Start with the one glance
         S(i)=exprnd(mu1);
         if S(i)>w, continue; end; % if it exceeds w then there are no zero glances and we continue
         L(i)=min(w-S(i),wblrnd(sigma0,k0)); % else we get the zero glance - it can be nc or rc
      end
      break; % if we get here then we got a zero glance
   end
end
end

%% Simulating random variables from the corresponding overshoot distribution 
function y = SimFromWeibullOvershoot(k0,sigma0)
    u=rand;
    y=sigma0*(gammaincinv(u,1/k0,'lower')^(1/k0));
end