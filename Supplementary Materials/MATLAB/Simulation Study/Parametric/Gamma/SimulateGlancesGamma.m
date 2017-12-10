%% Simulate censored glances under full-parametric Gamma model
%   INPUT:  mu1 - mean of the Exponential distribution, F1
%           k0 - shape parameter of the Gamma distribution, F0 
%           sigma0 - scale parameter of the Gamma distribution, F0
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
function [S, L]=SimulateGlancesGamma(N,k0,sigma0,mu1,w)
mu0=k0*sigma0;
p0=mu0/(mu0+mu1);
S=nan(1,N);
L=nan(1,N);
for i=1:N
   while true
      p=rand;
      if p<=p0 % Start with the zero glance 
         S(i)=0;
         L(i)=min(w,SimFromGammaOvershoot(k0,sigma0)); % it will be lc or dc
      else % Start with the one glance
         S(i)=exprnd(mu1);
         if S(i)>w, continue; end; % if it exceeds w then there are no zero glances and we continue
         L(i)=min(w-S(i),gamrnd(k0,sigma0)); % else we get the zero glance - it can be nc or rc
      end
      break; % if we get here then we got a zero glance
   end
end
end

%% Simulating random variables from the corresponding overshoot distribution 
function y = SimFromGammaOvershoot(k0,sigma0)
    C1=1/k0;
    C2=gamma(k0+1);
    L=@(x) (1-x*C1)*gammainc(x,k0,'upper')+(x)^k0*exp(-x)/C2;
    u=rand;
    y=sigma0*fzero(@(x) L(x)-u,[0 2*k0*max(-(1/k0*log(u*k0*gamma(k0))-log(2*k0)), max(1/sigma0,1))]);
end