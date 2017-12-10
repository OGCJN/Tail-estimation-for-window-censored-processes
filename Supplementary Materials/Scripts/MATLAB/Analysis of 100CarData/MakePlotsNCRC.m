%% Tail plots for non-censored and right-censored glances for Gamma and Weibull full-parametric models, and GP and Exponential semi-parametric models, 
%  see Figure 3 and Figure 4. Use 'Merged_200ms' data for glance durations and Merged_1000ms for task durations
%  Ver 1.0, 23 Jan 2014
%  Supplement to "Tail Estimation for Window Censored Processes"
name = 'Merged_200ms'; 
addpath(genpath('..\Estimators\'));

w=6;
u=2;

load(['Extracted Glance Data\S_' name '.mat']);
load(['Extracted Glance Data\L_' name '.mat']);

S(L==0)=[];
L(L==0)=[];

S_NCRC = S(S>0);
L_NCRC = L(S>0);

censored = (S_NCRC+L_NCRC==w);

disp(['There are a total of ' num2str(numel(L_NCRC)) ' nc+rc glances']);
disp(['Out of which ' num2str(sum(censored)) ' are rc glances']);

%% Gamma Distribution
figure('Name','Gamma');
ecdf(L_NCRC,'function','survivor','censoring',censored);
hold on;
xlabel('Off-road glance duration, $x$','interpreter','latex','fontsize',16);
ylabel('Residual tail function, $\bar{F}^r(x)$','interpreter','latex','fontsize',16);
title('NC+RC Glances', 'fontsize',14,'fontweight','bold');

[gamma, sigma] = EstimateGamma(S,L,w);
x = 0:0.05:w;
plot(x,1-cdf('Gamma',x,gamma,sigma),'r');

par = gamfit(L);

plot(x,1-cdf('gamma',x,par(1),par(2)),':k');

%saveas(gcf,['Gam_NCRC_' name '.fig']);
%saveas(gcf,['Gam_NCRC_' name '.eps']); 
disp('Gamma Distribution: ');
disp(['Standard method:        shape=' num2str(par(1)) ' scale=' num2str(par(2))]);
disp(['Rootzen & Zholud, 2013: shape=' num2str(gamma)  ' scale=' num2str(sigma)]);

%% Weibull Distribution
figure('Name','Weibull');
ecdf(L_NCRC,'function','survivor','censoring',censored);
set(gca,'FontSize',13);
hold on;
xlabel('Off-road glance duration, $x$','interpreter','latex','fontsize',16);
ylabel('Residual tail function, $\bar{F}^r(x)$','interpreter','latex','fontsize',16);
title('NC+RC Glances', 'fontsize',14,'fontweight','bold');

[gamma, sigma] = EstimateWeibull(S,L,w);
x = 0:0.05:w;
plot(x,1-cdf('Weibull',x,sigma,gamma),'r','LineWidth',2);

par = wblfit(L);

plot(x,1-cdf('Weibull',x,par(1),par(2)),':k','LineWidth',2);

%saveas(gcf,['Wbl_NCRC_' name '.fig']);
%saveas(gcf,['Wbl_NCRC_' name '.eps']); 
disp('Weibull Distribution: ');
disp(['Standard method:        shape=' num2str(par(2)) ' scale=' num2str(par(1))]);
disp(['Rootzen & Zholud, 2013: shape=' num2str(gamma)  ' scale=' num2str(sigma)]);

%% GP Distribution
figure('Name','GPD');
ecdf(L_NCRC(L_NCRC>u)-u,'function','survivor','censoring',censored(L_NCRC>u));
set(gca,'FontSize',13);
hold on;
xlabel('Off-road glance exceedances, $x$','interpreter','latex','fontsize',16);
ylabel('Tail function, $\bar{F}(x)$','interpreter','latex','fontsize',16);
title('NC+RC Glances', 'fontsize',14,'fontweight','bold');

[gamma, sigma] = EstimateCensoredSemiparametricGP(S,L,w,u);
x = 0:0.05:w-u;
plot(x,1-cdf('Generalized Pareto',x,gamma,sigma,0),'r','LineWidth',2);

par = gpfit(L(L>u)-u);

plot(x,1-cdf('Generalized Pareto',x,par(1),par(2)),':k','LineWidth',2);

%saveas(gcf,['GPD_NCRC_' name '.fig']);
%saveas(gcf,['GPD_NCRC_' name '.eps']); 
disp('Generalized Pareto Distribution: ');
disp(['Standard method:        shape=' num2str(par(1)) ' scale=' num2str(par(2))]);
disp(['Rootzen & Zholud, 2013: shape=' num2str(gamma)  ' scale=' num2str(sigma)]);

%% Exponential Distribution
figure('Name','Exponential');
ecdf(L_NCRC(L_NCRC>u)-u,'function','survivor','censoring',censored(L_NCRC>u));
hold on;
xlabel('Off-road glance duration, $x$','interpreter','latex','fontsize',16);
ylabel('Residual tail function, $\bar{F}^r(x)$','interpreter','latex','fontsize',16);
title('NC+RC Glances', 'fontsize',14,'fontweight','bold');

sigma = EstimateCensoredSemiparametricExponential(S,L,w,u);
x = 0:0.05:w-u;
plot(x,1-cdf('Exponential',x,sigma),'r');

par = expfit(L(L>u)-u);

plot(x,1-cdf('Exponential',x,par),':k');

%saveas(gcf,['Exp_NCRC_' name '.fig']);
%saveas(gcf,['Exp_NCRC_' name '.eps']); 
disp('Exponential Distribution: ');
disp(['Standard method:        mean=' num2str(par)]);
disp(['Rootzen & Zholud, 2013: mean=' num2str(sigma)]);
