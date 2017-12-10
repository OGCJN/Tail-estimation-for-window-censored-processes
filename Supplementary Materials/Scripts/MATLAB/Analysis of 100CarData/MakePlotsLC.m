%% Tail plots for left-censored glances for Gamma and Weibull full-parametric models, and GP and Exponential semi-parametric models, 
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

S_LC = S(S==0);
L_LC = L(S==0);

disp(['There are a total of ' num2str(numel(L_LC)) ' lc glances']);

%% Gamma Distribution
figure('Name','Gamma');
ecdf(L_LC,'function','survivor');
set(gca,'FontSize',13);
hold on;
xlabel('Off-road glance duration, $x$','interpreter','latex','fontsize',16);
ylabel('Residual tail function, $\bar{F}^r(x)$','interpreter','latex','fontsize',16);
title('LC Glances', 'fontsize',14,'fontweight','bold');

[gamma, sigma] = EstimateGamma(S,L,w);
x = 0:0.05:w;
Fbar = nan(1,numel(x));
for i=1:numel(x)
    Fbar(i) = gammainc(x(i)/sigma,gamma+1,'upper')-x(i)/(gamma*sigma)*gammainc(x(i)/sigma,gamma,'upper');
end

plot(x,Fbar,'r');

par = gamfit(L);

Fbar = nan(1,numel(x));
for i=1:numel(x)
    Fbar(i) = gammainc(x(i)/par(2),par(1)+1,'upper')-x(i)/(par(1)*par(2))*gammainc(x(i)/par(2),par(1),'upper');
end

plot(x,Fbar,':k');

%saveas(gcf,['Gam_LC_' name '.fig']);
%saveas(gcf,['Gam_LC_' name '.eps']); 
disp('Gamma Distribution: ');
disp(['Standard method:        shape=' num2str(par(1)) ' scale=' num2str(par(2))]);
disp(['Rootzen & Zholud, 2013: shape=' num2str(gamma)  ' scale=' num2str(sigma)]);

%% Weibull Distribution
figure('Name','Weibull');
ecdf(L_LC,'function','survivor');
set(gca,'FontSize',13);
hold on;
xlabel('Off-road glance duration, $x$','interpreter','latex','fontsize',16);
ylabel('Residual tail function, $\bar{F}^r(x)$','interpreter','latex','fontsize',16);
title('LC Glances', 'fontsize',14,'fontweight','bold');

[gamma, sigma] = EstimateWeibull(S,L,w);
x = 0:0.05:w;
Fbar = nan(1,numel(x));
for i=1:numel(x)
    Fbar(i) = gammainc((x(i)/sigma)^gamma,1/gamma,'lower');
end

plot(x,1-Fbar,'r','LineWidth',2);

par = wblfit(L);
Fbar = nan(1,numel(x));
for i=1:numel(x)
    Fbar(i) = gammainc((x(i)/par(1))^par(2),1/par(2),'lower');
end

plot(x,1-Fbar,':k','LineWidth',2);

%saveas(gcf,['Wbl_LC_' name '.fig']);
%saveas(gcf,['Wbl_LC_' name '.eps']); 
disp('Weibull Distribution: ');
disp(['Standard method:        shape=' num2str(par(2)) ' scale=' num2str(par(1))]);
disp(['Rootzen & Zholud, 2013: shape=' num2str(gamma)  ' scale=' num2str(sigma)]);

%% GP Distribution
figure('Name','GP');
ecdf(L_LC(L_LC>u)-u,'function','survivor');
set(gca,'FontSize',13);
hold on;
xlabel('Off-road glance exceedances, $x$','interpreter','latex','fontsize',16);
ylabel('Tail function, $\bar{F}(x)$','interpreter','latex','fontsize',16);
title('LC Glances', 'fontsize',14,'fontweight','bold');

[gamma, sigma] = EstimateCensoredSemiparametricGP(S,L,w,u);
x = 0:0.05:w-u;
Fbar = nan(1,numel(x));
for i=1:numel(x)
    Fbar(i) = max(1+gamma/sigma*x(i),0)^(1-1/gamma);
end

plot(x,Fbar,'r','LineWidth',2);

par = gpfit(L(L>u)-u);
Fbar = nan(1,numel(x));
for i=1:numel(x)
    Fbar(i) = max(1+par(1)/par(2)*x(i),0)^(1-1/par(1));
end

plot(x,Fbar,':k','LineWidth',2);

%saveas(gcf,['GPD_LC_' name '.fig']);
%saveas(gcf,['GPD_LC_' name '.eps']); 
disp('Generalized Pareto Distribution: ');
disp(['Standard method:        shape=' num2str(par(1)) ' scale=' num2str(par(2))]);
disp(['Rootzen & Zholud, 2013: shape=' num2str(gamma)  ' scale=' num2str(sigma)]);

%% Exponential Distribution
figure('Name','Exponential');
ecdf(L_LC(L_LC>u)-u,'function','survivor');
set(gca,'FontSize',13);
hold on;
xlabel('Off-road glance duration, $x$','interpreter','latex','fontsize',16);
ylabel('Residual tail function, $\bar{F}^r(x)$','interpreter','latex','fontsize',16);
title('LC Glances', 'fontsize',14,'fontweight','bold');

sigma = EstimateCensoredSemiparametricExponential(S,L,w,u);
x = 0:0.05:w-u;
plot(x,exp(-x/sigma),'r');

par = expfit(L(L>u)-u);

plot(x,exp(-x/par),':k');

%saveas(gcf,['Exp_LC_' name '.fig']);
%saveas(gcf,['Exp_LC_' name '.eps']); 
disp('Exponential Distribution: ');
disp(['Standard method:        mean=' num2str(par)]);
disp(['Rootzen & Zholud, 2013: mean=' num2str(sigma)]);
