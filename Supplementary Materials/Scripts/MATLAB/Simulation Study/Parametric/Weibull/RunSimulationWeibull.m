%% Simulation study for the full-parametric Weibull model.
% OBS: takes approximately a day to run.
%
%   Ver 1.0, 23 Jan 2014
%   Supplement to "Tail Estimation for Window Censored Processes"

w   = 6;     % Window length will be fixed
mu1 = 6;     % And F1 will be exponential with mean 6
m   = 10000; % Number of simulations is 10,000. 

Ns     = [50 250 1000];            % Sample sizes, N
k      = 1/2;                      % Shape of F0
sigmas = [0.2, 1, 5]/gamma(1+1/k); % Scale of F0

addpath(genpath('..\..\..\Estimators\'));

figure();
set(gcf,'units','pixels');
set(gcf,'position',[404         214        1300        1124]);

for i=1:numel(Ns)
    for j=1:numel(sigmas)
        k_hat     = zeros(1,m);
        sigma_hat = zeros(1,m);
        CI_k      = zeros(2,m);
        CI_sigma  = zeros(2,m);

        N    = Ns(i);
        sigma = sigmas(j);
        
        parfor p=1:m
            [S,R]=SimulateGlancesWeibull(N,k,sigma,mu1,w);                             % Simulation step
            [k_hat(p), sigma_hat(p), CI_k(:,p), CI_sigma(:,p)]=EstimateWeibull(S,R,w); % Estimation step
        end

        % We will make histograms of the estimates of k
        subplot(numel(Ns),2*numel(sigmas),(i-1)*2*numel(sigmas)+2*(j-1)+1);
        
        hist(k_hat);
        hold on;
        YL=get(gca,'YLim');
        plot([k k],[YL(1) YL(2)],'--r');
        
        % Compute errors and coverge probabilities
        bias_k = mean(k_hat)-k;
        std_k  = std(k_hat);
        rmse_k = sqrt(1/m*sum((k_hat-k).^2));
        CP_k   = sum((CI_k(1,:)<k)&(CI_k(2,:)>k))/m;
        
        % Finalize our plot
        plot([k+bias_k k+bias_k],[YL(1) YL(2)*0.99],'green','LineWidth',1);
        
        box on;
        title(['$k=$ ' num2str(k) ', $\sigma$= ' num2str(sigma) ', $N$= ' num2str(N)], 'interpreter','latex');
        xlabel('Estimates, $\hat{k}$','interpreter','latex');
        ylabel('Frequency','interpreter','latex');
        
        % ...and now, similarly, histograms of the estimates of sigma
        subplot(numel(Ns),2*numel(sigmas),(i-1)*2*numel(sigmas)+2*(j-1)+2);
        
        hist(sigma_hat);
        hold on;
        YL=get(gca,'YLim');
        plot([sigma sigma],[YL(1) YL(2)],'--r');
        
        % Compute errors and coverge probabilities
        bias_sigma = mean(sigma_hat)-sigma;
        std_sigma  = std(sigma_hat);
        rmse_sigma = sqrt(1/m*sum((sigma_hat-sigma).^2));
        CP_sigma   = sum((CI_sigma(1,:)<sigma)&(CI_sigma(2,:)>sigma))/m;
        
        % Finalize our plot
        plot([sigma+bias_sigma sigma+bias_sigma],[YL(1) YL(2)*0.99],'green','LineWidth',1);
        
        box on;
        title(['$k=$ ' num2str(k) ', $\sigma$= ' num2str(sigma) ', $N$= ' num2str(N)], 'interpreter','latex');
        xlabel('Estimates, $\hat{\sigma}$','interpreter','latex');
        ylabel('Frequency','interpreter','latex');       
        
        % Display results
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', Bias(k): ' num2str(bias_k)]);
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', STD(k) : ' num2str(std_k)]);
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', RMSE(k): ' num2str(rmse_k)]);
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', CP(k): '   num2str(CP_k)]);

        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', Bias(sigma): ' num2str(bias_sigma)]);
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', STD(sigma) : ' num2str(std_sigma)]);
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', RMSE(sigma): ' num2str(rmse_sigma)]);
        disp(['N = ' num2str(N) ', k= ' num2str(k) ', sigma= ' num2str(sigma) ', CP(sigma): '   num2str(CP_sigma)]);
        disp('----------------------------------------');
    end
end