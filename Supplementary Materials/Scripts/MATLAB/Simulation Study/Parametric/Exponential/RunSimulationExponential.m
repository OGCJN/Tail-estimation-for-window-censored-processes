%% Simulation study for the full-parametric Exponential model.
% OBS: takes approximately a day to run.
%
%   Ver 1.0, 23 Jan 2014
%   Supplement to "Tail Estimation for Window Censored Processes"

w   = 6;     % Window length will be fixed
mu1 = 6;     % And F1 will be exponential with mean 6
m   = 10000; % Number of simulations is 10,000. 

Ns   = [50 250 1000];  % Sample sizes, N
lambdas =  [0.2 1 5];  % Means of the exponential distribution F0

addpath(genpath('..\..\..\Estimators\'));

figure();
set(gcf,'units','pixels');
set(gcf,'position',[404         214        1300        1124]);

for i=1:numel(Ns)
    for j=1:numel(lambdas)
        lambda_hat=zeros(1,m);
        CI=zeros(2,m);
        N=Ns(i);
        lambda=lambdas(j);
        parfor k=1:m
            [S,L]=SimulateGlancesExponential(N,lambda,mu1,w);     % Simulation step
            [lambda_hat(k), CI(:,k)]=EstimateExponential(S,L,w);  % Estimation step
        end

        % We will make histograms of the estimates
        subplot(numel(Ns),numel(lambdas),(i-1)*numel(lambdas)+j)
        hist(lambda_hat);
        hold on;
        YL=get(gca,'YLim');
        plot([lambda lambda],[YL(1) YL(2)],'--r');
        
        % Compute errors and coverge probabilities
        bias_lambda=mean(lambda_hat)-lambda;
        std_lambda=std(lambda_hat);
        rmse_lambda=sqrt(1/m*sum((lambda_hat-lambda).^2));
        CP_lambda=sum((CI(1,:)<lambda)&(CI(2,:)>lambda))/m;
        
        % Finalize our plot
        plot([lambda+bias_lambda lambda+bias_lambda],[YL(1) YL(2)*0.99],'green','LineWidth',1);
        
        box on;
        title(['$\lambda=$ ' num2str(lambda) ', $N$= ' num2str(N)], 'interpreter','latex');
        xlabel('Estimates, $\hat{\lambda}$','interpreter','latex');
        ylabel('Frequency','interpreter','latex');
        
        % Display results
        disp(['N = ' num2str(N) ', Lambda= ' num2str(lambda) ', Bias: ' num2str(bias_lambda)]);
        disp(['N = ' num2str(N) ', Lambda= ' num2str(lambda) ', STD : ' num2str(std_lambda)]);
        disp(['N = ' num2str(N) ', Lambda= ' num2str(lambda) ', RMSE: ' num2str(rmse_lambda)]);
        disp(['N = ' num2str(N) ', Lambda= ' num2str(lambda) ', CP  : ' num2str(CP_lambda)]);
        disp('----------------------------------------');
    end
end