%% MACROECONOMETRICS - PROBLEM SET 1
% Complete MATLAB solution for Questions 1-6.
%
% This script follows the simple commands and matrix methods used in
% Professor Luca Sala's first MATLAB session.  In particular, it uses:
%   - randn and rng to generate reproducible Gaussian random variables;
%   - explicit for loops to show the recursions clearly;
%   - filter to implement AR and MA lag polynomials;
%   - X\Y, inv(X'*X), and cov(residuals) for OLS estimation;
%   - histogram and subplot for the Monte Carlo figures.
%
% Run this file from any folder.  It finds its own location and saves every
% table, figure, and data file in the same "code and outputs" directory.
%
% Important convention used below:
% An AR(1) with unconditional mean mu is written as
%
%       y_t = c + phi*y_(t-1) + epsilon_t,
%
% where c = mu*(1-phi).  Hence E(y_t) = mu when |phi| < 1.

clear;
clc;
close all;

%% Output directory and reproducibility information

% mfilename('fullpath') returns the full path of this script while it runs.
script_full_path = mfilename('fullpath');

% Remove the script's file name from that full path.  The remaining path is
% the "code" folder that contains both MATLAB files.
code_dir = fileparts(script_full_path);

% The main PS1 folder is one level above the "code" folder.
ps1_dir = fileparts(code_dir);

% Store all generated tables, figures, and data in a separate folder named
% "outputs" inside the main PS1 folder.
output_dir = fullfile(ps1_dir,'outputs');

% Create the outputs folder automatically the first time the script runs.
% If the folder already exists, MATLAB simply skips this block.
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Make the code folder visible to MATLAB.  This ensures that the separate
% arma_generator_ps1.m function can be found even when MATLAB's Current
% Folder is somewhere else.
addpath(code_dir);

% Record Command Window output.  Delete an older log first so that each run
% produces a clean, reproducible file rather than appending to an old run.
log_file = fullfile(output_dir,'PS1_command_window_log.txt');
if exist(log_file,'file')
    delete(log_file);
end
diary(log_file);

disp('Macroeconometrics Problem Set 1');
disp(['Output directory: ' output_dir]);
disp(' ');

%% Question 1(a): Generate an AR(1) with a for loop
% The requested process has T=500, E(y_t)=0, phi=0.4, and innovation
% variance sigma_e^2=0.2.  We choose y_0=0 so that the loop and filter have
% exactly the same starting condition.

T1 = 500;                 % Number of observations requested by PS1.
phi1 = 0.4;               % AR(1) coefficient.
sigma2_e1 = 0.2;          % Variance of the white-noise innovation.
mu1 = 0;                  % Unconditional mean of the process.
c1 = mu1*(1-phi1);        % Intercept implied by the unconditional mean.
y0_1 = 0;                 % Starting value y_0.

% Fix the random-number seed so that the results are reproducible.
rng(101,'twister');

% randn(T1,1) makes a T1-by-1 vector of standard-normal draws.  Multiplying
% by sqrt(sigma2_e1) gives innovations with variance sigma2_e1.
epsilon1 = sqrt(sigma2_e1)*randn(T1,1);

% Pre-allocate a T1-by-1 vector.  This is clearer and faster than enlarging
% the vector inside the loop.
y1_loop = zeros(T1,1);

% The variable y_previous stores the value dated one period before t.
y_previous = y0_1;

% Apply y_t = c + phi*y_(t-1) + epsilon_t for every observation.
for t = 1:T1
    y1_loop(t) = c1 + phi1*y_previous + epsilon1(t);
    y_previous = y1_loop(t);
end

%% Question 1(b): Generate the same AR(1) with filter
% filter(b,a,x) applies the lag-polynomial relation a(L)y_t=b(L)x_t.
% Therefore filter(1,[1 -phi1],epsilon1) solves
% (1-phi1*L)y_t=epsilon_t with a zero initial condition.

y1_filter = filter(1,[1 -phi1],epsilon1);

%% Question 1(c): Check that the two methods give the same output

% Use the same innovations and initial condition, then compare period by
% period.  Differences should be zero up to numerical rounding error.
difference1 = y1_loop-y1_filter;
max_abs_difference1 = max(abs(difference1));
rmse_difference1 = sqrt(mean(difference1.^2));

fprintf('QUESTION 1\n');
fprintf('Maximum absolute loop-filter difference: %.16g\n',max_abs_difference1);
fprintf('Root mean squared loop-filter difference: %.16g\n\n',rmse_difference1);

% Save the complete simulated paths so they can be inspected outside MATLAB.
question1_paths = table((1:T1)',epsilon1,y1_loop,y1_filter,difference1, ...
    'VariableNames',{'t','epsilon','y_loop','y_filter','difference'});
writetable(question1_paths,fullfile(output_dir,'question1_simulated_series.csv'));

% Save the two numerical comparison measures in a compact summary table.
question1_summary = table( ...
    ["Maximum absolute difference";"Root mean squared difference"], ...
    [max_abs_difference1;rmse_difference1], ...
    'VariableNames',{'Statistic','Value'});
writetable(question1_summary,fullfile(output_dir,'question1_summary.csv'));

% Plot the two paths together and then plot their numerical difference.
figure1 = figure('Visible','off','Color','w','Position',[100 100 1000 650]);
subplot(2,1,1);
plot(1:T1,y1_loop,'b-','LineWidth',1.1);
hold on;
plot(1:T1,y1_filter,'r--','LineWidth',1.0);
grid on;
xlabel('t');
ylabel('y_t');
title('Question 1: AR(1) generated by a loop and by filter');
legend('For loop','filter','Location','best');

subplot(2,1,2);
plot(1:T1,difference1,'k-','LineWidth',1.0);
grid on;
xlabel('t');
ylabel('y_t^{loop} - y_t^{filter}');
title('Period-by-period difference between the two methods');

exportgraphics(figure1,fullfile(output_dir,'question1_loop_vs_filter.png'), ...
    'Resolution',300);
close(figure1);

%% Question 2: Initial condition and burn-in for a stationary AR(1)
% The requested process has phi=0.6, innovation variance 0.4, and
% unconditional mean 3.  The starting condition y_0=20 is far above the
% mean.  Its expected effect after t periods is phi^t*(y_0-mu), which goes
% to zero because |phi|<1.

T2 = 500;                 % Number of observations retained at the end.
burn2 = 500;              % Initial observations discarded as burn-in.
phi2 = 0.6;               % AR coefficient.
sigma2_e2 = 0.4;          % Innovation variance.
mu2 = 3;                  % Unconditional mean.
c2 = mu2*(1-phi2);        % Intercept that makes E(y_t)=3.
y0_2 = 20;                % Starting value specified in PS1.

rng(202,'twister');

% Generate enough innovations for both the burn-in and the retained sample.
epsilon2_all = sqrt(sigma2_e2)*randn(T2+burn2,1);
y2_all = zeros(T2+burn2,1);
y_previous = y0_2;

% Generate the complete path from the deliberately distant starting value.
for t = 1:(T2+burn2)
    y2_all(t) = c2 + phi2*y_previous + epsilon2_all(t);
    y_previous = y2_all(t);
end

% Discard the first burn2 observations.  The remaining T2 observations are
% effectively independent of the arbitrary initial condition.
y2_stationary = y2_all(burn2+1:burn2+T2);

% The stationary variance of an AR(1) is sigma_e^2/(1-phi^2).
theoretical_variance2 = sigma2_e2/(1-phi2^2);

% The deterministic remaining effect of y_0 after the burn-in is tiny.
remaining_initial_effect2 = (phi2^burn2)*(y0_2-mu2);

fprintf('QUESTION 2\n');
fprintf('Mean of retained post-burn-in observations: %.6f\n',mean(y2_stationary));
fprintf('Variance of retained observations: %.6f\n',var(y2_stationary));
fprintf('Theoretical stationary variance: %.6f\n',theoretical_variance2);
fprintf('Initial-condition effect after burn-in: %.16g\n\n',remaining_initial_effect2);

% Save the retained sample path and its corresponding innovations.
question2_path = table((1:T2)',y2_stationary,epsilon2_all(burn2+1:end), ...
    'VariableNames',{'t','y_after_burnin','epsilon'});
writetable(question2_path,fullfile(output_dir,'question2_stationary_path.csv'));

% Save the quantities useful for discussing the initial-condition problem.
question2_summary = table( ...
    ["Unconditional mean";"Mean after burn-in"; ...
     "Theoretical stationary variance";"Sample variance after burn-in"; ...
     "Burn-in observations";"Remaining initial-condition effect"], ...
    [mu2;mean(y2_stationary);theoretical_variance2;var(y2_stationary); ...
     burn2;remaining_initial_effect2], ...
    'VariableNames',{'Statistic','Value'});
writetable(question2_summary,fullfile(output_dir,'question2_summary.csv'));

% The first panel emphasizes the rapid movement away from y_0=20.  The
% conditional expected path is mu+phi^t*(y_0-mu).  The second panel shows
% the observations retained after burn-in.
figure2 = figure('Visible','off','Color','w','Position',[100 100 1000 650]);
subplot(2,1,1);
initial_periods = (0:50)';
expected_initial_path2 = mu2+(phi2.^initial_periods)*(y0_2-mu2);
plot(initial_periods,[y0_2;y2_all(1:50)],'b-','LineWidth',1.1);
hold on;
plot(initial_periods,expected_initial_path2,'r--','LineWidth',1.4);
yline(mu2,'k:','Unconditional mean','LineWidth',1.1);
grid on;
xlabel('t');
ylabel('y_t');
title('Question 2: Effect of starting at y_0 = 20');
legend('Simulated path','Conditional expected path','Location','best');

subplot(2,1,2);
plot(1:T2,y2_stationary,'b-','LineWidth',1.0);
hold on;
yline(mu2,'k--','Unconditional mean','LineWidth',1.1);
grid on;
xlabel('Retained observation');
ylabel('y_t');
title(['Stationary realization after discarding ' num2str(burn2) ' observations']);

exportgraphics(figure2,fullfile(output_dir,'question2_initial_condition_and_burnin.png'), ...
    'Resolution',300);
close(figure2);

%% Question 3: Empirical distribution of OLS and t-test for phi=0
% Data-generating process:
%       y_t = phi*y_(t-1) + epsilon_t, phi=0.4.
% We choose Var(epsilon_t)=1, which PS1 leaves free.  We use 5,000 Monte
% Carlo replications, as in the previous official solution, and a burn-in
% so every retained sample is a proper stationary realization.

T3 = 250;                 % Observations in every retained sample.
burn3 = 300;              % Observations discarded before estimation.
phi3 = 0.4;               % True AR coefficient.
sigma2_e3 = 1;            % Chosen innovation variance.
MC3 = 5000;               % Number of Monte Carlo replications.
critical_value3 = 1.96;   % Two-sided 5% normal critical value.

rng(303,'twister');

% Pre-allocate one estimated coefficient and one t-statistic per replication.
phi_hat3 = zeros(MC3,1);
t_stat3 = zeros(MC3,1);

for simulation = 1:MC3
    % Generate T3+burn3 innovations, simulate with filter, and discard the
    % first burn3 observations to remove the zero starting condition.
    epsilon3 = sqrt(sigma2_e3)*randn(T3+burn3,1);
    y3_all = filter(1,[1 -phi3],epsilon3);
    y3 = y3_all(burn3+1:end);

    % The dependent variable contains y_2,...,y_T.
    Y3 = y3(2:end);

    % X3 has T3-1 rows and two columns.  The first column of ones includes
    % an intercept; the second column contains y_1,...,y_(T-1).
    X3 = [ones(T3-1,1) y3(1:end-1)];

    % Estimate beta=[intercept;phi] using the same backslash method used in
    % Professor Sala's regression_luca function.
    beta_hat3 = X3\Y3;

    % Compute residuals and, following regression_luca, use cov(residuals)
    % as the estimated innovation variance.
    residual3 = Y3-X3*beta_hat3;
    sigma2_hat3 = cov(residual3);

    % The estimated covariance matrix is sigma_hat^2*(X'X)^(-1).
    covariance_beta3 = inv(X3'*X3)*sigma2_hat3;
    standard_error_phi3 = sqrt(covariance_beta3(2,2));

    % Store the estimated AR coefficient and the t-test statistic for
    % H0: phi=0 against the two-sided alternative phi~=0.
    phi_hat3(simulation) = beta_hat3(2);
    t_stat3(simulation) = (beta_hat3(2)-0)/standard_error_phi3;
end

% The rejection frequency is the proportion of simulations with |t|>1.96.
reject3 = abs(t_stat3)>critical_value3;
rejection_rate3 = mean(reject3);

mean_phi_hat3 = mean(phi_hat3);
sd_phi_hat3 = std(phi_hat3);
bias_phi_hat3 = mean_phi_hat3-phi3;

fprintf('QUESTION 3\n');
fprintf('Monte Carlo replications: %d\n',MC3);
fprintf('Mean of phi_hat: %.6f\n',mean_phi_hat3);
fprintf('Standard deviation of phi_hat: %.6f\n',sd_phi_hat3);
fprintf('Bias of phi_hat: %.6f\n',bias_phi_hat3);
fprintf('Rejection rate for H0: phi=0 at 5%%: %.6f\n\n',rejection_rate3);

% Save every Monte Carlo estimate and test result.
question3_draws = table((1:MC3)',phi_hat3,t_stat3,reject3, ...
    'VariableNames',{'Simulation','phi_hat','t_statistic','reject_H0'});
writetable(question3_draws,fullfile(output_dir,'question3_monte_carlo_draws.csv'));

% Save the main results needed for the explanatory note.
question3_summary = table(T3,burn3,MC3,phi3,sigma2_e3,mean_phi_hat3, ...
    sd_phi_hat3,bias_phi_hat3,critical_value3,rejection_rate3, ...
    'VariableNames',{'T','BurnIn','Simulations','TruePhi','InnovationVariance', ...
    'MeanPhiHat','StdPhiHat','Bias','CriticalValue','RejectionRate'});
writetable(question3_summary,fullfile(output_dir,'question3_summary.csv'));

% Plot the empirical distribution and the large-sample normal approximation
% N(phi,(1-phi^2)/T).  The normal density is calculated manually to avoid
% requiring the Statistics and Machine Learning Toolbox.
figure3 = figure('Visible','off','Color','w','Position',[100 100 900 550]);
histogram3 = histogram(phi_hat3,35,'Normalization','pdf', ...
    'FaceColor',[0.35 0.60 0.85]);
hold on;
x_grid3 = linspace(min(phi_hat3),max(phi_hat3),300);
asymptotic_sd3 = sqrt((1-phi3^2)/T3);
normal_density3 = exp(-0.5*((x_grid3-phi3)/asymptotic_sd3).^2) ...
    /(sqrt(2*pi)*asymptotic_sd3);
normal_line3 = plot(x_grid3,normal_density3,'r-','LineWidth',1.6);
true_line3 = xline(phi3,'k--','LineWidth',1.2);
mean_line3 = xline(mean_phi_hat3,'b:','LineWidth',1.4);
grid on;
xlabel('Estimated phi');
ylabel('Density');
title('Question 3: Empirical distribution of the OLS estimator');
legend([histogram3 normal_line3 true_line3 mean_line3], ...
    {'Monte Carlo distribution','Asymptotic normal approximation', ...
    'True phi','Monte Carlo mean'},'Location','best');

exportgraphics(figure3,fullfile(output_dir,'question3_ols_distribution.png'), ...
    'Resolution',300);
close(figure3);

%% Question 4: OLS distribution for phi=0.9 as T increases
% For each requested sample size, run 1,000 simulations.  The same burn-in
% principle from Question 2 is used because phi=0.9 is persistent and a
% zero starting value would otherwise influence short samples.

phi4 = 0.9;
sigma2_e4 = 1;
sample_sizes4 = [50 100 200 1000];
MC4 = 1000;
burn4 = 500;

rng(404,'twister');

% Each column will contain the 1,000 estimates for one sample size.
phi_hat4_all = zeros(MC4,length(sample_sizes4));

% Pre-allocate the summary measures, one value for every sample size.
mean_phi_hat4 = zeros(length(sample_sizes4),1);
sd_phi_hat4 = zeros(length(sample_sizes4),1);
bias_phi_hat4 = zeros(length(sample_sizes4),1);
asymptotic_sd4 = zeros(length(sample_sizes4),1);

figure4 = figure('Visible','off','Color','w','Position',[100 100 1100 750]);

for size_index = 1:length(sample_sizes4)
    T4 = sample_sizes4(size_index);

    for simulation = 1:MC4
        epsilon4 = sqrt(sigma2_e4)*randn(T4+burn4,1);
        y4_all = filter(1,[1 -phi4],epsilon4);
        y4 = y4_all(burn4+1:end);

        % Estimate an intercept and phi with OLS, as in Question 3.
        Y4 = y4(2:end);
        X4 = [ones(T4-1,1) y4(1:end-1)];
        beta_hat4 = X4\Y4;
        phi_hat4_all(simulation,size_index) = beta_hat4(2);
    end

    % Summarize the empirical sampling distribution for this T.
    mean_phi_hat4(size_index) = mean(phi_hat4_all(:,size_index));
    sd_phi_hat4(size_index) = std(phi_hat4_all(:,size_index));
    bias_phi_hat4(size_index) = mean_phi_hat4(size_index)-phi4;
    asymptotic_sd4(size_index) = sqrt((1-phi4^2)/T4);

    % Put the four requested distributions in a simple 2-by-2 layout.
    subplot(2,2,size_index);
    histogram(phi_hat4_all(:,size_index),30,'Normalization','pdf', ...
        'FaceColor',[0.35 0.60 0.85]);
    hold on;
    true_line4 = xline(phi4,'k--','LineWidth',1.1);
    mean_line4 = xline(mean_phi_hat4(size_index),'r:','LineWidth',1.3);
    grid on;
    xlabel('Estimated phi');
    ylabel('Density');
    title(['T = ' num2str(T4) ', mean = ' num2str(mean_phi_hat4(size_index),'%.3f')]);

    % One compact legend is enough because the line meanings are identical
    % in all four panels.
    if size_index==1
        legend([true_line4 mean_line4],{'True phi','Monte Carlo mean'}, ...
            'Location','northwest');
    end
end

sgtitle('Question 4: OLS sampling distribution when phi = 0.9');
exportgraphics(figure4,fullfile(output_dir,'question4_distributions_by_sample_size.png'), ...
    'Resolution',300);
close(figure4);

% Save the four distributions themselves.
question4_draws = array2table(phi_hat4_all, ...
    'VariableNames',{'phi_hat_T50','phi_hat_T100','phi_hat_T200','phi_hat_T1000'});
writetable(question4_draws,fullfile(output_dir,'question4_monte_carlo_draws.csv'));

% Save the changing center, spread, and bias as T increases.
question4_summary = table(sample_sizes4',repmat(MC4,4,1),repmat(burn4,4,1), ...
    mean_phi_hat4,sd_phi_hat4,bias_phi_hat4,asymptotic_sd4, ...
    'VariableNames',{'T','Simulations','BurnIn','MeanPhiHat','StdPhiHat', ...
    'Bias','AsymptoticStd'});
writetable(question4_summary,fullfile(output_dir,'question4_summary.csv'));

fprintf('QUESTION 4\n');
disp(question4_summary);
disp(' ');

%% Question 5: Fit an AR(1) regression to data generated by an MA(1)
% The data-generating process is x_t=epsilon_t+theta*epsilon_(t-1), with
% theta=0.6.  The regression written in PS1 is x_t=a*x_(t-1)+v_t, with no
% intercept.  Its probability limit is the population linear-projection
% coefficient gamma(1)/gamma(0):
%
%             plim(a_hat) = theta/(1+theta^2).

theta5 = 0.6;
sigma2_e5 = 1;
sample_sizes5 = [50 100 200 1000];
MC5 = 1000;
theoretical_limit5 = theta5/(1+theta5^2);

rng(505,'twister');

a_hat5_all = zeros(MC5,length(sample_sizes5));
mean_a_hat5 = zeros(length(sample_sizes5),1);
sd_a_hat5 = zeros(length(sample_sizes5),1);
bias_from_limit5 = zeros(length(sample_sizes5),1);

figure5 = figure('Visible','off','Color','w','Position',[100 100 1100 750]);

for size_index = 1:length(sample_sizes5)
    T5 = sample_sizes5(size_index);

    for simulation = 1:MC5
        % Generate T5+1 innovations.  The extra first draw is epsilon_0, so
        % all T5 observations of x have the stationary MA(1) construction.
        epsilon5 = sqrt(sigma2_e5)*randn(T5+1,1);
        x5 = epsilon5(2:end)+theta5*epsilon5(1:end-1);

        % Estimate the exact regression stated in PS1: x_t on x_(t-1)
        % without an intercept.  X5 is therefore a single column vector.
        Y5 = x5(2:end);
        X5 = x5(1:end-1);
        a_hat5 = X5\Y5;
        a_hat5_all(simulation,size_index) = a_hat5;
    end

    mean_a_hat5(size_index) = mean(a_hat5_all(:,size_index));
    sd_a_hat5(size_index) = std(a_hat5_all(:,size_index));
    bias_from_limit5(size_index) = mean_a_hat5(size_index)-theoretical_limit5;

    subplot(2,2,size_index);
    histogram(a_hat5_all(:,size_index),30,'Normalization','pdf', ...
        'FaceColor',[0.35 0.60 0.85]);
    hold on;
    limit_line5 = xline(theoretical_limit5,'k--','LineWidth',1.1);
    mean_line5 = xline(mean_a_hat5(size_index),'r:','LineWidth',1.3);
    grid on;
    xlabel('Estimated a');
    ylabel('Density');
    title(['T = ' num2str(T5) ', mean = ' num2str(mean_a_hat5(size_index),'%.3f')]);

    if size_index==1
        legend([limit_line5 mean_line5],{'Theoretical limit','Monte Carlo mean'}, ...
            'Location','northwest');
    end
end

sgtitle('Question 5: AR(1) coefficient estimated from MA(1) data');
exportgraphics(figure5,fullfile(output_dir,'question5_distributions_by_sample_size.png'), ...
    'Resolution',300);
close(figure5);

question5_draws = array2table(a_hat5_all, ...
    'VariableNames',{'a_hat_T50','a_hat_T100','a_hat_T200','a_hat_T1000'});
writetable(question5_draws,fullfile(output_dir,'question5_monte_carlo_draws.csv'));

question5_summary = table(sample_sizes5',repmat(MC5,4,1),mean_a_hat5, ...
    sd_a_hat5,bias_from_limit5,repmat(theoretical_limit5,4,1), ...
    'VariableNames',{'T','Simulations','MeanAHat','StdAHat', ...
    'BiasFromTheoreticalLimit','TheoreticalLimit'});
writetable(question5_summary,fullfile(output_dir,'question5_summary.csv'));

% A second figure makes convergence of the Monte Carlo mean especially clear.
figure6 = figure('Visible','off','Color','w','Position',[100 100 850 500]);
plot(sample_sizes5,mean_a_hat5,'bo-','LineWidth',1.5,'MarkerFaceColor','b');
hold on;
yline(theoretical_limit5,'k--','theta/(1+theta^2)','LineWidth',1.2);
grid on;
xlabel('Sample size T');
ylabel('Monte Carlo mean of a-hat');
title('Question 5: Convergence of the mean OLS estimate');

exportgraphics(figure6,fullfile(output_dir,'question5_mean_convergence.png'), ...
    'Resolution',300);
close(figure6);

fprintf('QUESTION 5\n');
fprintf('Theoretical probability limit theta/(1+theta^2): %.9f\n', ...
    theoretical_limit5);
disp(question5_summary);
disp(' ');

%% Question 6: Test the general ARMA(p,q) for-loop function
% The separate file arma_generator_ps1.m contains the requested function.
% Here we demonstrate an ARMA(1,2), save its output, and verify that giving
% lag-polynomial coefficients or the equivalent roots produces the same
% series when MATLAB receives the same random-number seed.

T6 = 500;
burn6 = 500;
sigma2_e6 = 0.4;

% A(L)=1-0.8L and B(L)=1+0.4L+0.2L^2.
AR_polynomial6 = [1 -0.8];
MA_polynomial6 = [1 0.4 0.2];

rng(606,'twister');
[y6_coefficients,epsilon6,AR_used6,MA_used6] = arma_generator_ps1( ...
    T6,sigma2_e6,AR_polynomial6,MA_polynomial6,'coefficients',burn6);

% roots returns the roots of the supplied polynomial.  Reusing seed 606
% makes the shocks identical, so the two parameterizations can be checked
% against one another observation by observation.
AR_roots6 = roots(AR_polynomial6);
MA_roots6 = roots(MA_polynomial6);

rng(606,'twister');
[y6_roots,epsilon6_roots] = arma_generator_ps1( ...
    T6,sigma2_e6,AR_roots6,MA_roots6,'roots',burn6);

max_parameterization_difference6 = max(abs(y6_coefficients-y6_roots));
max_shock_difference6 = max(abs(epsilon6-epsilon6_roots));

% Independent validation against MATLAB's filter command.
% burn=0 is necessary here because filter also starts from zero initial conditions.
% Resetting the seed ensures that the loop and filter receive exactly the same innovations.
validation_seed6 = 707;

rng(validation_seed6,'twister');
[y6_loop_check,e6_loop_check] = arma_generator_ps1( ...
    T6,sigma2_e6,AR_polynomial6,MA_polynomial6,'coefficients',0);

rng(validation_seed6,'twister');
e6_filter_check = sqrt(sigma2_e6)*randn(T6,1);
y6_filter_check = filter( ...
    MA_polynomial6,AR_polynomial6,e6_filter_check);

max_loop_filter_shock_difference6 = ...
    max(abs(e6_loop_check-e6_filter_check));

max_loop_filter_series_difference6 = ...
    max(abs(y6_loop_check-y6_filter_check));

fprintf('QUESTION 6\n');
fprintf('Maximum coefficient-root output difference: %.16g\n', ...
    max_parameterization_difference6);
fprintf('Maximum shock difference under the same seed: %.16g\n', ...
    max_shock_difference6);
fprintf('Maximum loop-filter shock difference: %.16g\n', ...
    max_loop_filter_shock_difference6);
fprintf('Maximum loop-filter series difference: %.16g\n\n', ...
    max_loop_filter_series_difference6);

question6_series = table( ...
    (1:T6)',y6_coefficients,epsilon6,y6_roots, ...
    'VariableNames', ...
    {'t','y_from_coefficients','epsilon','y_from_roots'});

writetable(question6_series, ...
    fullfile(output_dir,'question6_arma_example.csv'));

question6_summary = table( ...
    T6,burn6,sigma2_e6,validation_seed6, ...
    string(mat2str(AR_used6)),string(mat2str(MA_used6)), ...
    max_parameterization_difference6,max_shock_difference6, ...
    max_loop_filter_series_difference6, ...
    max_loop_filter_shock_difference6, ...
    'VariableNames',{ ...
    'T','BurnIn','InnovationVariance','ValidationSeed', ...
    'ARPolynomial','MAPolynomial', ...
    'MaxCoefficientRootSeriesDifference', ...
    'MaxCoefficientRootShockDifference', ...
    'MaxLoopFilterSeriesDifference', ...
    'MaxLoopFilterShockDifference'});

writetable(question6_summary, ...
    fullfile(output_dir,'question6_summary.csv'));

figure7 = figure('Visible','off','Color','w','Position',[100 100 1000 650]);
subplot(2,1,1);
plot(1:T6,y6_coefficients,'b-','LineWidth',1.0);
grid on;
xlabel('t');
ylabel('y_t');
title('Question 6: Example ARMA(1,2) realization');

subplot(2,1,2);
plot(1:T6,epsilon6,'k-','LineWidth',0.9);
grid on;
xlabel('t');
ylabel('epsilon_t');
title('White-noise innovations returned by the function');

exportgraphics(figure7,fullfile(output_dir,'question6_arma_example.png'), ...
    'Resolution',300);
close(figure7);

%% Save one MATLAB data file containing all major numerical results
% CSV files are convenient for tables; the MAT file preserves the numerical
% arrays and MATLAB tables exactly for later analysis.

save(fullfile(output_dir,'PS1_all_results.mat'), ...
    'question1_paths','question1_summary', ...
    'question2_path','question2_summary', ...
    'question3_draws','question3_summary', ...
    'question4_draws','question4_summary', ...
    'question5_draws','question5_summary', ...
    'question6_series','question6_summary');

disp('All Problem Set 1 calculations completed successfully.');
disp('All tables, figures, and the MAT file were saved in the output directory.');

% Stop recording the Command Window only after the success message.
diary off;
