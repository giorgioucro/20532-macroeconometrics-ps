%% 20532 Macroeconometrics - Problem Set 2, Fall 2026
% Requirements: MATLAB (R2021a or newer) and Statistics and Machine Learning
% Toolbox (prctile, norminv, chi2inv, finv, tinv, fcdf, chi2cdf).
clear; clc; close all; format compact;
output_dir = fileparts(mfilename('fullpath'));
log_file = fullfile(output_dir,'PS2_command_window_log.txt');
if isfile(log_file), delete(log_file); end
diary(log_file);
fprintf('20532 Macroeconometrics - PS2, Fall 2026\nMATLAB %s\n',version);
fprintf('Outputs: %s\n',output_dir);
% Confidence level 95%
alpha = 0.05;

%% Question 1(a-d): pure random walk and unit-root inference
% There are T=250 observed levels y_1,...,y_T, initialized at y_0=0.
% No burn-in is used: a random walk has no stationary initial distribution.
% All four parts use the SAME driftless DGP, with innovation variance 1.
rng(2053201,'twister');
T1 = 250; N1 = 10000; sigma2_1 = 1; n1 = T1-1;
phi_no_constant = zeros(N1,1);
phi_with_constant = zeros(N1,1);
df_coefficients = zeros(N1,2); % columns: intercept alpha and rho=phi-1
df_t = zeros(N1,1); joint_F = zeros(N1,1); joint_Wald = zeros(N1,1);
for replication = 1:N1
    y = cumsum(sqrt(sigma2_1)*randn(T1,1));
    y_lag = y(1:end-1); dy = diff(y);

    % (a) OLS in the specified no-intercept levels regression.
    % Also retain the intercept regression slope to connect (a) with (b).
    % These are different estimators.
    phi_no_constant(replication) = y_lag\y(2:end);
    X = [ones(n1,1), y_lag];
    b = X\dy;
    residual = dy-X*b;
    SSR_u = residual'*residual;
    sigma2_hat = SSR_u/(n1-2);
    covariance_b = sigma2_hat*((X'*X)\eye(2));
    df_coefficients(replication,:) = b';
    phi_with_constant(replication) = 1+b(2);

    % (b) DF t statistic for rho=0 in dy=alpha+rho*y_lag+error.
    % The intercept is estimated even though its DGP value is zero.
    df_t(replication) = b(2)/sqrt(covariance_b(2,2));

    % (d) joint null: alpha=0 AND rho=0 (q=2 restrictions).
    % Under this null fitted dy is zero, not its sample mean. Hence the
    % restricted SSR is sum(dy.^2). W=q*F; chi-square criticals apply to W.
    SSR_r = dy'*dy;
    joint_F(replication) = ((SSR_r-SSR_u)/2)/sigma2_hat;
    joint_Wald(replication) = b'*(covariance_b\b);
end
assert(max(abs(joint_Wald-2*joint_F))<1e-8,'Q1: W must equal 2F.');

estimator = ["OLS slope, no intercept"; "OLS slope, intercept"];
phi_draws = [phi_no_constant,phi_with_constant];
q1_estimators = table(estimator,mean(phi_draws)',std(phi_draws)', ...
    median(phi_draws)',prctile(phi_draws,2.5)',prctile(phi_draws,97.5)', ...
    'VariableNames',{'Estimator','Mean','SD','Median','P025','P975'});
disp('Q1(a): estimator distributions'); disp(q1_estimators);
writetable(q1_estimators,fullfile(output_dir,'question1_estimators.csv'));

% (c) Hamilton, Table B.6, p. 763, Case 2, T=250: intercept, no trend.
% Compare BOTH tails as well as the usual lower-tail significance levels.
percentile = [1;2.5;5;10;90;95;97.5;99];
tabulated_t = [-3.46;-3.14;-2.88;-2.57;-0.42;-0.06;0.24;0.62];
empirical_t = prctile(df_t,percentile); empirical_t = empirical_t(:);
q1_percentiles = table(percentile,empirical_t,tabulated_t, ...
    empirical_t-tabulated_t,'VariableNames', ...
    {'Percentile','Empirical_t','DF_intercept_T250','Empirical_minus_table'});
disp('Q1(c): Dickey-Fuller t percentiles'); disp(q1_percentiles);
writetable(q1_percentiles,fullfile(output_dir,'question1_df_percentiles.csv'));

% A one-sided 5% Normal test uses -1.645, not the two-sided value -1.96.
% The non-standard joint F is the DF Phi_1 statistic (Enders, p. 207).
% Hamilton Table B.7, p. 764, Case 2, T=250 gives 4.63 at 5%.
normal_critical = norminv(alpha);
chi2_critical = chi2inv(1-alpha,2);
usual_F_critical = finv(1-alpha,2,n1-2);
phi1_critical = 4.63;
test_name = ["Q1b: t, standard Normal";"Q1c: t, Dickey-Fuller"; ...
    "Q1d: W=2F, chi-square(2)";"Q1d: F, conventional F(2,247)"; ...
    "Q1d: F, Dickey-Fuller Phi_1"];
critical_value = [normal_critical;-2.88;chi2_critical;usual_F_critical;phi1_critical];
rejection_rate = [mean(df_t<normal_critical);mean(df_t<-2.88); ...
    mean(joint_Wald>chi2_critical);mean(joint_F>usual_F_critical); ...
    mean(joint_F>phi1_critical)];
mc_standard_error = sqrt(rejection_rate.*(1-rejection_rate)/N1);
q1_rejections = table(test_name,critical_value,rejection_rate,mc_standard_error, ...
    'VariableNames',{'Test','Critical_value','Rejection_frequency','MC_SE'});
disp('Q1(b,d): rejection frequencies under the true joint null'); disp(q1_rejections);
writetable(q1_rejections,fullfile(output_dir,'question1_rejections.csv'));
% An additional check of the F distribution, using published critical values
% rather than estimating a cutoff from the same draws being evaluated.
upper_tail_probability = [0.10;0.05;0.01];
tabulated_phi1 = [3.81;4.63;6.52];
empirical_phi1 = prctile(joint_F,100*(1-upper_tail_probability));
q1_joint_percentiles = table(upper_tail_probability,empirical_phi1(:), ...
    tabulated_phi1,'VariableNames',{'Upper_tail_probability','Empirical_F','DF_Phi1_T250'});
writetable(q1_joint_percentiles,fullfile(output_dir,'question1_joint_percentiles.csv'));
centered_t = df_t-mean(df_t);
t_skewness = mean(centered_t.^3)/mean(centered_t.^2)^1.5;
fprintf('Q1(b): t mean %.4f, median %.4f, skewness %.4f.\n', ...
    mean(df_t),median(df_t),t_skewness);
fprintf(['The DF t is shifted left and is not exactly symmetric; compare its\n' ...
    'tails and skewness. Normal inference substantially over-rejects.\n' ...
    'The joint test rejects a composite null; rejection alone need not\n' ...
    'establish stationarity (a drift can also violate the joint null).\n']);

fig = figure('Color','w','Position',[80 80 1250 400]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; hold on;
histogram(phi_no_constant,50,'Normalization','pdf','DisplayName','No intercept');
histogram(phi_with_constant,50,'Normalization','pdf','DisplayName','Intercept');
xline(1,'k--','HandleVisibility','off'); grid on;
xlabel('OLS slope'); ylabel('Density'); title('Q1(a): pure random walk');
legend('Location','northwest');
nexttile; hold on;
histogram(df_t,50,'Normalization','pdf','DisplayName','Empirical DF t');
grid_t = linspace(-5,3,300);
plot(grid_t,exp(-grid_t.^2/2)/sqrt(2*pi),'k-','DisplayName','N(0,1)');
xline(normal_critical,'r--','DisplayName','Normal 5%');
xline(-2.88,'b--','DisplayName','DF 5%');
xlabel('t statistic'); title('Q1(b,c): intercept DF test'); grid on;
legend('Location','northeast');
nexttile; hold on;
histogram(joint_F,50,'Normalization','pdf','DisplayName','Empirical Phi_1');
grid_f = linspace(0,max(joint_F),300);
% chi-square(2)/2 has density exp(-f).
plot(grid_f,exp(-grid_f),'k-','DisplayName','Conventional chi-square(2)/2');
xline(chi2_critical/2,'r--','DisplayName','Conventional 5%');
xline(phi1_critical,'b--','DisplayName','Phi_1 5%');
xlabel('Joint F statistic'); title('Q1(d): two restrictions'); grid on;
legend('Location','northeast');
exportgraphics(fig,fullfile(output_dir,'question1_unit_root_distributions.png'),'Resolution',180);
close(fig);
save(fullfile(output_dir,'question1_simulations.mat'),'T1','N1','sigma2_1', ...
    'phi_no_constant','phi_with_constant','df_coefficients','df_t', ...
    'joint_F','joint_Wald','t_skewness','q1_estimators','q1_percentiles', ...
    'q1_rejections','q1_joint_percentiles');

%% Question 2: the four cases in Enders, pp. 195-199
% In every case estimate y_t=a+b*z_t+e_t with conventional OLS standard
% errors. Vary T to distinguish convergence from a finite-sample accident.
% Cases 1-3 have unrelated series and null slope zero. In case 4 the true
% cointegrating slope is ONE: testing b=0 measures power, not test size.
% The conventional t around b=1 is stored separately; cointegration alone
% does not guarantee valid conventional t inference.
rng(2053202,'twister');
sample_sizes2 = [100 500 1000]; N2 = 5000; burn2 = 500;
case_names = ["1: both I(0)";"2: y I(1), z I(0)"; ...
    "3: independent I(1)";"4: cointegrated I(1)"];
true_slope2 = [0 0 0 1];
intercept2 = zeros(N2,4,numel(sample_sizes2));
slope2 = intercept2; t_zero2 = intercept2; t_true2 = intercept2; r_squared2 = intercept2;
q2_rows = cell(4*numel(sample_sizes2),14); row = 0;
for sample_index = 1:numel(sample_sizes2)
    T2 = sample_sizes2(sample_index);
    for case_index = 1:4
        for replication = 1:N2
            switch case_index
                case 1
                    % White noise is I(0). Independent white noise ensures
                    % a clean benchmark with valid classical t inference.
                    y = randn(T2,1); z = randn(T2,1);
                case 2
                    % Enders' orientation: nonstationary dependent variable,
                    % stationary regressor. Only the AR needs a burn-in.
                    y = cumsum(randn(T2,1));
                    stationary_z = filter(1,[1 -0.7],randn(T2+burn2,1));
                    z = stationary_z(burn2+1:end);
                case 3
                    % Independent random walks: any apparent relationship
                    % is spurious; neither series shares the other's shocks.
                    y = cumsum(randn(T2,1)); z = cumsum(randn(T2,1));
                case 4
                    % Enders' common random walk plus independent noise.
                    % y-z is stationary although y and z individually are I(1).
                    common_trend = cumsum(randn(T2,1));
                    y = common_trend+randn(T2,1); z = common_trend+randn(T2,1);
            end
            X = [ones(T2,1),z]; b = X\y;
            residual = y-X*b;
            SSR = residual'*residual;
            covariance_b = (SSR/(T2-2))*((X'*X)\eye(2));
            slope_se = sqrt(covariance_b(2,2));
            intercept2(replication,case_index,sample_index) = b(1);
            slope2(replication,case_index,sample_index) = b(2);
            t_zero2(replication,case_index,sample_index) = b(2)/slope_se;
            t_true2(replication,case_index,sample_index) = (b(2)-true_slope2(case_index))/slope_se;
            r_squared2(replication,case_index,sample_index) = 1-SSR/sum((y-mean(y)).^2);
        end
        slope_draw = slope2(:,case_index,sample_index);
        t0_draw = t_zero2(:,case_index,sample_index);
        ttrue_draw = t_true2(:,case_index,sample_index);
        r2_draw = r_squared2(:,case_index,sample_index);
        critical_t = tinv(1-alpha/2,T2-2);
        row = row+1;
        q2_rows(row,:) = {case_names(case_index),T2,N2,true_slope2(case_index), ...
            mean(slope_draw),std(slope_draw),prctile(slope_draw,2.5), ...
            prctile(slope_draw,97.5),std(t0_draw),mean(abs(t0_draw)>critical_t), ...
            mean(abs(ttrue_draw)>critical_t),mean(r2_draw), ...
            prctile(r2_draw,2.5),prctile(r2_draw,97.5)};
    end
end
q2_summary = cell2table(q2_rows,'VariableNames',{'Case','T','Replications', ...
    'True_slope','Mean_slope','SD_slope','Slope_P025','Slope_P975', ...
    'SD_t_zero','Reject_slope_zero','Reject_true_slope','Mean_R2','R2_P025','R2_P975'});
disp('Q2: four Enders cases; rejection of zero in case 4 is power'); disp(q2_summary);
writetable(q2_summary,fullfile(output_dir,'question2_spurious_summary.csv'));
fprintf(['Q2: independent I(1) series can produce large t statistics and R2,\n' ...
    'while the slope fails to concentrate at zero as T grows. Mixed orders\n' ...
    'can give misleading t tests even with small R2. The cointegrated case\n' ...
    'has a meaningful long-run slope, but usual t inference can still fail.\n']);

% One compact figure compares the requested three objects in all four cases.
% The middle column is centered on the TRUE slope (one in case 4).
plot_sample_index = find(sample_sizes2==500,1);
fig = figure('Color','w','Position',[80 40 1200 1050]);
tiledlayout(4,3,'Padding','compact','TileSpacing','compact');
for case_index = 1:4
    nexttile;
    histogram(slope2(:,case_index,plot_sample_index),50,'Normalization','pdf');
    xline(true_slope2(case_index),'r--'); grid on;
    title(case_names(case_index),'Interpreter','none'); xlabel('OLS slope'); ylabel('Density');
    nexttile;
    histogram(t_true2(:,case_index,plot_sample_index),50,'Normalization','pdf');
    xline([-tinv(.975,498),tinv(.975,498)],'r--'); grid on;
    xlabel('t statistic around true slope');
    nexttile;
    histogram(r_squared2(:,case_index,plot_sample_index),50,'Normalization','pdf');
    if case_index<=2
        xlim([0 1.05*max(r_squared2(:,case_index,plot_sample_index))]);
    else
        xlim([0 1]);
    end
    grid on; xlabel('R^2');
end
sgtitle(sprintf('Q2: T=500, %d replications per case',N2));
exportgraphics(fig,fullfile(output_dir,'question2_four_cases.png'),'Resolution',180);
close(fig);
save(fullfile(output_dir,'question2_simulations.mat'),'sample_sizes2','N2','burn2', ...
    'case_names','true_slope2','intercept2','slope2','t_zero2','t_true2', ...
    'r_squared2','q2_summary');

%% Question 3: non-fundamental DGP, VAR(4) and Cholesky IRFs
% Exact DGP in the PS:
% x_t = eta_t + epsilon_(t-2)
% y_t = beta/(1-beta)*eta_t + beta^2/(1-beta)*epsilon_t + beta*epsilon_(t-1).
% Shocks are independent, with variances 1 and 0.8. Generate two genuine
% presample shocks, then retain exactly 500 observations (no MA burn-in).
rng(2053203,'twister');
beta3 = 0.6; shock_covariance3 = diag([1 0.8]);
T3 = 500; N3 = 1000; lags3 = 4; horizon3 = 20;
shock_scale3 = diag(sqrt(diag(shock_covariance3)));
true_irf3 = zeros(2,2,horizon3+1);
true_irf3(:,:,1) = [1 0; beta3/(1-beta3) beta3^2/(1-beta3)]*shock_scale3;
true_irf3(:,:,2) = [0 0;0 beta3]*shock_scale3;
true_irf3(:,:,3) = [0 1;0 0]*shock_scale3;
% Every later true response is zero. Both estimated and true responses are
% to ONE STANDARD DEVIATION shocks; no ex-post rescaling masks disagreement.
irf_draws3 = zeros(2,2,horizon3+1,N3);
for replication = 1:N3
    eta = randn(T3+2,1);
    epsilon = sqrt(0.8)*randn(T3+2,1);
    x = eta(3:end)+epsilon(1:end-2);
    y = beta3/(1-beta3)*eta(3:end)+beta3^2/(1-beta3)*epsilon(3:end) ...
        +beta3*epsilon(2:end-1);
    [B3,U3,Sigma3,X3,Y3] = estimate_var([x y],lags3);
    P3 = chol(Sigma3,'lower');
    A3 = [B3(2:end,:)';eye(2*(lags3-1)),zeros(2*(lags3-1),2)];

    % The hint is equivalently implemented by two recursive equations.
    % Verify their equivalence to reduced-form OLS + Cholesky once.
    if replication==1
        recursive_x = X3\Y3(:,1);
        recursive_y = [Y3(:,1),X3]\Y3(:,2);
        a2 = recursive_y(1);
        contemporaneous3 = [1 0;-a2 1];
        C3 = contemporaneous3\eye(2);
        recursive_residual = [Y3(:,1)-X3*recursive_x, ...
            Y3(:,2)-[Y3(:,1),X3]*recursive_y];
        % Use the SAME reduced-form degrees of freedom for both covariance
        % representations; a different df convention would change scale.
        recursive_variances = sum(recursive_residual.^2,1)/(size(X3,1)-size(X3,2));
        recursive_impact3 = C3*diag(sqrt(recursive_variances));
        recursive_lags = [recursive_x(2:end)';recursive_y(3:end)'];
        recursive_reduced_lags = C3*recursive_lags;
        recursive_equivalence_error3 = max([ ...
            norm(recursive_impact3-P3,'fro'), ...
            norm(recursive_reduced_lags-B3(2:end,:)','fro')]);
        assert(recursive_equivalence_error3<1e-10,'Q3 recursive equivalence failed.');
    end

    % Let J=[I_2,0,...,0]. Structural IRFs are J*A^h*J'*P.
    % Thus the shock-loading matrix D has P in its top block. Using I_2
    % there gives reduced-form responses unless P is applied separately.
    % D is (8 x 2), and A^h*D is (8 x 2), not (8 x 8).
    response_state = [P3;zeros(2*(lags3-1),2)];
    for h = 0:horizon3
        irf_draws3(:,:,h+1,replication) = response_state(1:2,:);
        response_state = A3*response_state;
    end
end
mean_irf3 = mean(irf_draws3,4);
lower_irf3 = prctile(irf_draws3,2.5,4);
upper_irf3 = prctile(irf_draws3,97.5,4);
% det C(z) = -beta/(1-beta)*(z-1)*(z+beta).
% The root -0.6 INSIDE the unit circle proves non-fundamentalness. The root
% 1 on the boundary additionally rules out an absolutely summable inverse;
% a boundary root alone need not imply non-fundamentalness (Sala, p. 23).
determinant_coefficients3 = [-beta3/(1-beta3),beta3,beta3^2/(1-beta3)];
determinant_roots3 = roots(determinant_coefficients3);
fprintf('Q3: det C(z) roots: %.6f, %.6f\n',determinant_roots3);
fprintf('Recursive-equation / Cholesky equivalence error: %.3g\n',recursive_equivalence_error3);
fprintf(['Q3: triangular impact restrictions are compatible with the DGP but\n' ...
    'cannot recover non-fundamental shocks from present and past data.\n' ...
    'A finite VAR also approximates the dynamics. Increasing T alone does\n' ...
    'not solve structural identification. Bands are pointwise sampling\n' ...
    'percentiles over independent datasets, not confidence bands for a mean.\n']);

variable_names3 = ["x" "y"]; shock_names3 = ["eta" "epsilon"];
q3_rows = cell(4*(horizon3+1),8); row = 0;
horizons3 = 0:horizon3;
fig = figure('Color','w','Position',[80 80 1100 720]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for variable = 1:2
    for shock = 1:2
        truth = reshape(true_irf3(variable,shock,:),1,[]);
        mu = reshape(mean_irf3(variable,shock,:),1,[]);
        lo = reshape(lower_irf3(variable,shock,:),1,[]);
        hi = reshape(upper_irf3(variable,shock,:),1,[]);
        nexttile; hold on;
        fill([horizons3 fliplr(horizons3)],[lo fliplr(hi)],[.82 .88 .97], ...
            'EdgeColor','none','DisplayName','2.5-97.5% sampling percentiles');
        plot(horizons3,mu,'b-','LineWidth',1.7,'DisplayName','Estimated mean');
        plot(horizons3,truth,'k--','LineWidth',1.7,'DisplayName','Analytical truth');
        yline(0,':','HandleVisibility','off'); grid on;
        title(variable_names3(variable)+" response to "+shock_names3(shock));
        xlabel('Horizon (0 = impact)'); ylabel('Response to one-SD shock');
        if variable==1 && shock==1, legend('Location','best'); end
        for h = 0:horizon3
            row = row+1;
            q3_rows(row,:) = {variable_names3(variable),shock_names3(shock),h, ...
                truth(h+1),mu(h+1),lo(h+1),hi(h+1),mu(h+1)-truth(h+1)};
        end
    end
end
sgtitle(sprintf('Q3: VAR(4), T=%d, N=%d; Cholesky ordering [x,y]',T3,N3));
exportgraphics(fig,fullfile(output_dir,'question3_irf_comparison.png'),'Resolution',180);
close(fig);
q3_irfs = cell2table(q3_rows,'VariableNames',{'Response','Shock','Horizon', ...
    'True_IRF','Mean_estimated_IRF','P025','P975','Mean_minus_true'});
writetable(q3_irfs,fullfile(output_dir,'question3_irf_summary.csv'));
disp('Q3: impact and first two horizons'); disp(q3_irfs(q3_irfs.Horizon<=2,:));
save(fullfile(output_dir,'question3_simulations.mat'),'beta3','shock_covariance3', ...
    'T3','N3','lags3','horizon3','true_irf3','irf_draws3','mean_irf3', ...
    'lower_irf3','upper_irf3','determinant_coefficients3','determinant_roots3', ...
    'contemporaneous3','C3','recursive_impact3','recursive_equivalence_error3','q3_irfs');

%% Question 4: Romer-Romer shocks and conditional Granger causality
data_file4 = fullfile(output_dir,'Romer_Romer.xlsx');
if ~isfile(data_file4)
    ancestor = output_dir;
    while ~isfile(data_file4)
        candidate = fullfile(ancestor,'MATLAB','PS2','Romer_Romer.xlsx');
        if isfile(candidate), data_file4 = candidate; break; end
        parent = fileparts(ancestor);
        if strcmp(parent,ancestor), break; end
        ancestor = parent;
    end
end
assert(isfile(data_file4),'Place Romer_Romer.xlsx alongside this script.');
data_table4 = readtable(data_file4,'Sheet',1,'VariableNamingRule','preserve');
names4 = string(data_table4.Properties.VariableNames(2:5));
assert(isequal(names4,["inflation","unemployment","ffr","Romer Romer quarterly"]), ...
    'Unexpected column headers: verify the current course workbook.');
series4 = data_table4{:,2:5};
% Excel formatting extends past the actual sample. Remove only wholly
% blank rows; never silently splice together observations across gaps.
keep4 = ~all(ismissing(series4),2);
series4 = series4(keep4,:); dates4 = data_table4{keep4,1};
if isnumeric(dates4), dates4 = datetime(dates4,'ConvertFrom','excel'); end
if ~isdatetime(dates4), dates4 = datetime(dates4); end
assert(all(isfinite(series4),'all') && all(~isnat(dates4)), 'Missing observations in Q4.');
quarter_index4 = year(dates4)*4+ceil(month(dates4)/3)-1;
assert(all(diff(quarter_index4)==1),'Quarterly sample must be contiguous and ordered.');
quarters4 = string(year(dates4))+"Q"+string(ceil(month(dates4)/3));
assert(size(series4,1)==112 && quarters4(1)=="1969Q1" && quarters4(end)=="1996Q4", ...
    'The bundled current dataset should span 1969Q1-1996Q4 (112 quarters).');
lags4 = 4;
[B4,U4,Sigma4,X4,Y4] = estimate_var(series4,lags4);
effective_n4 = size(X4,1); parameters4 = size(X4,2);
df4 = effective_n4-parameters4; inverse_XX4 = (X4'*X4)\eye(parameters4);
fprintf('Q4: original variable names: %s\n',strjoin(names4,', '));
fprintf('Data: %s through %s (%d quarters); estimation: %s through %s.\n', ...
    quarters4(1),quarters4(end),size(series4,1),quarters4(lags4+1),quarters4(end));
fprintf('VAR(4) + intercept: n=%d, parameters/equation=%d, residual df=%d.\n', ...
    effective_n4,parameters4,df4);

% All tests condition on the SAME four-variable VAR. The columns of X are
% [1, inflation_l1, unemployment_l1, ffr_l1, RR_l1, inflation_l2,...].
% Six individual-direction tests plus a joint test of all macro lags in
% the RR equation address whether macro information predicts the shocks.
causes4 = {4,4,4,1,2,3,1:3}; effects4 = [1 2 3 4 4 4 4];
q4_rows = cell(7,12);
for test = 1:7
    causes = causes4{test}; equation = effects4(test);
    restrictions = [];
    for lag = 1:lags4
        restrictions = [restrictions,1+(lag-1)*4+causes]; %#ok<AGROW>
    end
    q = numel(restrictions);
    restricted_coefficients = B4(restrictions,equation);
    covariance_restrictions = Sigma4(equation,equation)*inverse_XX4(restrictions,restrictions);
    W = restricted_coefficients'*(covariance_restrictions\restricted_coefficients);
    F = W/q;
    p_F = fcdf(F,q,df4,'upper'); p_chi2 = chi2cdf(W,q,'upper');

    % Independent restricted/unrestricted SSR check catches lag-index or
    % degrees-of-freedom errors in the Wald implementation.
    restricted_X = X4; restricted_X(:,restrictions) = [];
    restricted_residual = Y4(:,equation)-restricted_X*(restricted_X\Y4(:,equation));
    SSR_r = restricted_residual'*restricted_residual;
    SSR_u = U4(:,equation)'*U4(:,equation);
    F_from_SSR = ((SSR_r-SSR_u)/q)/(SSR_u/df4);
    assert(abs(F-F_from_SSR)<1e-8,'Q4: Wald and SSR F statistics disagree.');
    q4_rows(test,:) = {strjoin(names4(causes),' + '),names4(equation),q,df4,W, ...
        p_chi2,F,p_F,p_F<.01,p_F<.05,p_F<.10,p_chi2<.05};
end
q4_tests = cell2table(q4_rows,'VariableNames',{'Cause','Effect','Restrictions', ...
    'Residual_df','Wald_chi2','P_chi2','F','P_F','Reject_F_1pct', ...
    'Reject_F_5pct','Reject_F_10pct','Reject_chi2_5pct'});
disp('Q4: H0 = all listed cause lags are zero in the effect equation'); disp(q4_tests);
writetable(q4_tests,fullfile(output_dir,'question4_granger_tests.csv'));
fprintf(['Q4: Granger causality measures incremental predictive content, not\n' ...
    'structural causation. Failure to reject macro -> RR is consistent with\n' ...
    'the intended forecast purging, but does not prove shock exogeneity.\n' ...
    'Conventional F and chi-square references are the course approximation\n' ...
    'for this dynamic VAR, not an exact finite-sample distribution.\n']);
save(fullfile(output_dir,'question4_var_results.mat'),'names4','series4','dates4', ...
    'quarters4','lags4','effective_n4','parameters4','df4','B4','U4','Sigma4','q4_tests');

fprintf('\nAll Fall 2026 PS2 questions completed. Outputs saved beside this script.\n');
diary off;

%% Local function: transparent equation-by-equation VAR OLS
function [B,U,Sigma,X,Y] = estimate_var(data,p)
% B has one column per equation. Each lag block contains all variables in
% their original order. Regressors are identical across equations, so X\Y
% is exactly the same as estimating the equations separately (Sala's style).
    [T,n] = size(data);
    X = ones(T-p,1+n*p);
    for lag = 1:p
        X(:,2+(lag-1)*n:1+lag*n) = data(p+1-lag:T-lag,:);
    end
    Y = data(p+1:end,:);
    B = X\Y;
    U = Y-X*B;
    Sigma = (U'*U)/(size(X,1)-size(X,2));
end
