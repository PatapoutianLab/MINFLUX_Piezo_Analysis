
% Code to analyze the variance and 95% confidence interval of inter-blade distance data

close all
clear
clc

% Input inter-blade distance data into two matrices

A = [];
B = [];

%%

% Calculate variances and degrees of freedom
var_A = var(A);
var_B = var(B);
df_A = length(A) - 1;
df_B = length(B) - 1;

% Calculate confidence intervals for variances
alpha = 0.05; % significance level
CI_A = df_A * var_A ./ chi2inv([1-alpha/2, alpha/2], df_A);
CI_B = df_B * var_B ./ chi2inv([1-alpha/2, alpha/2], df_B);

% Calculate F-statistic and p-value
F = var_A / var_B;
pval = fcdf(F, df_A, df_B, 'upper');

% Plot error bars for confidence intervals
figure;
errorbar(1, var_A, var_A - CI_A(1), CI_A(2) - var_A, 'bo', 'LineWidth', 2);
hold on;
errorbar(2, var_B, var_B - CI_B(1), CI_B(2) - var_B, 'ro', 'LineWidth', 2);
xlim([0.5, 2.5]);
xticks([1, 2]);
xticklabels({'Dataset A', 'Dataset B'});
ylabel('Variance');
title('95% Confidence intervals for variances');

% Add significance asterisk if CIs don't overlap
if ~(CI_A(2) < CI_B(1) || CI_B(2) < CI_A(1))
text(1.5, max(var_A, var_B), '*', 'FontSize', 24, 'HorizontalAlignment', 'center');
end

% Print F-statistic and p-value
fprintf('F-statistic: %f\n', F);
fprintf('p-value: %f\n', pval);

