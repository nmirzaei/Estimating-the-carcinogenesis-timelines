function out = aicbic_msce_poisson_brc(theta, T, Age, yobs, cell_num, N0, bm, Menarche_Age, varargin)
% Poisson likelihood AIC/BIC for MSCE model
%
% Model observable used: y(2,:) at Age

p = inputParser;
p.addParameter('k', numel(theta));
p.parse(varargin{:});
opt = p.Results;

Age  = Age(:);

Tmax = min(T,max(Age));
tt_plot = 0:min(T,Tmax);
x_ic = [1 0 1 0 1 0];

[~,Yage] = ode15s(@(t,x) hazardfunc_multi_malignant_cells_SimpleBirth_brc( ...
    t, x, theta, cell_num, N0, bm,Menarche_Age), tt_plot, x_ic);

lambda = Yage(:,2);  % predicted intensity
IDX1 = find(Age==tt_plot(end));
IDX2 = find(tt_plot==min(Age));
yobs_ = yobs(1:IDX1);
lambda_=lambda(IDX2:end);
% ensure positivity
lambda_ = max(lambda_, 1e-8);

% Poisson log-likelihood
logL = sum(yobs_ .* log(lambda_) - lambda_ - gammaln(yobs_+1));

n = numel(yobs_);
k = opt.k;

AIC = 2*k - 2*logL;
BIC = k*log(n) - 2*logL;

out = struct('logL',logL,'AIC',AIC,'BIC',BIC,'lambda',lambda);
end