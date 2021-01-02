function ci = nlparci(beta,resid,dummy,Jacobian)

J = Jacobian;
Sigma = [];
alpha = 0.05;


% Remove missing values.
resid = resid(:);
missing = isnan(resid);
if ~isempty(missing)
    resid(missing) = [];
end

n = length(resid);
p = numel(beta);
v = n-p;

J(missing,:) = [];
if size(J,1)~=n || size(J,2)~=p
    error('Jacobian Size Mismatch');
end

% Approximation when a column is zero vector
temp = find(max(abs(J)) == 0);
if ~isempty(temp)
    J(:,temp) = sqrt(eps(class(J)));
end

% Calculate covariance matrix
[dummy,R] = qr(J,0);
Rinv = R\eye(size(R));
diag_info = sum(Rinv.*Rinv,2);

rmse = norm(resid) / sqrt(v);
se = sqrt(diag_info) * rmse;

% Calculate confidence interval
delta = se * tinv(1-alpha/2,v);
ci = [(beta(:) - delta) (beta(:) + delta)];
   