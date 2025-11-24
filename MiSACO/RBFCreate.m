function para = RBFCreate(ax, ay, kernel)
[N, D] = size(ax);
% Normlization
xmin = min(ax, [], 1);
xmax = max(ax, [], 1);
ymin = min(ay, [], 1);
ymax = max(ay, [], 1);
ax = 2./(repmat(xmax - xmin, N, 1)) .* (ax - repmat(xmin, N, 1)) - 1;
ay = 2./(repmat(ymax - ymin, N, 1)) .* (ay - repmat(ymin, N, 1)) - 1;

% xxh = ax * ax';
% r = sqrt(repmat(diag(xxh)', N, 1) + repmat(diag(xxh), 1, N) - 2 * xxh);
r = dist(ax, ax');
switch kernel
    case 'gaussian'
%         Phi = 1 / sqrt(2 * pi) * exp(- 1/2 * r.^2);
        Phi = radbas(sqrt(-log(.5))*r);
    case 'cubic'
        Phi = r.^3;
end
P = [ones(N, 1), ax];
A = [Phi, P; P' zeros(D + 1, D + 1)];
b = [ay; zeros(D + 1, size(ay, 2))];

% if rcond(A) <= eps
%     para.alpha = Phi \ ay;
% else
    theta = A \ b;
    
    para.alpha = theta(1 : N, :);
    para.beta = theta(N + 1 : end, :);
% end
para.xmin = xmin;
para.xmax = xmax;
para.ymin = ymin;
para.ymax = ymax;
para.nodes = ax;    % normlization
para.kernel = kernel;
para.Phi = Phi;
para.ax = ax;
para.ay = ay;
end