function y = RBFInterp(x, para)
ax = para.nodes;
nx = size(x, 1);
np = size(ax, 1);    % np: the size of data set

xmin = para.xmin;
xmax = para.xmax;
ymin = para.ymin;
ymax = para.ymax;
% normalization
x = 2./(repmat(xmax - xmin, nx, 1)) .* (x - repmat(xmin, nx, 1)) - 1;

% r = sqrt(abs(repmat(diag(ax * ax')', nx, 1) + repmat(diag(x * x'), 1, np) - 2 * (x * ax')));
r = dist(x, ax');
switch para.kernel
    case 'gaussian'
%         Phi = 1 / sqrt(2 * pi) * exp(- 1/2 * r.^2);
        Phi = radbas(sqrt(-log(.5))*r);
    case 'cubic'
        Phi = r.^3;
end

y = Phi * para.alpha + [ones(nx, 1), x] * para.beta;
% renormalization
y = repmat(ymax - ymin, nx, 1)./2 .* (y + 1) + repmat(ymin, nx, 1);

% switch para.kernel
%     case 'gaussian'
% %         sigma2(i, :) = 1 / sqrt(2 * pi) - Phi(i, :)*inv(para.Phi)*(Phi(i, :))';
%         sigma2 = 1 / sqrt(2 * pi) - Phi*(para.Phi\(Phi'));
%     case 'cubic'
% %         for i = 1 : size(y, 1)
% %             sigma2(i, :) = -Phi(i, :)*inv(para.Phi)*(Phi(i, :))';
% %         end
%         sigma2 =  - Phi*(para.Phi\(Phi'));
% end
% sigma2 = abs(diag(sigma2));

end