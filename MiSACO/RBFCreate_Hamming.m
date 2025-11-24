function para = RBFCreate_Hamming(ax, ay, len_c, up, dn, kernel)
[N, D] = size(ax);

ax_nc = ax(:,1:D-len_c);
ax_nc = (ax_nc - dn)./(up - dn) - 0.5; 
ax_c  = ax(:,D-len_c+1:D);

ymin = min(ay, [], 1);
ymax = max(ay, [], 1);
ay = 2./(repmat(ymax - ymin, N, 1)) .* (ay - repmat(ymin, N, 1)) - 1;

ax = [ax_nc,ax_c];

for i = 1:N 
    dis_gower(i,:) = sum( ax_c~= ax_c(i,:) ,2)';
    dis_euc(i,:) = sum( abs(ax_nc - ax_nc(i,:)),2)';
end

r = (dis_gower + dis_euc);
switch kernel
    case 'gaussian'
        Phi = radbas(sqrt(-log(.5))*r);
    case 'cubic'
        Phi = r.^3;
end
theta = (Phi'*Phi)\(Phi'*ay);
para.alpha = theta(1 : N, :);
para.ymin = ymin;
para.ymax = ymax;
para.nodes = ax;
para.kernel = kernel;
para.Phi = Phi;
end