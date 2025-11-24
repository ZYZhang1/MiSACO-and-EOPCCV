function OffDec = ACO_MV_generate(Parent,len_r,len_c,k,l,rank_v,m,up_r,dn_r)
q = 0.05099;
kesi  = 0.6795;
w = ( 1/(q*k*sqrt(2*pi)) )*exp( -((rank_v-1).^2)/(2*(q*k)^2) );
p_ro = w/sum(w);

v_r = Parent(:,1:len_r);
v_c = Parent(:,len_r+1:len_r+len_c);
for i = 1:m
    v_r_generate(i,:) = ACO_RO(v_r,p_ro,len_r,k,kesi);
    v_r_generate(i,:) = Repair(v_r_generate(i,:),up_r,dn_r);
       
    v_c_generate(i,:) = ACO_C(v_c,l,len_c,q,w);
end

OffDec = [v_r_generate,v_c_generate];
end

function [n_x] = ACO_C(x,l,len_c,q,w)
for j = 1:len_c
    [pl] = Cal_pl(x(:,j),l(:,j),w,q);
    idx_gc = Select(pl);            

    n_x(1,j) = idx_gc;
end 
end

function [out] = ACO_RO(x,p,len,k,kesi)
global zzz BBB v000 xxxxx idx_grrr
idx_gr = Select(p);
[z,B,v0] = Rot(x,idx_gr,len);
% z = x;
idx_grrr = idx_gr;
zzz = z;
BBB = B;
v000 = v0;
xxxxx = x;
for j = 1:len
    mu = z(idx_gr,j);
    sigma = kesi*sum( abs( z(idx_gr,j)-z(:,j) ) )/(k-1);
           
    n_z(1,j) = mu + sigma*randn(1);     
end
out = n_z*B'+v0;
end

function [z_r,B,v0] = Rot(v_r,idx_gr,len_r)
if (sum(sum(v_r - v_r(idx_gr,:))) ~= 0)&&(len_r>1)
    B = VCH(v_r,v_r(idx_gr,:));
else
    B = eye(len_r);
end

if rank(B) ~= len_r
%     B = B;
% else
    B = eye(len_r);
end

z_r = (v_r - v_r(idx_gr,:))*B;
% z_r = v_r*B;
v0 = v_r(idx_gr,:);
end

function i = Select(p)
p_sel = cumsum(p);
R     = rand;
i = 1;
while p_sel(i) < R
    i = i+1;
end
end

function out = Repair(x,up,dn)
x = (x>=up).*max(dn,2*up-x) + (x<up).*x;
x = (x<=dn).*min(up,2*dn-x) + (x>dn).*x;
out = x;
end

function [out] = Cal_pl(vc,l,w,q)

for i = 1:l
    idx_l = ( vc==i );
    u(i) = sum( idx_l );

    if isempty(w(idx_l))
       wjl(i) = 0; 
    else
       wjl(i) = max(w(idx_l));
    end
end

eta = 100*sum( u==0 );
for i = 1:l
    if (eta>0)&&(u(i)>0)
        wl(i) = wjl(i)/u(i) + q/eta;
    elseif (eta==0)&&(u(i)>0)
        wl(i) = wjl(i)/u(i);
    elseif (eta>0)&&(u(i)==0)
        wl(i) = q/eta;
    end
end

out = wl/sum(wl);
end

function out = VCH(s,sl)
[~,n] = size(s);
for i = 1:n
    ds = sqrt( sum( (sl(:,i:n) - s(:,i:n)).^2, 2) );
    p  = ds.^4 / sum(ds.^4);
    idx = Select(p);
    A(i,:) = sl - s(idx,:);
    s(idx,:) = [];
end
if max(max(A))<1e-5
    B = Gram_Schmidt_process(rand(n));
else
    B = Gram_Schmidt_process(A');
end
out = B;
end
