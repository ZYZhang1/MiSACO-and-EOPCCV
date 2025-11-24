function [v_r, v_o, v_c, f, rank_out] = Initialize(Problem)


up_r = Par.up_r;
dn_r = Par.dn_r;
up_o = Par.up_o;
dn_o = Par.dn_o;
l = Par.l;
k = Par.k;
len_r = Par.len_r;
len_o = Par.len_o;
len_c = Par.len_c;
func_num = Par.func_num;
shift = Par.shift;

v_r = repmat((up_r - dn_r), k, 1) .* lhsdesign(k, len_r) + repmat(dn_r, k, 1);
v_o = repmat((up_o - dn_o), k, 1) .* lhsdesign(k, len_o) + repmat(dn_o, k, 1);
v_c = ceil(repmat(l, k, 1) .* lhsdesign(k, len_c));

[f] = Fitness_func(v_r, v_o, v_c,func_num,shift);

[~,rank_v] = sort(f);
[~,rank_out] = sort(rank_v);