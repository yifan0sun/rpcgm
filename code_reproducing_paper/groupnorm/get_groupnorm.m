
function gx = get_groupnorm(x,groups)
n = length(x);
if isnan(sum(x)) || isinf(sum(x))
    gx = nan;
    return
end
cvx_begin
for k = 1:length(groups)
    eval(sprintf('variable s%d(length(groups{%d}))',k,k))

    if k == 1
        t = norm(s1);
        bigx(groups{1},1) = s1;
        
bigx(n,1) = 0;
    else
        eval(sprintf('sk = s%d;',k));
        t = t + norm(sk);
        bigx(groups{k}) = bigx(groups{k}) + sk;
    end
end

minimize t
subject to
    bigx == x;

cvx_end
gx = cvx_optval;

end