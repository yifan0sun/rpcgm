function [s,nu] = lmo(z,atom_type)
    if strcmp(atom_type,'tv')
        z = z - mean(z);
        u = cumsum(z);
        u = u(1:end-1);
        [nu,sk] = max(abs(u));
        s = z*0;
        s(1:sk) = sign(u(sk));
        s(sk+1:end) = -sign(u(sk));
    
    
    end
end