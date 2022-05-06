function [s,nu] = lmo(z,atom_type)
    if strcmp(atom_type,'1norm')
    [~,sk] = max(abs(z));
    s = z*0;
    s(sk) = sign(z(sk));
    nu = abs(z(sk));
    end
end