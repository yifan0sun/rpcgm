function [s,nu] = lmo(z,atom_type,groups)
    if strcmp(atom_type,'groupnorm')
        nz = zeros(length(groups),1);
        for k = 1:length(groups)
            nz(k) = norm(z(groups{k}));
        end
    [nu,sk] = max(nz);
    s = z*0;
    s(groups{sk}) = z(groups{sk});
    s = s / norm(s);
    end
end