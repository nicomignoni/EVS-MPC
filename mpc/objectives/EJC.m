function EJC = EJC(k, d, r, c, p, l)
    EJC = sum(p.C(k,:).*(sum(squeeze(l.psi(:,k,:)).*d)) + r.*c.sell, "all");
end

