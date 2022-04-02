function qos_EJC = qos_EJC(k, T, beta, gamma, K, d, r, b, Max, init)
    qos_EJC = sum(beta(k).*((Max.b(k) - b)/K - r - sum(d)).^2) ...
            - sum(gamma(k).*(b - [init.b.C(k) b(1:T-1)]));
end

