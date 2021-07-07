function Alm = AssociatedLegendre(l,m,x)
Alm = 0;
    for r = 0 : floor(1/2 * l - 1/2 * abs(m))
        Alm = Alm + (-1)^r * nchoosek(l - 2 * r, abs(m)) * nchoosek(l, r) * ...
            nchoosek(2 * l - 2 * r, l) * x.^(l - 2 * r - abs(m));
    end
    Alm = (1 - x.^2).^(abs(m) / 2) .* (factorial(abs(m)) / 2^l * Alm);
end