function Anm = AssociatedLaguerre(n,m,x)
Anm = 0;
    for i = 0 : n
        Anm = Anm + factorial(m + n) * nchoosek(m + n, n - i) / factorial(i) * (-x).^i;
    end
end