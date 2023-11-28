    @kwdef struct Flatte1405
        m0::Float64
        Γ0::Float64
        L::Int
        l::Int
    end
    function (lineshape::Flatte1405)(σ)
        (; m0, Γ0, L, l) = lineshape
        l != 0 && error("Not intended to be calle with l!=0, l=$l")
        _q = q(σ)
        _q0 = q(m0^2)
        m1, m2 = mK, mp
        mπ, mΣ = 0.14, 1.197
        _p = breakup(σ, m1^2, m2^2)
        _p′, _p0′ = breakup(σ, mπ^2, mΣ^2), breakup(m0^2, mπ^2, mΣ^2)
        Γ1 = Γ0 * (_p / _p0′) * m0 / sqrt(σ)
        Γ2 = Γ0 * (_p′ / _p0′) * m0 / sqrt(σ)
        Γ = Γ1 + Γ2
        R5 = 5
        R1 = 1.5 # 1/GeV
        return BW(σ, m0, Γ) *
               (_q / _q0)^L * h(_q * R5, L) / h(_q0 * R5, L)
    end
