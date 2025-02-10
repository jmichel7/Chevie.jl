chevieset(:imp, :UnipotentCharacters, function (p, q, r)
        local uc, cusp, f, l, ci, seteig, s, extra, addextra
        if !(q in [1, p])
            return false
        end
        uc = Dict{Symbol, Any}(:charSymbols => (chevieget(:imp, :CharSymbols))(p, q, r))
        uc[:a] = map(valuation_gendeg_symbol, uc[:charSymbols])
        uc[:A] = map(degree_gendeg_symbol, uc[:charSymbols])
        ci = (chevieget(:imp, :CharInfo))(p, q, r)
        if q == 1
            cusp = gapSet(map((S->begin
                                map(length, S) - Minimum(map(length, S))
                            end), uc[:charSymbols]))
            cusp = map((x->begin
                            map((y->begin
                                        0:y - 1
                                    end), x)
                        end), cusp)
            SortBy(cusp, ranksymbol)
            uc[:harishChandra] = map(function (c,)
                        local cr, res
                        cr = ranksymbol(c)
                        res = Dict{Symbol, Any}(:levi => 1:cr)
                        if cr < r
                            res[:parameterExponents] = [map(length, c)]
                        else
                            res[:parameterExponents] = []
                        end
                        res[:parameterExponents] = Append(res[:parameterExponents], fill(0, max(0, (1 + r) - (2 + cr))) + 1)
                        if r == cr
                            res[:relativeType] = Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0)
                        else
                            res[:relativeType] = Dict{Symbol, Any}(:series => "ST", :indices => 1 + cr:r, :rank => r - cr, :p => p, :q => 1)
                        end
                        res[:eigenvalue] = E(24, -2 * (p ^ 2 - 1) * div(Sum(c, length), p)) * E(2p, Sum(0:p - 1, (i->begin
                                                -((i ^ 2 + p * i)) * length(c[i + 1])
                                            end)))
                        res[:charNumbers] = map((x->begin
                                        Position(uc[:charSymbols], symbol_partition_tuple(x, map(length, c)))
                                    end), map((x->begin
                                            map(partβ, x)
                                        end), ((chevieget(:imp, :CharSymbols))(p, 1, r - cr))[1:length(partition_tuples(r - cr, p))]))
                        res[:cuspidalName] = ImprimitiveCuspidalName(c)
                        return res
                    end, cusp)
            uc[:b] = uc[:a] * 0
            uc[:B] = uc[:a] * 0
            (uc[:b])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:b]
            (uc[:B])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:B]
            uc[:families] = map((y->begin
                            MakeFamilyImprimitive(y, uc)
                        end), CollectBy(uc[:charSymbols], (x->begin
                                tally(Concatenation(x))
                            end)))
            SortBy(uc[:families], (x->begin
                        x[:charNumbers]
                    end))
            if r == 1
                l = map(function (S,)
                            local p
                            p = Position(S, [])
                            if p == false
                                return 1
                            else
                                return (-1) ^ p
                            end
                        end, (uc[:charSymbols])[((uc[:families])[2])[:charNumbers]])
                ((uc[:families])[2])[:fourierMat] = ((uc[:families])[2])[:fourierMat] ^ DiagonalMat(l)
                uc[:cyclicparam] = map(function (s,)
                            if count((x->begin
                                                x == 1
                                            end), Flat(s)) == 1
                                return [1]
                            else
                                s = deepcopy(s)
                                l = PositionProperty(s, (p->begin
                                                1 in p
                                            end))
                                s[l] = []
                                return [PositionProperty(s, (p->begin
                                                        1 in p
                                                    end)) - 1, l - 1]
                            end
                        end, uc[:charSymbols])
            elseif r == 2 && p == 3
                ((uc[:families])[4])[:fourierMat] = ((uc[:families])[4])[:fourierMat] ^ DiagonalMat(-1, 1, 1)
                ((uc[:families])[1])[:fourierMat] = ((uc[:families])[1])[:fourierMat] ^ DiagonalMat(1, -1, -1, 1, 1, 1, 1, 1, 1)
            end
            return uc
        elseif p == q
            uc[:families] = []
            for f = CollectBy(1:length(uc[:charSymbols]), (i->begin
                                tally(Concatenation(fullsymbol((uc[:charSymbols])[i])))
                            end))
                if length(gapSet(map(fullsymbol, (uc[:charSymbols])[f]))) > 1
                    push!(uc[:families], Dict{Symbol, Any}(:charNumbers => f))
                else
                    uc[:families] = Append(uc[:families], map((x->begin
                                        Family("C1", [x])
                                    end), f))
                end
            end
            SortBy(uc[:families], (x->begin
                        x[:charNumbers]
                    end))
            uc[:harishChandra] = map((l->begin
                            Dict{Symbol, Any}(:charNumbers => l)
                        end), CollectBy(1:length(uc[:charSymbols]), function (i,)
                            local s, l
                            s = fullsymbol((uc[:charSymbols])[i])
                            l = map(length, s)
                            return [Sum(s, (x->begin
                                                Sum(partβ(x))
                                            end)), l - Minimum(l)]
                        end))
            SortBy(uc[:harishChandra], (x->begin
                        x[:charNumbers]
                    end))
            extra = []
            for f = uc[:harishChandra]
                addextra = false
                s = fullsymbol((uc[:charSymbols])[(f[:charNumbers])[1]])
                l = r - Sum(s, (x->begin
                                    Sum(partβ(x))
                                end))
                f[:levi] = 1:l
                s = map(length, s)
                s = s - Minimum(s)
                f[:eigenvalue] = E(24, -2 * (p ^ 2 - 1) * div(Sum(s), p)) * E(2p, Sum(0:p - 1, (i->begin
                                        -((i ^ 2 + p * i)) * s[i + 1]
                                    end)))
                if l == r
                    f[:relativeType] = Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0)
                    f[:parameterExponents] = []
                    if length(f[:charNumbers]) == 2
                        addextra = true
                    end
                elseif l == 0
                    f[:relativeType] = Dict{Symbol, Any}(:series => "ST", :indices => 1:r, :rank => r, :p => p, :q => q)
                    f[:parameterExponents] = fill(0, max(0, (1 + r) - 1)) + 1
                else
                    f[:relativeType] = Dict{Symbol, Any}(:series => "ST", :indices => l + 1:r, :rank => r - l, :p => p, :q => 1)
                    f[:parameterExponents] = Concatenation([s], fill(0, max(0, (1 + ((r - l) - 1)) - 1)) + 1)
                end
                s = map((x->begin
                                0:x - 1
                            end), s)
                f[:cuspidalName] = ImprimitiveCuspidalName(s)
                if addextra
                    s = deepcopy(f[:charNumbers])
                    f[:charNumbers] = s[[1]]
                    f = deepcopy(f)
                    f[:charNumbers] = s[[2]]
                    push!(f[:cuspidalName], '2')
                    push!(extra, f)
                end
            end
            uc[:harishChandra] = Append(uc[:harishChandra], extra)
            for f = uc[:families]
                f[:eigenvalues] = map((i->begin
                                (First(uc[:harishChandra], (s->(i in s[:charNumbers];))))[:eigenvalue]
                            end), f[:charNumbers])
            end
            uc[:b] = fill(0, max(0, (1 + length(uc[:charSymbols])) - 1))
            uc[:B] = fill(0, max(0, (1 + length(uc[:charSymbols])) - 1))
            (uc[:b])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:b]
            (uc[:B])[((uc[:harishChandra])[1])[:charNumbers]] = ci[:B]
            uc[:families] = map((x->begin
                            MakeFamilyImprimitive((uc[:charSymbols])[x[:charNumbers]], uc)
                        end), uc[:families])
            if [p, q, r] == [3, 3, 3]
                uc[:curtis] = [1, 2, 3, 7, 8, 10, 4, 5, 9, 6, -12, -11]
            end
            return uc
        end
    end)
chevieset(:imp, :InitHeckeBasis, function (p, q, r, H)
        if q != 1 || r != 1
            error("implemented only for G(d,1,1)")
        end
        H[:mul] = function (x, y)
                local H, res, ops, W, temp, i, xi, temp1, j, e, pol, d
                H = hecke(y)
                if !(IsRec(x)) || (!(haskey(x, :hecke)) || !(haskey(x, :elm)))
                    if x == x * 0
                        return HeckeElt(H, y[:basis], [], [])
                    else
                        return HeckeElt(H, y[:basis], y[:elm], y[:coeff] * (x * H[:unit]))
                    end
                end
                if !(IsIdentical(H, hecke(x)))
                    error("not elements of the same algebra")
                end
                ops = H[:operations]
                pol = Coefficients(Product((H[:parameter])[1], (u->begin
                                    Mvp("xxx") - u
                                end)), "xxx")
                d = length(pol) - 1
                if x[:basis] != y[:basis]
                    return (basis(H, "T"))(x) * (basis(H, "T"))(y)
                elseif x[:basis] == "T"
                    W = Group(H)
                    res = HeckeElt(H, x[:basis], [], [])
                    for i = 1:length(x[:elm])
                        temp = (x[:coeff])[i] * y
                        xi = (x[:elm])[i]
                        for i = 1:length(xi)
                            temp1 = HeckeElt(H, x[:basis], [], [])
                            for j = 1:length(temp[:elm])
                                e = length((temp[:elm])[j])
                                if e + 1 < d
                                    push!(temp1[:elm], map((i->begin
                                                    1
                                                end), 1:e + 1))
                                    push!(temp1[:coeff], (temp[:coeff])[j])
                                else
                                    temp1[:elm] = Append(temp1[:elm], map((i->begin
                                                        fill(0, max(0, (1 + i) - 1)) + 1
                                                    end), 0:d - 1))
                                    temp1[:coeff] = Append(temp1[:coeff], -(pol[1:d]) * (temp[:coeff])[j])
                                end
                            end
                            temp = temp1
                            CollectCoefficients(temp)
                        end
                        res = res + temp
                    end
                    return res
                else
                    return (basis(H, x[:basis]))((basis(H, "T"))(x) * (basis(H, "T"))(y))
                end
            end
        H[:inverse] = function (h,)
                local H, d, pol
                if length(h[:elm]) != 1
                    error("inverse implemented only for single T_w")
                end
                H = hecke(h)
                pol = Coefficients(Product((H[:parameter])[1], (u->begin
                                    Mvp("xxx") - u
                                end)), "xxx")
                d = length(pol) - 1
                return (h[:coeff])[1] ^ -1 * (basis(H, "T"))(map((i->begin
                                            fill(0, max(0, (1 + i) - 1)) + 1
                                        end), 0:d - 1), -(pol[2:d + 1]) // pol[1]) ^ length((h[:elm])[1])
            end
        return true
    end)
