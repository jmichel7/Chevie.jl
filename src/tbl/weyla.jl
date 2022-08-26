chevieset(:A, :CartanMat, function (n,)
        local a, i
        a = IdentityMat(n)
        for i = 1:n
            (a[i])[i] = 2
            if i < n
                (a[i])[i + 1] = -1
            end
            if i > 1
                (a[i])[i - 1] = -1
            end
        end
        return a
    end)
chevieset(:A, :ReflectionDegrees, (n->begin
            2:n + 1
        end))
chevieset(:A, :ReflectionName, function (r, option)
        local o
        if haskey(option, :arg)
            return SPrint("\"A\",", r)
        elseif haskey(option, :TeX)
            if haskey(option, :Au)
                o = ["Z_2", "S_3", "S_4", "S_5"]
                return o[r]
            else
                return SPrint("A_", TeXBracket(r))
            end
        elseif haskey(option, :Au)
            o = ["Z2", "S3", "S4", "S5"]
            return o[r]
        else
            return SPrint("A", r)
        end
    end)
chevieset(:A, :GeneratingRoots, function (l,)
        local r, i
        r = map((i->begin
                        fill(0, max(0, (1 + (l + 1)) - 1))
                    end), 1:l)
        for i = 1:l
            (r[i])[[i, i + 1]] = [1, -1]
        end
        return r
    end)
chevieset(:A, :ParabolicRepresentatives, function (l, s)
        return (chevieget(:imp, :ParabolicRepresentatives))(1, 1, l, s)
    end)
chevieset(:A, :WordClass, function (pi,)
        local w, i, l, r
        w = []
        i = 0
        for l = pi
            r = mod(l, 2)
            w = Append(w, i + Concatenation(1:3 - 1:(l - 1) - r, 2:4 - 2:(l + r) - 2))
            i = i + l
        end
        return w
    end)
chevieset(:A, :ClassInfo, function (n,)
        local res
        res = Dict{Symbol, Any}(:classparams => partitions(n + 1))
        res[:classnames] = map(joindigits, res[:classparams])
        res[:classtext] = map(chevieget(:A, :WordClass), res[:classparams])
        res[:classes] = map((pi->begin
                        factorial(n + 1) // ((CharTableSymmetric[:centralizers])[1])(n, pi)
                    end), res[:classparams])
        res[:orders] = map(Lcm, res[:classparams])
        return res
    end)
chevieset(:A, :NrConjugacyClasses, (n->begin
            npartitions(n + 1)
        end))
chevieset(:A, :WeightInfo, function (n,)
        local M, i
        M = IdentityMat(n)
        for i = 1:n - 1
            (M[i])[n] = -i
        end
        return Dict{Symbol, Any}(:minusculeWeights => 1:n, :decompositions => map((i->begin
                                [(n + 1) - i]
                            end), 1:n), :moduli => [n + 1], :chosenAdaptedBasis => M)
    end)
chevieset(:A, :CharParams, (n->begin
            partitions(n + 1)
        end))
chevieset(:A, :LowestPowerFakeDegree, (p->begin
            p * (0:length(p) - 1)
        end))
chevieset(:A, :HighestPowerFakeDegree, (p->begin
            (Sum(p) * (Sum(p) - 1)) // 2 - (chevieget(:A, :LowestPowerFakeDegree))(conjugate_partition(p))
        end))
chevieset(:A, :CharInfo, function (n,)
        local res
        res = Dict{Symbol, Any}(:charparams => partitions(n + 1))
        res[:charnames] = map(joindigits, res[:charparams])
        res[:extRefl] = map((i->begin
                        Position(res[:charparams], Concatenation([(n + 1) - i], fill(0, max(0, (1 + i) - 1)) + 1))
                    end), 0:n)
        res[:b] = map((p->begin
                        p * (0:length(p) - 1)
                    end), res[:charparams])
        res[:B] = map(chevieget(:A, :HighestPowerFakeDegree), res[:charparams])
        res[:a] = res[:b]
        res[:A] = res[:B]
        return res
    end)
chevieset(:A, :PoincarePolynomial, function (n, param)
        return Product(1:n, (i->begin
                        Sum(0:i, (k->begin
                                    (-((param[1])[1]) // (param[1])[2]) ^ k
                                end))
                    end))
    end)
chevieset(:A, :SchurElement, function (n, alpha, param, sqrtparam)
        local i, j, lambda, res, q
        q = -((param[1])[1]) // (param[1])[2]
        lambda = βset(alpha)
        res = q ^ binomial(length(lambda), 3)
        for i = lambda
            for j = 0:i - 1
                if j in lambda
                    res = res // q ^ j
                else
                    res = res * Sum(0:(i - j) - 1, (e->begin
                                        q ^ e
                                    end))
                end
            end
        end
        return res
    end)
chevieset(:A, :FactorizedSchurElement, function (arg...,)
        return (chevieget(:imp, :FactorizedSchurElement))(1, 1, arg[1] + 1, [arg[2]], arg[3], [])
    end)
chevieset(:A, :HeckeRepresentation, function (n, param, sqrtparam, i)
        local H
        H = hecke(CoxeterGroup("A", n), -((param[1])[1]) // (param[1])[2])
        return -((param[1])[2]) * Spechtmodel(H, (partitions(n + 1))[i])
    end)
chevieset(:A, :Representation, function (n, i)
        return ((chevieget(:imp, :Representation))(1, 1, n + 1, i))[2:n + 1]
    end)
chevieset(:A, :FakeDegree, function (n, p, q)
        return (chevieget(:A, :PoincarePolynomial))(Sum(p) - 1, [[q, -1]]) // (chevieget(:A, :SchurElement))(Sum(p) - 1, p, [[q, -1]], [])
    end)
chevieset(:A, :DecompositionMatrix, function (l, p)
        return [[1:npartitions(l + 1), MatrixDecompositionMatrix(DecompositionMatrix(Specht(p, p), l + 1))]]
    end)
chevieset(:A, :UnipotentCharacters, function (l,)
        local ci
        ci = (chevieget(:A, :CharInfo))(l)
        return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:levi => [], :relativeType => Dict{Symbol, Any}(:series => "A", :indices => 1:l, :rank => l), :parameterExponents => fill(0, max(0, (1 + l) - 1)) + 1, :cuspidalName => "", :eigenvalue => 1, :charNumbers => 1:length(ci[:charparams]))], :families => map((i->begin
                                Family("C1", [i])
                            end), 1:length(ci[:charparams])), :charParams => ci[:charparams], :charSymbols => map((x->begin
                                [βset(x)]
                            end), ci[:charparams]), :a => ci[:a], :A => ci[:A])
    end)
chevieset(:A, :Ennola, function (l,)
        if l == 1
            return SPerm([-1, 2])
        else
            return SPerm()
        end
    end)
chevieset(:A, :Invariants, function (n,)
        local m
        m = ((CHEVIE[:A])[:GeneratingRoots])(n)
        return map((i->begin
                        function (arg...,)
                            local v
                            v = arg * m
                            return Sum(arrangements(1:n + 1, i), (a->begin
                                            Product(v[a])
                                        end))
                        end
                    end), 2:n + 1)
    end)
chevieset(:A, :UnipotentClasses, function (n, p)
        local uc, i, j, cl, d, ss, partition2parab
        uc = Dict{Symbol, Any}(:classes => map((p->begin
                                Dict{Symbol, Any}(:parameter => p)
                            end), partitions(n + 1)), :springerSeries => Concatenation(map((d->begin
                                    map((i->begin
                                                Dict{Symbol, Any}(:relgroup => CoxeterGroup("A", (n + 1) // d - 1), :Z => [E(d, i)], :levi => filter((i->begin
                                                                    mod(i, d) != 0
                                                                end), 1:n + 1), :locsys => [])
                                            end), prime_residues(d))
                                end), divisors(n + 1))))
        ss = (z->begin
                    First(uc[:springerSeries], (x->begin
                                x[:Z] == [z]
                            end))
                end)
        partition2parab(p) = begin
                local res, c, pa, i
                res = []
                c = 1
                for pa = p
                    for i = 1:pa - 1
                        push!(res, c)
                        c = c + 1
                    end
                    c = c + 1
                end
                return res
            end
        for i = 1:length(uc[:classes])
            cl = (uc[:classes])[i]
            p = cl[:parameter]
            d = gcd(p)
            cl[:name] = joindigits(p)
            cl[:Au] = crg(d, 1, 1)
            cl[:balacarter] = Concatenation(map((i->begin
                                Sum(p[1:i - 1]) + (1:p[i] - 1)
                            end), 1:length(p)))
            p = Concatenation(map((x->begin
                                1 - x:(3 - x) - (1 - x):x - 1
                            end), p))
            sort!(p)
            cl[:dynkin] = map((i->begin
                            p[i + 1] - p[i]
                        end), 1:length(p) - 1)
            cl[:red] = []
            p = 1
            for j = tally(cl[:parameter])
                cl[:red] = Append(cl[:red], p:(p + j[2]) - 2)
                p = p + j[2]
            end
            cl[:red] = ReflectionSubgroup(CoxeterGroup("A", p - 2), cl[:red])
            cl[:AuAction] = ExtendedReflectionGroup(cl[:red], [IdentityMat(rank(cl[:red]))])
            cl[:rep] = partition2parab(cl[:parameter])
            if d == 2
                push!((ss(1))[:locsys], [i, 2])
                push!((ss(-1))[:locsys], [i, 1])
            else
                for j = 0:d - 1
                    push!((ss(E(d, j)))[:locsys], [i, j + 1])
                end
            end
        end
        for ss = (uc[:springerSeries])[2:length(uc[:springerSeries])]
            ss[:hc] = 0
        end
        uc[:orderClasses] = hasse(Poset(map((x->begin
                                map((y->begin
                                            dominates(y[:parameter], x[:parameter])
                                        end), uc[:classes])
                            end), uc[:classes])))
        return uc
    end)
