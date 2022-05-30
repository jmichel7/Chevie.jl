
chevieset(Symbol("2A"), :WordsClassRepresentatives, function (arg...,)
        local n, part, guesslongest, redw, l, w0, p
        n = arg[1]
        if length(arg) > 1
            part = arg[2]
        else
            part = partitions(n + 1)
        end
        redw(n, w) = begin
                local l, i
                l = []
                while true
                    i = PositionProperty(1:n, (j->begin
                                    j ^ w > (j + 1) ^ w
                                end))
                    if i == false
                        return l
                    end
                    push!(l, i)
                    w = Perm(i, i + 1) * w
                end
                return l
            end
        guesslongest(p) = begin
                local x, off, i
                p = Concatenation(filter((i->begin
                                    mod(i, 2) == 0
                                end), p), filter((i->begin
                                    i != 1 && mod(i, 2) == 1
                                end), p))
                x = Perm()
                off = 0
                for i = p
                    x = x * Product(map((j->begin
                                            Perm(l[off + 1], l[off + j])
                                        end), 2:i))
                    off = off + i
                end
                return x
            end
        l = []
        w0 = Perm()
        for p = 1:div(n + 1, 2)
            l = Append(l, [p, (n - p) + 2])
            w0 = w0 * Perm(p, (n - p) + 2)
        end
        if mod(n, 2) == 0
            push!(l, div(n, 2) + 1)
        end
        return map((p->begin
                        redw(n, guesslongest(p) * w0)
                    end), part)
    end)
chevieset(Symbol("2A"), :ClassInfo, function (n,)
        local res
        res = (chevieget(:A, :ClassInfo))(n)
        res[:classtext] = (chevieget(Symbol("2A"), :WordsClassRepresentatives))(n, res[:classparams])
        delete!(res, :orders)
        return res
    end)
chevieset(Symbol("2A"), :NrConjugacyClasses, (n->begin
            npartitions(n + 1)
        end))
chevieset(Symbol("2A"), :CharParams, (n->begin
            partitions(n + 1)
        end))
chevieset(Symbol("2A"), :CharInfo, (n->begin
            (chevieget(:A, :CharInfo))(n)
        end))
chevieset(Symbol("2A"), :PhiFactors, (n->begin
            map((x->begin
                        (-1) ^ x
                    end), 2:n + 1)
        end))
chevieset(Symbol("2A"), :HeckeRepresentation, function (n, param, sqrtparam, i)
        local H, res, W, p
        W = CoxeterGroup("A", n)
        H = hecke(W, -((param[1])[1]) // (param[1])[2])
        p = (partitions(n + 1))[i]
        res = Dict{Symbol, Any}(:gens => Spechtmodel(H, p))
        res[:F] = Product((res[:gens])[LongestCoxeterWord(W)]) // GetRoot((HeckeCentralMonomials(H))[i]) * (-1) ^ (chevieget(:A, :LowestPowerFakeDegree))(p)
        return res
    end)
chevieset(Symbol("2A"), :Representation, function (n, i)
        return (chevieget(Symbol("2A"), :HeckeRepresentation))(n, map((x->begin
                            [1, -1]
                        end), 1:n), 1:n * 0 + 1, i)
    end)
PartitionTwoCoreQuotient(d, p) = begin
        local x
        x = symbol_partition_tuple(reverse(p), -d)
        return partβ(gapSet(Concatenation(2 * x[1], 2 * x[2] + 1)))
    end
chevieset(Symbol("2A"), :UnipotentCharacters, function (l,)
        local uc, d, k, s, i, r
        uc = (chevieget(:A, :UnipotentCharacters))(l)
        uc[:charSymbols] = map((i->begin
                        [i]
                    end), (chevieget(:A, :CharParams))(l))
        uc[:almostHarishChandra] = uc[:harishChandra]
        ((uc[:almostHarishChandra])[1])[:relativeType] = Dict{Symbol, Any}(:orbit => [Dict{Symbol, Any}(:series => "A", :indices => 1:l, :rank => l)], :twist => Product(1:div(l, 2), (i->begin
                                Perm(i, (l + 1) - i)
                            end)), :rank => l)
        uc[:harishChandra] = []
        d = 0
        while (d * (d + 1)) // 2 <= l + 1
            k = (l + 1) - div(d * (d + 1), 2)
            if mod(k, 2) == 0
                r = div(k, 2)
                s = Dict{Symbol, Any}(:levi => r + 1:l - r, :relativeType => Dict{Symbol, Any}(:series => "B", :indices => r:(r - 1) - r:1, :rank => r), :eigenvalue => (-1) ^ (Product(d + (-1:2)) // 8))
                if d == 0
                    (s[:relativeType])[:cartanType] = 1
                end
                if r != 0
                    s[:parameterExponents] = Concatenation([2d + 1], fill(0, max(0, (1 + r) - 2)) + 2)
                else
                    s[:parameterExponents] = []
                end
                if k < l
                    if l - k < 10
                        s[:cuspidalName] = SPrint("{}^2A_", l - k, "")
                    else
                        s[:cuspidalName] = SPrint("{}^2A_{", l - k, "}")
                    end
                else
                    s[:cuspidalName] = ""
                end
                s[:charNumbers] = map((a->begin
                                Position(uc[:charSymbols], [PartitionTwoCoreQuotient(d, a)])
                            end), (chevieget(:B, :CharParams))(r))
                FixRelativeType(s)
                push!(uc[:harishChandra], s)
            end
            d = d + 1
        end
        for i = 1:length(uc[:families])
            if 0 != mod((uc[:a])[i] + (uc[:A])[i], 2)
                (uc[:families])[i] = Family("C'1", ((uc[:families])[i])[:charNumbers])
            end
        end
        return uc
    end)