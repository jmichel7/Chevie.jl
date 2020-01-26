
chevieset(:D, :CartanMat, function (n,)
        #= none:10 =#
        local a, m
        #= none:12 =#
        if n < 3
            #= none:12 =#
            m = 3
        else
            #= none:13 =#
            m = n
        end
        #= none:15 =#
        a = (chevieget(:A, :CartanMat))(m)
        #= none:17 =#
        (a[1:3])[1:3] = [[2, 0, -1], [0, 2, -1], [-1, -1, 2]]
        #= none:19 =#
        return (a[1:n])[1:n]
    end)
chevieset(:D, :Size, function (arg...,)
        #= none:4 =#
        return 2 ^ (arg[1] - 1) * factorial(arg[1])
    end)
chevieset(:D, :PrintDiagram, function (r, indices, title)
        #= none:3 =#
        local i, s
        #= none:5 =#
        print(title, " ", indices[1], "\n")
        #= none:6 =#
        s = pad("", length(title) + 1)
        #= none:8 =#
        print(s, " \\\n", s, "  ", indices[3])
        #= none:10 =#
        for i = 4:r
            #= none:10 =#
            print(" - ", indices[i])
        end
        #= none:11 =#
        print("\n")
        #= none:13 =#
        print(s, " /\n", s, indices[2], "\n")
    end)
chevieset(:D, :GeneratingRoots, function (l,)
        #= none:3 =#
        local r, rts, i
        #= none:5 =#
        rts = []
        #= none:7 =#
        for i = 1:l - 1
            #= none:7 =#
            r = 0 * (1:l)
            #= none:8 =#
            r[[i, i + 1]] = [1, -1]
            #= none:9 =#
            push!(rts, r)
        end
        #= none:11 =#
        r = 0 * (1:l)
        #= none:12 =#
        r[[l - 1, l]] = [1, 1]
        #= none:13 =#
        push!(rts, r)
        #= none:15 =#
        return reverse(rts)
    end)
chevieset(:D, :WeightInfo, function (n,)
        #= none:3 =#
        local res
        #= none:5 =#
        if mod(n, 2) == 1
            #= none:5 =#
            return Dict{Symbol, Any}(:minusculeWeights => [1, 2, n], :decompositions => [[1], [3], [2]], :moduli => [4])
        else
            #= none:8 =#
            return Dict{Symbol, Any}(:minusculeWeights => [1, 2, n], :decompositions => [[1, 0], [0, 1], [1, 1]], :moduli => [2, 2])
        end
    end)
chevieset(:D, :ParabolicRepresentatives, function (l, s)
        #= none:4 =#
        return (chevieget(:imp, :ParabolicRepresentatives))(2, 2, l, s)
    end)
chevieset(:D, :WordsClassRepresentatives, function (arg...,)
        #= none:5 =#
        local n, param, res, w, i, pi, l, r
        #= none:7 =#
        n = arg[1]
        #= none:9 =#
        if length(arg) == 2
            #= none:9 =#
            param = map((a->begin
                            #= none:9 =#
                            map(copy, a)
                        end), arg[2])
        else
            #= none:11 =#
            param = partition_tuples(n, 2)
        end
        #= none:13 =#
        res = []
        #= none:17 =#
        for pi = param
            #= none:19 =#
            if pi[2] == '+'
                #= none:19 =#
                pi[2] = []
            end
            #= none:21 =#
            if IsList(pi[2]) && mod(length(pi[2]), 2) == 0
                #= none:22 =#
                w = []
                #= none:23 =#
                i = 1
                #= none:27 =#
                for l = reverse(pi[2])
                    #= none:30 =#
                    if i == 1
                        #= none:30 =#
                        w = Append(w, 2:(i + l) - 1)
                    else
                        #= none:32 =#
                        w = Append(w, i:(i - 1) - i:3)
                        #= none:33 =#
                        w = Append(w, 1:(i + l) - 1)
                    end
                    #= none:35 =#
                    i = i + l
                end
                #= none:39 =#
                for l = pi[1]
                    #= none:40 =#
                    r = mod(l, 2)
                    #= none:42 =#
                    w = Append(w, i + Concatenation(1:3 - 1:(l - 1) - r, 2:4 - 2:(l + r) - 2))
                    #= none:44 =#
                    i = i + l
                end
                #= none:48 =#
                if w != [] && w[1] == 2
                    #= none:48 =#
                    w[1] = 1
                end
                #= none:53 =#
                if pi[2] == [] && all((x->begin
                                    #= none:53 =#
                                    mod(x, 2) == 0
                                end), pi[1])
                    #= none:54 =#
                    push!(res, w)
                    #= none:55 =#
                    w = copy(w)
                    #= none:56 =#
                    w[1] = 2
                end
                #= none:58 =#
                push!(res, w)
            end
        end
        #= none:62 =#
        return res
    end)
chevieset(:D, :ClassInfo, function (n,)
        #= none:20 =#
        local res
        #= none:22 =#
        res = (chevieget(:imp, :ClassInfo))(2, 2, n)
        #= none:24 =#
        res[:classparams] = map(function (x,)
                    #= none:25 =#
                    if length(x) == 2
                        #= none:25 =#
                        return x
                    end
                    #= none:27 =#
                    if x[3] == 0
                        #= none:27 =#
                        return [x[1], '+']
                    else
                        #= none:28 =#
                        return [x[1], '-']
                    end
                end, res[:classparams])
        #= none:31 =#
        res[:classtext] = (chevieget(:D, :WordsClassRepresentatives))(n, res[:classparams])
        #= none:33 =#
        return res
    end)
chevieset(:D, :NrConjugacyClasses, function (n,)
        #= none:4 =#
        if mod(n, 2) == 1
            #= none:4 =#
            return NrPartitionTuples(n, 2) // 2
        else
            #= none:6 =#
            return (NrPartitionTuples(n, 2) + 3 * NrPartitions(n // 2)) // 2
        end
    end)
chevieset(:D, :CharInfo, (n->begin
            #= none:7 =#
            (chevieget(:imp, :CharInfo))(2, 2, n)
        end))
chevieset(:D, :CharName, function (arg...,)
        #= none:5 =#
        return PartitionTupleToString(arg[2])
    end)
chevieset(:D, :gensMODA, [nothing, nothing, nothing, [[perm"(1,2)(7,8)", perm"(3,4)(5,6)", perm"(2,3)(6,7)", perm"(3,5)(4,6)"], [[4], [nothing, nothing, 2]], [[2], [1, nothing, 1]]], nothing, [[perm"(1,2)(8,11)(12,14)(15,17)(16,18)(19,21)(22,25)(31,32)", perm"(3,4)(5,6)(7,9)(10,13)(20,23)(24,26)(27,28)(29,30)", perm"(2,3)(6,8)(9,12)(13,16)(17,20)(21,24)(25,27)(30,31)", perm"(3,5)(4,6)(12,15)(14,17)(16,19)(18,21)(27,29)(28,30)", perm"(5,7)(6,9)(8,12)(11,14)(19,22)(21,25)(24,27)(26,28)", perm"(7,10)(9,13)(12,16)(14,18)(15,19)(17,21)(20,24)(23,26)"], [[16], [4, nothing, 6], [1, nothing, nothing, nothing, 5]], [[12], [2, nothing, 6], [nothing, 2, nothing, nothing, 4]]], nothing, [[perm"(1,2)(8,11)(12,15)(16,20)(17,21)(22,26)(23,27)(28,33)(29,34)(30,35)(36,41)(37,42)(43,50)(44,51)(52,59)(60,68)(61,69)(70,77)(78,85)(79,86)(87,92)(88,93)(94,99)(95,100)(96,101)(102,106)(103,107)(108,112)(109,113)(114,117)(118,121)(127,128)", perm"(3,4)(5,6)(7,9)(10,13)(14,18)(19,24)(25,31)(32,38)(39,45)(40,46)(47,53)(48,54)(49,55)(56,62)(57,63)(58,64)(65,71)(66,72)(67,73)(74,80)(75,81)(76,82)(83,89)(84,90)(91,97)(98,104)(105,110)(111,115)(116,119)(120,122)(123,124)(125,126)", perm"(2,3)(6,8)(9,12)(13,17)(18,23)(20,25)(24,30)(26,32)(33,39)(34,40)(41,48)(42,49)(50,57)(51,58)(53,61)(59,67)(62,70)(68,76)(71,78)(72,79)(80,87)(81,88)(89,95)(90,96)(97,103)(99,105)(104,109)(106,111)(112,116)(117,120)(121,123)(126,127)", perm"(3,5)(4,6)(12,16)(15,20)(17,22)(21,26)(23,29)(27,34)(30,37)(35,42)(39,47)(45,53)(48,56)(54,62)(57,65)(58,66)(63,71)(64,72)(67,75)(73,81)(76,84)(82,90)(87,94)(92,99)(95,102)(100,106)(103,108)(107,112)(109,114)(113,117)(123,125)(124,126)", perm"(5,7)(6,9)(8,12)(11,15)(22,28)(26,33)(29,36)(32,39)(34,41)(37,44)(38,45)(40,48)(42,51)(46,54)(49,58)(55,64)(65,74)(71,80)(75,83)(78,87)(81,89)(84,91)(85,92)(88,95)(90,97)(93,100)(96,103)(101,107)(114,118)(117,121)(120,123)(122,124)", perm"(7,10)(9,13)(12,17)(15,21)(16,22)(20,26)(25,32)(31,38)(36,43)(41,50)(44,52)(48,57)(51,59)(54,63)(56,65)(58,67)(62,71)(64,73)(66,75)(70,78)(72,81)(77,85)(79,88)(86,93)(91,98)(97,104)(103,109)(107,113)(108,114)(112,117)(116,120)(119,122)", perm"(10,14)(13,18)(17,23)(21,27)(22,29)(26,34)(28,36)(32,40)(33,41)(38,46)(39,48)(45,54)(47,56)(52,60)(53,62)(59,68)(61,70)(67,76)(69,77)(73,82)(75,84)(81,90)(83,91)(88,96)(89,97)(93,101)(95,103)(100,107)(102,108)(106,112)(111,116)(115,119)", perm"(14,19)(18,24)(23,30)(27,35)(29,37)(34,42)(36,44)(40,49)(41,51)(43,52)(46,55)(48,58)(50,59)(54,64)(56,66)(57,67)(62,72)(63,73)(65,75)(70,79)(71,81)(74,83)(77,86)(78,88)(80,89)(85,93)(87,95)(92,100)(94,102)(99,106)(105,111)(110,115)"], [[64], [16, nothing, 24], [nothing, nothing, 32], [4, nothing, nothing, nothing, 20], [nothing, nothing, nothing, nothing, nothing, nothing, 16]], [[56], [12, nothing, 24], [6, nothing, 28], [2, 4, nothing, nothing, 18], [1, nothing, 3, nothing, nothing, nothing, 14]]]])
chevieset(:D, :ClassParameter, function (n, w)
        #= none:5 =#
        local x, i, res, mark, cyc, j, tmp, gens
        #= none:8 =#
        x = Perm()
        #= none:10 =#
        for i = w
            #= none:11 =#
            if i == 1
                #= none:11 =#
                x = x * (Perm(1, n + 2))(2, n + 1)
            else
                #= none:12 =#
                x = x * (Perm(i - 1, i))((i - 1) + n, i + n)
            end
        end
        #= none:15 =#
        res = [[], []]
        #= none:17 =#
        mark = 1:n
        #= none:19 =#
        for i = 1:n
            #= none:20 =#
            if mark[i] != 0
                #= none:21 =#
                cyc = CyclePermInt(x, i)
                #= none:23 =#
                if i + n in cyc
                    #= none:23 =#
                    push!(res[2], length(cyc) // 2)
                else
                    #= none:25 =#
                    push!(res[1], length(cyc))
                end
                #= none:27 =#
                for j = cyc
                    #= none:28 =#
                    if j > n
                        #= none:28 =#
                        mark[j - n] = 0
                    else
                        #= none:30 =#
                        mark[j] = 0
                    end
                end
            end
        end
        #= none:34 =#
        if res[2] == [] && all((i->begin
                            #= none:34 =#
                            mod(i, 2) == 0
                        end), res[1])
            #= none:40 =#
            if !((CHEVIE.R("gensMODA", "D"))[n] !== nothing)
                #= none:41 =#
                tmp = CoxeterGroup("D", n)
                #= none:43 =#
                gens = PermCosetsSubgroup(tmp, ReflectionSubgroup(tmp, 2:n))
                #= none:46 =#
                tmp = (chevieget(:D, :ClassInfo))(n)
                #= none:48 =#
                tmp = (tmp[:classtext])[Filtered(1:length(tmp[:classnames]), (i->(#= none:48 =#
                                    '+' in (tmp[:classnames])[i] || '-' in (tmp[:classnames])[i])))]
                #= none:51 =#
                tmp = map((a->begin
                                #= none:51 =#
                                CycleStructurePerm(Product(gens[a]))
                            end), tmp)
                #= none:53 =#
                (chevieget(:D, :gensMODA))[n] = [gens, tmp[2 * (1:length(tmp) // 2) - 1], tmp[2 * (1:length(tmp) // 2)]]
            end
            #= none:57 =#
            tmp = CycleStructurePerm(Product((((chevieget(:D, :gensMODA))[n])[1])[w]))
            #= none:59 =#
            if tmp in ((chevieget(:D, :gensMODA))[n])[2] && !tmp in ((chevieget(:D, :gensMODA))[n])[3]
                #= none:61 =#
                res[2] = '+'
            elseif #= none:63 =# !tmp in ((chevieget(:D, :gensMODA))[n])[2] && tmp in ((chevieget(:D, :gensMODA))[n])[3]
                #= none:65 =#
                res[2] = '-'
            end
        end
        #= none:69 =#
        sort!(res[1])
        #= none:71 =#
        if IsList(res[2])
            #= none:72 =#
            sort!(res[2])
            #= none:74 =#
            return [reverse(res[1]), reverse(res[2])]
        else
            #= none:77 =#
            return [reverse(res[1]), res[2]]
        end
    end)
chevieset(:D, :FactorizedSchurElement, function (arg...,)
        #= none:3 =#
        local p, i, n
        #= none:5 =#
        p = arg[2]
        #= none:6 =#
        n = arg[1]
        #= none:8 =#
        if p[2] in "+-"
            #= none:8 =#
            p = [p[1], p[1]]
        end
        #= none:10 =#
        return (chevieget(:imp, :FactorizedSchurElement))(2, 2, n, p, arg[3], [])
    end)
chevieset(:D, :HeckeRepresentation, function (arg...,)
        #= none:3 =#
        local p, i, n
        #= none:5 =#
        i = arg[4]
        #= none:6 =#
        n = arg[1]
        #= none:8 =#
        p = (((chevieget(:D, :CharInfo))(n))[:charparams])[i]
        #= none:10 =#
        if p[length(p)] == 0
            #= none:10 =#
            i = i + 1
        elseif #= none:11 =# p[length(p)] == 1
            #= none:11 =#
            i = i - 1
        end
        #= none:13 =#
        return (chevieget(:imp, :HeckeRepresentation))(2, 2, n, arg[2], [], i)
    end)
chevieset(:D, :Representation, function (n, i)
        #= none:3 =#
        local p
        #= none:5 =#
        p = (((chevieget(:D, :CharInfo))(n))[:charparams])[i]
        #= none:7 =#
        if p[length(p)] == 0
            #= none:7 =#
            i = i + 1
        elseif #= none:8 =# p[length(p)] == 1
            #= none:8 =#
            i = i - 1
        end
        #= none:10 =#
        return (chevieget(:imp, :Representation))(2, 2, n, i)
    end)
chevieset(:D, :PoincarePolynomial, function (n, para)
        #= none:7 =#
        local q
        #= none:9 =#
        q = -((para[1])[1]) // (para[1])[2]
        #= none:11 =#
        return Sum(0:n - 1, (k->begin
                            #= none:11 =#
                            q ^ k
                        end)) * Product(1:n - 1, (i->begin
                            #= none:11 =#
                            (q ^ i + 1) * Sum(0:i - 1, (k->begin
                                            #= none:11 =#
                                            q ^ k
                                        end))
                        end))
    end)
chevieset(:D, :symbolcharparam, (c->begin
            #= none:4 =#
            SymbolPartitionTuple(c, 0)
        end))
chevieset(:D, :Invariants, function (n,)
        #= none:3 =#
        local m
        #= none:5 =#
        m = (chevieget(:imp, :GeneratingRoots))(2, 2, n)
        #= none:7 =#
        return map((f->begin
                        #= none:7 =#
                        function (arg...,)
                            #= none:8 =#
                            return ApplyFunc(f, arg * m)
                        end
                    end), ((CHEVIE[:imp])[:Invariants])(2, 2, n))
    end)
chevieset(:D, :CycPolGenericDegree, (c->begin
            #= none:13 =#
            CycPolGenericDegreeSymbol(SymbolPartitionTuple(c, 0))
        end))
chevieset(:D, :SchurElement, function (n, phi, q, sqrtparam)
        #= none:4 =#
        return (chevieget(:D, :PoincarePolynomial))(n, q) // Value((chevieget(:D, :CycPolGenericDegree))(phi), -((q[1])[1]) // (q[1])[2])
    end)
chevieset(:D, :FakeDegree, function (n, c, q)
        #= none:8 =#
        return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c, 0)), q)
    end)
chevieset(:D, :UnipotentCharacters, function (rank,)
        #= none:3 =#
        local uc, symbols, r, d, s
        #= none:5 =#
        uc = Dict{Symbol, Any}(:harishChandra => [], :charSymbols => [])
        #= none:7 =#
        for d = 4 * (0:RootInt(div(rank, 4), 2))
            #= none:8 =#
            r = div(d ^ 2, 4)
            #= none:10 =#
            s = Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "B", :indices => 1 + r:rank, :rank => rank - r), :levi => 1:r, :eigenvalue => (-1) ^ div(d + 1, 4), :parameterExponents => Concatenation([d], fill(0, max(0, (1 + rank) - (2 + r))) + 1))
            #= none:14 =#
            if r < 10
                #= none:14 =#
                s[:cuspidalName] = SPrint("D_", r, "")
            else
                #= none:16 =#
                s[:cuspidalName] = SPrint("D_{", r, "}")
            end
            #= none:18 =#
            if d == 0
                #= none:19 =#
                (s[:relativeType])[:series] = "D"
                #= none:21 =#
                s[:cuspidalName] = ""
                #= none:23 =#
                (s[:parameterExponents])[1] = 1
            end
            #= none:25 =#
            push!(uc[:harishChandra], s)
            #= none:27 =#
            symbols = BDSymbols(rank, d)
            #= none:29 =#
            s[:charNumbers] = (1:length(symbols)) + length(uc[:charSymbols])
            #= none:31 =#
            FixRelativeType(s)
            #= none:33 =#
            uc[:charSymbols] = Append(uc[:charSymbols], symbols)
        end
        #= none:35 =#
        uc[:a] = map(LowestPowerGenericDegreeSymbol, uc[:charSymbols])
        #= none:37 =#
        uc[:A] = map(HighestPowerGenericDegreeSymbol, uc[:charSymbols])
        #= none:39 =#
        uc[:families] = FamiliesClassical(uc[:charSymbols])
        #= none:41 =#
        return uc
    end)
chevieset(:D, :ReflectionDegrees, (n->begin
            #= none:3 =#
            Concatenation(2 * (1:n - 1), [n])
        end))
