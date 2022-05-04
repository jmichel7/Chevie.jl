
chevieset(:G2, :CartanMat, function (arg...,)
        local a, type_
        if length(arg) > 0
            type_ = arg[1]
        else
            type_ = 1
        end
        return [[2, -type_], [-3 // type_, 2]]
    end)
chevieset(:G2, :ReflectionName, function (arg...,)
        local i, opt, type_
        if length(arg) == 1
            return "G2(?)"
        end
        type_ = arg[2]
        opt = arg[1]
        if type_ == 1
            if haskey(opt, :TeX)
                return "G_2"
            elseif haskey(opt, :arg)
                return "\"G\",2"
            else
                return "G2"
            end
        elseif type_ == root(3)
            if haskey(opt, :TeX)
                return "G_{\\hbox{sym}2}"
            elseif haskey(opt, :arg)
                return "\"Gsym\",2"
            else
                return "Gsym2"
            end
        elseif haskey(opt, :TeX)
            return SPrint("G_2(", Format(type_ ^ 2 // 3, opt), ")")
        elseif haskey(opt, :arg)
            return SPrint("\"G\",", 2, ",", Format(type_ ^ 2 // 3, opt))
        else
            return SPrint("G2(", Format(type_ ^ 2 // 3, opt), ")")
        end
    end)
chevieset(:G2, :ParabolicRepresentatives, (s->begin
            (chevieget(:imp, :ParabolicRepresentatives))(6, 6, 2, s)
        end))
chevieset(:G2, :GeneratingRoots, [[1, -1, 0], [-2, 1, 1]])
chevieset(:G2, :HyperplaneRepresentatives, [1, 2])
chevieset(:G2, :Size, 12)
chevieset(:G2, :ReflectionDegrees, [2, 6])
chevieset(:G2, :NrConjugacyClasses, 6)
chevieset(:G2, :CharInfo, function ()
        local res
        res = Dict{Symbol, Any}(:charparams => [[1, 0], [1, 6], [1, 3, 1], [1, 3, 2], [2, 1], [2, 2]], :extRefl => [1, 5, 2], :a => [0, 6, 1, 1, 1, 1], :A => [0, 6, 5, 5, 5, 5])
        res[:b] = map((x->begin
                        x[2]
                    end), res[:charparams])
        res[:B] = [0, 6, 3, 3, 5, 4]
        res[:spaltenstein] = ["1", "\\varepsilon", "\\varepsilon_l", "\\varepsilon_c", "\\theta'", "\\theta''"]
        return res
    end)
chevieset(:G2, :ClassNames, ["A_0", "\\tilde A_1", "A_1", "G_2", "A_2", "A_1+\\tilde A_1"])
chevieset(:G2, :ClassInfo, Dict{Symbol, Any}(:classtext => [[], [2], [1], [1, 2], [1, 2, 1, 2], [1, 2, 1, 2, 1, 2]], :classnames => chevieget(:G2, :ClassNames), :classparams => chevieget(:G2, :ClassNames), :orders => [1, 2, 2, 6, 3, 2], :classes => [1, 3, 3, 2, 2, 1]))
chevieset(:G2, :PowerMaps, [nothing, [1, 1, 1, 5, 5, 1], [1, 2, 3, 6, 1, 6]])
chevieset(:G2, :sparseFakeDegrees, [[1, 0], [1, 6], [1, 3], [1, 3], [1, 1, 1, 5], [1, 2, 1, 4]])
chevieset(:G2, :ClassParameter, (w->begin
            (chevieget(:G2, :ClassNames))[PositionProperty([[[]], [[2], [1, 2, 1], [2, 1, 2, 1, 2]], [[1], [2, 1, 2], [1, 2, 1, 2, 1]], [[2, 1], [1, 2]], [[2, 1, 2, 1], [1, 2, 1, 2]], [[1, 2, 1, 2, 1, 2]]], (x->(w in x;)))]
        end))
chevieset(:G2, :squv, function (para, sqrtpara)
        local u, v
        u = Product(para[1])
        v = Product(para[2])
        if u == v
            return u
        elseif u == v ^ 3
            return -(v ^ 2)
        elseif v == u ^ 3
            return -(u ^ 2)
        elseif sqrtpara[1] !== nothing && sqrtpara[2] !== nothing
            return sqrtpara[1] * sqrtpara[2]
        else
            return GetRoot(u * v, 2, "Hecke(G2)")
        end
    end)
chevieset(:G2, :HeckeCharTable, function (para, sqrtpara)
        local x, y, z, t, tbl, f1, f2, one
        x = (para[1])[1]
        y = (para[1])[2]
        z = (para[2])[1]
        t = (para[2])[2]
        one = (x * y * z * t) ^ 0
        f1(u, v) = begin
                return [1, v, u, v * u, v ^ 2 * u ^ 2, v ^ 3 * u ^ 3] * one
            end
        f2(x, y, z, t, eps) = begin
                local squv
                squv = eps * (chevieget(:G2, :squv))(para, sqrtpara)
                return [2, z + t, x + y, -squv, -x * y * z * t, 2 * squv ^ 3] * one
            end
        tbl = Dict{Symbol, Any}(:identifier => "H(G2)", :parameter => [[x, y], [z, t]], :size => 12, :powermap => chevieget(:G2, :PowerMaps), :irreducibles => [f1(x, z), f1(y, t), f1(y, z), f1(x, t), f2(x, y, z, t, 1), f2(x, y, z, t, -1)], :irredinfo => chevieget(:G2, :IrredInfo))
        Inherit(tbl, chevieget(:G2, :ClassInfo))
        tbl[:centralizers] = map((x->begin
                        tbl[:size] // x
                    end), tbl[:classes])
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        return tbl
    end)
chevieset(:G2, :HeckeRepresentation, function (para, sqrtpara, i)
        local one, squv, x, y, z, t
        one = Product(para[1]) ^ 0 * Product(para[2]) ^ 0
        x = (para[1])[1]
        y = (para[1])[2]
        z = (para[2])[1]
        t = (para[2])[2]
        if i == 1
            return [[[x]], [[z]]] * one
        elseif i == 2
            return [[[y]], [[t]]] * one
        elseif i == 3
            return [[[y]], [[z]]] * one
        elseif i == 4
            return [[[x]], [[t]]] * one
        else
            squv = (chevieget(:G2, :squv))(para, sqrtpara)
            if i == 6
                squv = -squv
            end
            return [[[y, -1], [0, x]], [[z, 0], [squv + y * z + x * t, t]]] * one
        end
    end)
chevieset(:G2, :Representation, function (i,)
        local para
        return (chevieget(:G2, :HeckeRepresentation))([[1, -1], [1, -1]], [1, 1], i)
    end)
chevieset(:G2, :PoincarePolynomial, function (param,)
        local u, v
        u = -((param[1])[1]) // (param[1])[2]
        v = -((param[2])[1]) // (param[2])[2]
        return (1 + u) * (v + 1) * (1 + u * v + u ^ 2 * v ^ 2)
    end)
chevieset(:G2, :SchurModels, Dict{Symbol, Any}(:f1 => Dict{Symbol, Any}(:vcyc => [[[1, -1, 0, 0], 1], [[0, 0, 1, -1], 1], [[1, -1, 1, -1], 3]]), :f2 => Dict{Symbol, Any}(:coeff => -2, :root => [1, -1, 1, -1] // 2, :factor => [-1, 1, 0, 0], :vcyc => [[[0, 0, 0, 0, 1], 3], [[0, 0, -1, 1, 1], 3]])))
chevieset(:G2, :SchurData, [Dict{Symbol, Any}(:name => "f1", :order => [1, 2, 3, 4]), Dict{Symbol, Any}(:name => "f1", :order => [2, 1, 4, 3]), Dict{Symbol, Any}(:name => "f1", :order => [2, 1, 3, 4]), Dict{Symbol, Any}(:name => "f1", :order => [1, 2, 4, 3]), Dict{Symbol, Any}(:name => "f2", :order => [1, 2, 3, 4], :rootPower => -1), Dict{Symbol, Any}(:name => "f2", :order => [1, 2, 3, 4], :rootPower => 1)])
chevieset(:G2, :SchurElement, function (phi, para, sqrtpara)
        local u, v, squv, p
        u = -((para[1])[1]) // (para[1])[2]
        v = -((para[2])[1]) // (para[2])[2]
        p = Position(((chevieget(:G2, :CharInfo))())[:charparams], phi)
        if p == 1
            return (1 + u) * (v + 1) * (u ^ 2 * v ^ 2 + u * v + 1)
        elseif p == 2
            return (((1 + u) * (v + 1) * (u ^ 2 * v ^ 2 + u * v + 1)) // u ^ 3) // v ^ 3
        elseif p == 3
            return ((u ^ 2 + v ^ 2 + u * v) * (1 + u) * (v + 1)) // u ^ 3
        elseif p == 4
            return ((u ^ 2 + v ^ 2 + u * v) * (1 + u) * (v + 1)) // v ^ 3
        end
        squv = ((chevieget(:G2, :squv))(para, sqrtpara) // (para[1])[2]) // (para[2])[2]
        if p == 6
            squv = -squv
        end
        return 2 * (u * v) ^ -1 * (u * v + 1 + squv) * ((u + v) - squv)
    end)
chevieset(:G2, :UnipotentCharacters, function ()
        return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "G", :indices => 1:2, :rank => 2), :levi => [], :parameterExponents => [1, 1], :charNumbers => 1:6, :eigenvalue => 1, :cuspidalName => ""), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [10], :eigenvalue => E(3, 2), :cuspidalName => "G_2[\\zeta_3^2]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [7], :eigenvalue => -1, :cuspidalName => "G_2[-1]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [9], :eigenvalue => E(3), :cuspidalName => "G_2[\\zeta_3]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [8], :eigenvalue => 1, :cuspidalName => "G_2[1]")], :families => [Family("S3", [5, 6, 4, 3, 8, 7, 9, 10], Dict{Symbol, Any}(:ennola => -5)), Family("C1", [1]), Family("C1", [2])], :a => [0, 6, 1, 1, 1, 1, 1, 1, 1, 1], :A => [0, 6, 5, 5, 5, 5, 5, 5, 5, 5], :charSymbols => [[[0], [0], [0], [0], [0], [2]], [[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [1, 2]], [[0], [0], [1], 2, 0], [[0], [0], [1], 2, 1], [[0], [0], [0], [0], [1], [1]], [[0], [0], [0], [1], [0], [1]], [[0, 1], [0], [0, 1], [], [0], []], [[0, 1], [0, 1], [0], [], [], [0]], [[0, 1], [0], [0], [0, 1], [], []], [[0, 1], [0, 1], [], [0], [0], []]])
    end)
chevieset(:G2, :Invariants, [function (x, y)
            return -3 * x * y + 3 * x ^ 2 + y ^ 2
        end, function (x, y)
            return (((x ^ 2 * y ^ 4 - 6 * x ^ 3 * y ^ 3) + 13 * x ^ 4 * y ^ 2) - 12 * x ^ 5 * y) + 4 * x ^ 6
        end])
chevieset(:G2, :Discriminant, function ()
        return function (x, y)
                return 4 * x ^ 3 * y - 27 * y ^ 2
            end
    end)
chevieset(:G2, :UnipotentClasses, function (p, type_)
        local uc, Z, c
        if p == 0
            p = 1
        end
        Z = (n->begin
                    ComplexReflectionGroup(n, 1, 1)
                end)
        uc = Dict{Symbol, Any}(:classes => [Dict{Symbol, Any}(:name => "1", :succ => ["A1"], :dynkin => [0, 0], :balacarter => [], :red => CoxeterGroup("G", 2)), Dict{Symbol, Any}(:name => "A_1", :succ => ["~A1"], :dynkin => [1, 0], :balacarter => [1], :red => Z(2)), Dict{Symbol, Any}(:name => "\\tilde A_1", :succ => ["G2(a1)"], :dynkin => [0, 1], :balacarter => [2], :red => Z(2 - (gcd(p, 3) - 1) // 2)), Dict{Symbol, Any}(:name => "G_2(a_1)", :succ => ["G2"], :dynkin => [2, 0], :balacarter => [1, -2], :Au => CoxeterGroup("A", 2 - (gcd(p, 3) - 1) // 2)), Dict{Symbol, Any}(:name => "G_2", :succ => [], :dynkin => [2, 2], :Au => Z(gcd(p, 6)), :balacarter => [1, 2])], :springerSeries => [Dict{Symbol, Any}(:relgroup => CoxeterGroup("G", 2), :levi => "", :Z => [], :locsys => [[5, 1], [1, 1], [4, 2], [2, 1], [4, 3], [3, 1]]), Dict{Symbol, Any}(:relgroup => CoxeterGroup(), :levi => [1, 2], :Z => [], :locsys => [[4, 1]], :hc => 5)])
        if p == 2
            (((uc[:springerSeries])[1])[:locsys])[1] = [5, 2]
            push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup(), :levi => [1, 2], :Z => [], :locsys => [[5, 1]], :hc => 3))
        elseif p == 3
            push!(uc[:classes], Dict{Symbol, Any}(:name => "(\\tilde A_1)_3", :succ => ["~A1"], :dimBu => 3, :red => Z(2), :Au => CoxeterGroup()))
            push!(((uc[:classes])[1])[:succ], "(~A1)3")
            ((uc[:classes])[3])[:dimBu] = 2
            delete!((uc[:classes])[3], :dynkin)
            (((uc[:springerSeries])[1])[:locsys])[[3, 5]] = [[6, 1], [4, 2]]
            for c = [2, 3]
                push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup(), :levi => [1, 2], :Z => [], :locsys => [[5, c]], :hc => 2c - 2))
            end
        end
        uc[:orderClasses] = map((c->begin
                        map((n->begin
                                    PositionProperty(uc[:classes], (c->begin
                                                (UnipotentClassOps[:Name])(c) == n
                                            end))
                                end), c[:succ])
                    end), uc[:classes])
        for c = uc[:classes]
            delete!(c, :succ)
            if !(haskey(c, :red))
                c[:red] = Z(1)
            end
            if !(haskey(c, :Au))
                c[:Au] = Z(1)
            end
            c[:AuAction] = ExtendedReflectionGroup(c[:red], map((x->begin
                                IdentityMat(rank(c[:red]))
                            end), 1:semisimplerank(c[:Au])))
        end
        return uc
    end)
chevieset(:G2, :KLeftCellRepresentatives, [Dict{Symbol, Any}(:character => [1], :duflo => [1, 2], :reps => ""), Dict{Symbol, Any}(:character => [2], :duflo => [7, 8], :reps => ""), Dict{Symbol, Any}(:character => [3, 5, 6], :duflo => [5, 8], :reps => [[6, 10], [12, 3]]), Dict{Symbol, Any}(:character => [4, 5, 6], :duflo => [7, 3], :reps => [[5, 10], [12, 4]])])