
chevieset(:G2, :CartanMat, function (arg...,)
        #= none:9 =#
        local a, type_
        #= none:11 =#
        if length(arg) > 0
            #= none:11 =#
            type_ = arg[1]
        else
            #= none:12 =#
            type_ = 1
        end
        #= none:14 =#
        a = [[2, -1], [-3, 2]]
        #= none:16 =#
        (a[1])[2] = -type_
        #= none:17 =#
        (a[2])[1] = 3 // (a[1])[2]
        #= none:19 =#
        return a
    end)
chevieset(:G2, :PrintDiagram, function (indices, title, type_)
        #= none:4 =#
        print(title, " ", indices[1])
        #= none:6 =#
        if type_ == 1
            #= none:6 =#
            print(" >>> ")
        elseif #= none:8 =# type_ == ER(3)
            #= none:8 =#
            print(" ==6== ")
        else
            #= none:10 =#
            print(" ?6? ")
        end
        #= none:12 =#
        print(indices[2], " \n")
    end)
chevieset(:G2, :ReflectionName, function (arg...,)
        #= none:3 =#
        local i, opt, type_
        #= none:5 =#
        if length(arg) == 1
            #= none:5 =#
            return "G2(?)"
        end
        #= none:7 =#
        type_ = arg[2]
        #= none:8 =#
        opt = arg[1]
        #= none:10 =#
        if type_ == 1
            #= none:11 =#
            if haskey(opt, :TeX)
                #= none:11 =#
                return "G_2"
            elseif #= none:13 =# haskey(opt, :arg)
                #= none:13 =#
                return "\"G\",2"
            else
                #= none:15 =#
                return "G2"
            end
        elseif #= none:17 =# type_ == ER(3)
            #= none:18 =#
            if haskey(opt, :TeX)
                #= none:18 =#
                return "G_{\\hbox{sym}2}"
            elseif #= none:20 =# haskey(opt, :arg)
                #= none:20 =#
                return "\"Gsym\",2"
            else
                #= none:22 =#
                return "Gsym2"
            end
        elseif #= none:24 =# haskey(opt, :TeX)
            #= none:24 =#
            return SPrint("G_2(", Format(type_ ^ 2 // 3, opt), ")")
        elseif #= none:26 =# haskey(opt, :arg)
            #= none:26 =#
            return SPrint("\"G\",", 2, ",", Format(type_ ^ 2 // 3, opt))
        else
            #= none:28 =#
            return SPrint("G2(", Format(type_ ^ 2 // 3, opt), ")")
        end
    end)
chevieset(:G2, :ParabolicRepresentatives, (s->begin
            #= none:4 =#
            (chevieget(:imp, :ParabolicRepresentatives))(6, 6, 2, s)
        end))
chevieset(:G2, :GeneratingRoots, [[1, -1, 0], [-2, 1, 1]])
chevieset(:G2, :HyperplaneRepresentatives, [1, 2])
chevieset(:G2, :Size, 12)
chevieset(:G2, :ReflectionDegrees, [2, 6])
chevieset(:G2, :NrConjugacyClasses, 6)
chevieset(:G2, :CharInfo, function ()
        #= none:3 =#
        local res
        #= none:5 =#
        res = Dict{Symbol, Any}(:charparams => [[1, 0], [1, 6], [1, 3, 1], [1, 3, 2], [2, 1], [2, 2]], :extRefl => [1, 5, 2], :a => [0, 6, 1, 1, 1, 1], :A => [0, 6, 5, 5, 5, 5])
        #= none:10 =#
        res[:b] = map((x->begin
                        #= none:10 =#
                        x[2]
                    end), res[:charparams])
        #= none:12 =#
        res[:B] = [0, 6, 3, 3, 5, 4]
        #= none:15 =#
        res[:spaltenstein] = ["1", "\\varepsilon", "\\varepsilon_l", "\\varepsilon_c", "\\theta'", "\\theta''"]
        #= none:18 =#
        return res
    end)
chevieset(:G2, :ClassNames, ["A_0", "\\tilde A_1", "A_1", "G_2", "A_2", "A_1+\\tilde A_1"])
chevieset(:G2, :ClassInfo, Dict{Symbol, Any}(:classtext => [[], [2], [1], [1, 2], [1, 2, 1, 2], [1, 2, 1, 2, 1, 2]], :classnames => chevieget(:G2, :ClassNames), :classparams => chevieget(:G2, :ClassNames), :orders => [1, 2, 2, 6, 3, 2], :classes => [1, 3, 3, 2, 2, 1]))
chevieset(:G2, :PowerMaps, [nothing, [1, 1, 1, 5, 5, 1], [1, 2, 3, 6, 1, 6]])
chevieset(:G2, :sparseFakeDegrees, [[1, 0], [1, 6], [1, 3], [1, 3], [1, 1, 1, 5], [1, 2, 1, 4]])
chevieset(:G2, :ClassParameter, (w->begin
            #= none:9 =#
            (chevieget(:G2, :ClassNames))[PositionProperty([[[]], [[2], [1, 2, 1], [2, 1, 2, 1, 2]], [[1], [2, 1, 2], [1, 2, 1, 2, 1]], [[2, 1], [1, 2]], [[2, 1, 2, 1], [1, 2, 1, 2]], [[1, 2, 1, 2, 1, 2]]], (x->(#= none:10 =#
                            w in x)))]
        end))
chevieset(:G2, :squv, function (para, sqrtpara)
        #= none:10 =#
        local u, v
        #= none:12 =#
        u = Product(para[1])
        #= none:13 =#
        v = Product(para[2])
        #= none:15 =#
        if u == v
            #= none:15 =#
            return u
        elseif #= none:17 =# u == v ^ 3
            #= none:17 =#
            return -(v ^ 2)
        elseif #= none:19 =# v == u ^ 3
            #= none:19 =#
            return -(u ^ 2)
        elseif #= none:21 =# sqrtpara[1] !== nothing && sqrtpara[2] !== nothing
            #= none:22 =#
            return sqrtpara[1] * sqrtpara[2]
        else
            #= none:24 =#
            return GetRoot(u * v, 2, "Hecke(G2)")
        end
    end)
chevieset(:G2, :HeckeCharTable, function (para, sqrtpara)
        #= none:4 =#
        local x, y, z, t, tbl, f1, f2, one
        #= none:6 =#
        x = (para[1])[1]
        #= none:7 =#
        y = (para[1])[2]
        #= none:8 =#
        z = (para[2])[1]
        #= none:9 =#
        t = (para[2])[2]
        #= none:11 =#
        one = (x * y * z * t) ^ 0
        #= none:13 =#
        f1 = function (u, v)
                #= none:13 =#
                return [1, v, u, v * u, v ^ 2 * u ^ 2, v ^ 3 * u ^ 3] * one
            end
        #= none:16 =#
        f2 = function (x, y, z, t, eps)
                #= none:16 =#
                local squv
                #= none:18 =#
                squv = eps * (chevieget(:G2, :squv))(para, sqrtpara)
                #= none:20 =#
                return [2, z + t, x + y, -squv, -x * y * z * t, 2 * squv ^ 3] * one
            end
        #= none:23 =#
        tbl = Dict{Symbol, Any}(:identifier => "H(G2)", :parameter => [[x, y], [z, t]], :size => 12, :powermap => chevieget(:G2, :PowerMaps), :irreducibles => [f1(x, z), f1(y, t), f1(y, z), f1(x, t), f2(x, y, z, t, 1), f2(x, y, z, t, -1)], :irredinfo => chevieget(:G2, :IrredInfo))
        #= none:29 =#
        Inherit(tbl, chevieget(:G2, :ClassInfo))
        #= none:31 =#
        tbl[:centralizers] = map((x->begin
                        #= none:31 =#
                        tbl[:size] // x
                    end), tbl[:classes])
        #= none:33 =#
        tbl = ((CHEVIE[:compat])[:MakeCharacterTable])(tbl)
        #= none:35 =#
        return tbl
    end)
chevieset(:G2, :HeckeRepresentation, function (para, sqrtpara, i)
        #= none:4 =#
        local one, squv, x, y, z, t
        #= none:5 =#
        one = Product(para[1]) ^ 0 * Product(para[2]) ^ 0
        #= none:7 =#
        x = (para[1])[1]
        #= none:8 =#
        y = (para[1])[2]
        #= none:9 =#
        z = (para[2])[1]
        #= none:10 =#
        t = (para[2])[2]
        #= none:12 =#
        if i == 1
            #= none:12 =#
            return [[[x]], [[z]]] * one
        elseif #= none:14 =# i == 2
            #= none:14 =#
            return [[[y]], [[t]]] * one
        elseif #= none:16 =# i == 3
            #= none:16 =#
            return [[[y]], [[z]]] * one
        elseif #= none:18 =# i == 4
            #= none:18 =#
            return [[[x]], [[t]]] * one
        else
            #= none:21 =#
            squv = (chevieget(:G2, :squv))(para, sqrtpara)
            #= none:23 =#
            if i == 6
                #= none:23 =#
                squv = -squv
            end
            #= none:25 =#
            return [[[y, -1], [0, x]], [[z, 0], [squv + y * z + x * t, t]]] * one
        end
    end)
chevieset(:G2, :Representation, function (i,)
        #= none:3 =#
        local para
        #= none:5 =#
        return (chevieget(:G2, :HeckeRepresentation))([[1, -1], [1, -1]], [1, 1], i)
    end)
chevieset(:G2, :PoincarePolynomial, function (param,)
        #= none:9 =#
        local u, v
        #= none:11 =#
        u = -((param[1])[1]) // (param[1])[2]
        #= none:12 =#
        v = -((param[2])[1]) // (param[2])[2]
        #= none:14 =#
        return (1 + u) * (v + 1) * (1 + u * v + u ^ 2 * v ^ 2)
    end)
chevieset(:G2, :SchurModels, Dict{Symbol, Any}(:f1 => Dict{Symbol, Any}(:vcyc => [[[1, -1, 0, 0], 1], [[0, 0, 1, -1], 1], [[1, -1, 1, -1], 3]]), :f2 => Dict{Symbol, Any}(:coeff => -2, :root => [1, -1, 1, -1] // 2, :factor => [-1, 1, 0, 0], :vcyc => [[[0, 0, 0, 0, 1], 3], [[0, 0, -1, 1, 1], 3]])))
chevieset(:G2, :SchurData, [Dict{Symbol, Any}(:name => "f1", :order => [1, 2, 3, 4]), Dict{Symbol, Any}(:name => "f1", :order => [2, 1, 4, 3]), Dict{Symbol, Any}(:name => "f1", :order => [2, 1, 3, 4]), Dict{Symbol, Any}(:name => "f1", :order => [1, 2, 4, 3]), Dict{Symbol, Any}(:name => "f2", :order => [1, 2, 3, 4], :rootPower => -1), Dict{Symbol, Any}(:name => "f2", :order => [1, 2, 3, 4], :rootPower => 1)])
chevieset(:G2, :SchurElement, function (phi, para, sqrtpara)
        #= none:7 =#
        local u, v, squv, p
        #= none:9 =#
        u = -((para[1])[1]) // (para[1])[2]
        #= none:10 =#
        v = -((para[2])[1]) // (para[2])[2]
        #= none:12 =#
        p = Position(((chevieget(:G2, :CharInfo))())[:charparams], phi)
        #= none:14 =#
        if p == 1
            #= none:14 =#
            return (1 + u) * (v + 1) * (u ^ 2 * v ^ 2 + u * v + 1)
        elseif #= none:16 =# p == 2
            #= none:16 =#
            return (((1 + u) * (v + 1) * (u ^ 2 * v ^ 2 + u * v + 1)) // u ^ 3) // v ^ 3
        elseif #= none:18 =# p == 3
            #= none:18 =#
            return ((u ^ 2 + v ^ 2 + u * v) * (1 + u) * (v + 1)) // u ^ 3
        elseif #= none:20 =# p == 4
            #= none:20 =#
            return ((u ^ 2 + v ^ 2 + u * v) * (1 + u) * (v + 1)) // v ^ 3
        end
        #= none:22 =#
        squv = ((chevieget(:G2, :squv))(para, sqrtpara) // (para[1])[2]) // (para[2])[2]
        #= none:24 =#
        if p == 6
            #= none:24 =#
            squv = -squv
        end
        #= none:26 =#
        return 2 * (u * v) ^ -1 * (u * v + 1 + squv) * ((u + v) - squv)
    end)
chevieset(:G2, :UnipotentCharacters, function ()
        #= none:4 =#
        return Dict{Symbol, Any}(:harishChandra => [Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "G", :indices => 1:2, :rank => 2), :levi => [], :parameterExponents => [1, 1], :charNumbers => 1:6, :eigenvalue => 1, :cuspidalName => ""), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [10], :eigenvalue => E(3, 2), :cuspidalName => "G_2[\\zeta_3^2]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [7], :eigenvalue => -1, :cuspidalName => "G_2[-1]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [9], :eigenvalue => E(3), :cuspidalName => "G_2[\\zeta_3]"), Dict{Symbol, Any}(:relativeType => Dict{Symbol, Any}(:series => "A", :indices => [], :rank => 0), :levi => 1:2, :parameterExponents => [], :charNumbers => [8], :eigenvalue => 1, :cuspidalName => "G_2[1]")], :families => [Family("S3", [5, 6, 4, 3, 8, 7, 9, 10], Dict{Symbol, Any}(:ennola => -5)), Family("C1", [1]), Family("C1", [2])], :a => [0, 6, 1, 1, 1, 1, 1, 1, 1, 1], :A => [0, 6, 5, 5, 5, 5, 5, 5, 5, 5], :charSymbols => [[[0], [0], [0], [0], [0], [2]], [[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [1, 2]], [[0], [0], [1], 2, 0], [[0], [0], [1], 2, 1], [[0], [0], [0], [0], [1], [1]], [[0], [0], [0], [1], [0], [1]], [[0, 1], [0], [0, 1], [], [0], []], [[0, 1], [0, 1], [0], [], [], [0]], [[0, 1], [0], [0], [0, 1], [], []], [[0, 1], [0, 1], [], [0], [0], []]])
    end)
chevieset(:G2, :Invariants, [function (x, y)
            #= none:4 =#
            return -3 * x * y + 3 * x ^ 2 + y ^ 2
        end, function (x, y)
            #= none:6 =#
            return (((x ^ 2 * y ^ 4 - 6 * x ^ 3 * y ^ 3) + 13 * x ^ 4 * y ^ 2) - 12 * x ^ 5 * y) + 4 * x ^ 6
        end])
chevieset(:G2, :Discriminant, function ()
        #= none:3 =#
        return function (x, y)
                #= none:3 =#
                return 4 * x ^ 3 * y - 27 * y ^ 2
            end
    end)
chevieset(:G2, :UnipotentClasses, function (p, type_)
        #= none:4 =#
        local uc, Z, c
        #= none:6 =#
        if p == 0
            #= none:6 =#
            p = 1
        end
        #= none:7 =#
        Z = (n->begin
                    #= none:7 =#
                    ComplexReflectionGroup(n, 1, 1)
                end)
        #= none:9 =#
        uc = Dict{Symbol, Any}(:classes => [Dict{Symbol, Any}(:name => "1", :succ => ["A1"], :dynkin => [0, 0], :balacarter => [], :red => CoxeterGroup("G", 2)), Dict{Symbol, Any}(:name => "A_1", :succ => ["~A1"], :dynkin => [1, 0], :balacarter => [1], :red => Z(2)), Dict{Symbol, Any}(:name => "\\tilde A_1", :succ => ["G2(a1)"], :dynkin => [0, 1], :balacarter => [2], :red => Z(2 - (gcd(p, 3) - 1) // 2)), Dict{Symbol, Any}(:name => "G_2(a_1)", :succ => ["G2"], :dynkin => [2, 0], :balacarter => [1, -2], :Au => CoxeterGroup("A", 2 - (gcd(p, 3) - 1) // 2)), Dict{Symbol, Any}(:name => "G_2", :succ => [], :dynkin => [2, 2], :Au => Z(gcd(p, 6)), :balacarter => [1, 2])], :springerSeries => [Dict{Symbol, Any}(:relgroup => CoxeterGroup("G", 2), :levi => "", :Z => [], :locsys => [[5, 1], [1, 1], [4, 2], [2, 1], [4, 3], [3, 1]]), Dict{Symbol, Any}(:relgroup => CoxeterGroup(), :levi => [1, 2], :Z => [], :locsys => [[4, 1]], :parameter => [8])])
        #= none:23 =#
        if p == 2
            #= none:23 =#
            (((uc[:springerSeries])[1])[:locsys])[1] = [5, 2]
            #= none:25 =#
            push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup(), :levi => [1, 2], :Z => [], :locsys => [[5, 1]]))
        elseif #= none:28 =# p == 3
            #= none:29 =#
            push!(uc[:classes], Dict{Symbol, Any}(:name => "(\\tilde A_1)_3", :succ => ["~A1"], :dimBu => 3, :red => Z(2), :Au => CoxeterGroup()))
            #= none:32 =#
            push!(((uc[:classes])[1])[:succ], "(~A1)3")
            #= none:34 =#
            ((uc[:classes])[3])[:dimBu] = 2
            #= none:35 =#
            delete!((uc[:classes])[3], :dynkin)
            #= none:37 =#
            (((uc[:springerSeries])[1])[:locsys])[[3, 5]] = [[6, 1], [4, 2]]
            #= none:39 =#
            for c = [2, 3]
                #= none:39 =#
                push!(uc[:springerSeries], Dict{Symbol, Any}(:relgroup => CoxeterGroup(), :levi => [1, 2], :Z => [], :locsys => [[5, c]]))
            end
        end
        #= none:43 =#
        uc[:orderClasses] = map((c->begin
                        #= none:43 =#
                        map((n->begin
                                    #= none:44 =#
                                    PositionProperty(uc[:classes], (c->begin
                                                #= none:44 =#
                                                (UnipotentClassOps[:Name])(c) == n
                                            end))
                                end), c[:succ])
                    end), uc[:classes])
        #= none:46 =#
        for c = uc[:classes]
            #= none:46 =#
            delete!(c, :succ)
            #= none:48 =#
            if !(haskey(c, :red))
                #= none:48 =#
                c[:red] = Z(1)
            end
            #= none:50 =#
            if !(haskey(c, :Au))
                #= none:50 =#
                c[:Au] = Z(1)
            end
            #= none:52 =#
            c[:AuAction] = ExtendedReflectionGroup(c[:red], map((x->begin
                                #= none:53 =#
                                IdentityMat(Rank(c[:red]))
                            end), 1:SemisimpleRank(c[:Au])))
        end
        #= none:55 =#
        return uc
    end)
chevieset(:G2, :KLeftCellRepresentatives, [Dict{Symbol, Any}(:character => [1], :duflo => [1, 2], :reps => ""), Dict{Symbol, Any}(:character => [2], :duflo => [7, 8], :reps => ""), Dict{Symbol, Any}(:character => [3, 5, 6], :duflo => [5, 8], :reps => [[6, 10], [12, 3]]), Dict{Symbol, Any}(:character => [4, 5, 6], :duflo => [7, 3], :reps => [[5, 10], [12, 4]])])
