chevieset(["G24", "G25", "G26", "G27", "G29", "G31", "G32", "G33", "G34", "H3", "H4", "2E6", "2F4", "3D4", "E6", "E7", "E8", "F4", "G2"], :IrredInfo, function (t,)
        local ci
        ci = (chevieget(t, :CharInfo))()
        return map(function (x, y)
                    return Dict{Symbol, Any}(:charparam => x, :charname => y)
                end, ci[:charparams], ci[:charnames])
    end)
chevieset(["G24", "G25", "G27", "G29", "G31", "G32", "G33", "G34", "E6", "E7", "E8", "2E6", "2F4", "3D4", "H3", "H4"], :ReflectionName, (t->begin
            function (option,)
                local i, o
                i = ["G24", "G25", "G27", "G29", "G31", "G32", "G33", "G34", "E6", "E7", "E8", "2E6", "2F4", "3D4", "H3", "H4"]
                o = ["G_{24}", "G_{25}", "G_{27}", "G_{29}", "G_{31}", "G_{32}", "G_{33}", "G_{34}", "E_6", "E_7", "E_8", "{}^2E_6", "{}^2F_4", "{}^3D_4", "H_3", "H_4"]
                if haskey(option, :TeX)
                    return o[Position(i, t)]
                else
                    return t
                end
            end
        end))
chevieset(["D", "2A", "2D"], :ReflectionName, (t->begin
            function (r, option)
                local i, o
                i = ["D", "2A", "2D"]
                o = ["D", "{}^2A", "{}^2D"]
                if haskey(option, :arg)
                    return SPrint(FormatGAP(t), ",", r)
                elseif haskey(option, :TeX)
                    return SPrint(o[Position(i, t)], "_", TeXBracket(r))
                else
                    return SPrint(t, r)
                end
            end
        end))
chevieset(["3D4", "E6", "2E6", "E7", "E8", "F4", "2F4", "G2", "H3", "H4"], :CharTable, (t->begin
            function ()
                local res, rank
                rank = Position("12345678", t[length(t)])
                res = (chevieget(t, :HeckeCharTable))(map((x->begin
                                    [1, -1]
                                end), 1:rank), map((x->begin
                                    1
                                end), 1:rank))
                ((CHEVIE[:compat])[:ChangeIdentifier])(res, SPrint("W(", t, ")"))
                return res
            end
        end))
chevieset(["G24", "G27", "G29", "G33", "G34", "H3", "H4", "E6", "E7", "E8"], :PoincarePolynomial, (t->begin
            function (q,)
                return Product(chevieget(t, :ReflectionDegrees), (x->begin
                                Sum(0:x - 1, (y->begin
                                            (-((q[1])[1]) // (q[1])[2]) ^ y
                                        end))
                            end))
            end
        end))
chevieset(["G24", "G25", "G26", "G27", "G29"], :Representation, (t->begin
            function (i,)
                local para
                para = chevieget(t, :EigenvaluesGeneratingReflections)
                para = map((x->begin
                                map((j->begin
                                            E(1 // x, j)
                                        end), 0:1 // x - 1)
                            end), para)
                return (chevieget(t, :HeckeRepresentation))(para, [], i)
            end
        end))
chevieset(["G2", "F4", "H3", "E6", "G24", "G25", "G26", "G27", "G29", "G31", "G32", "G33", "G34"], :SemisimpleRank, function (t,)
        local r
        r = chevieget(t, :GeneratingRoots)
        if r isa Function
            r = r()
        end
        return length(r[1])
    end)
chevieset(["A", "B", "D"], :SemisimpleRank, (t->begin
            r->begin
                    r
                end
        end))
chevieset(["3D4", "G2", "F4", "2F4", "H3", "E6", "G24", "G25", "G26", "G27", "G29", "G32", "G33", "G34"], :FakeDegree, (t->begin
            function (phi, q)
                local f
                f = (chevieget(t, :sparseFakeDegrees))[Position(((chevieget(t, :CharInfo))())[:charparams], phi)]
                return Sum(1:3 - 1:length(f) - 1, (i->begin
                                f[i] * q ^ f[i + 1]
                            end))
            end
        end))
chevieset(["H4", "E7", "E8", "G31"], :FakeDegree, (t->begin
            function (phi, q)
                local f, res
                f = (chevieget(t, :cycpolfakedegrees))[Position(((chevieget(t, :CharInfo))())[:charparams], phi)]
                if IsList(f[1])
                    res = ValuePol(f[1], q ^ 2)
                else
                    res = f[1]
                end
                f = copy(f)
                f[1] = 1
                return res * Value(CycPol(f), q)
            end
        end))
chevieset(["H4", "E7", "E8", "G31"], :HighestPowerFakeDegrees, (t->begin
            function ()
                return map(function (f,)
                            local res
                            if IsList(f[1])
                                res = (2 * length(f[1]) + f[2]) - 2
                            else
                                res = f[2]
                            end
                            return res + Sum(f[3:length(f)], phi)
                        end, chevieget(t, :cycpolfakedegrees))
            end
        end))
chevieset(["E6", "G32", "G33", "G34", "G2", "F4", "H3", "G24", "G25", "G26", "G27", "G29"], :HighestPowerFakeDegrees, (t->begin
            function ()
                return map((x->begin
                                x[length(x)]
                            end), chevieget(t, :sparseFakeDegrees))
            end
        end))
chevieset(["G2", "F4", "H3", "H4", "G24", "G25", "G26", "G27", "G29", "E6", "E7", "E8", "G31", "G32", "G33", "G34"], :LowestPowerFakeDegrees, (t->begin
            function ()
                return map((x->begin
                                x[2]
                            end), chevieget(t, :sparseFakeDegrees))
            end
        end))
chevieset(["G24", "G27", "G29", "G33", "G34", "H3", "H4", "E6", "E7", "E8"], :HighestPowerGenericDegrees, (t->begin
            function ()
                local N
                N = Sum(chevieget(t, :ReflectionDegrees), (x->begin
                                x - 1
                            end))
                return map((x->begin
                                N - degree(CycPol(x))
                            end), chevieget(t, :CycPolSchurElements))
            end
        end))
chevieset(["G24", "G27", "G29", "G33", "G34", "H3", "H4", "E6", "E7", "E8"], :LowestPowerGenericDegrees, (t->begin
            function ()
                return map((x->begin
                                -(x[2])
                            end), chevieget(t, :CycPolSchurElements))
            end
        end))
chevieset(["F4", "G2", "G25", "G26"], :DecompositionMatrix, (t->begin
            function (p,)
                local T, m
                T = (chevieget(t, :CharTable))()
                T[:name] = T[:identifier]
                m = DecompositionMatrix(mod(T, p))
                return map((c->begin
                                [c[1], (m[c[1]])[c[2]]]
                            end), BlocksMat(m))
            end
        end))
chevieset(["G24", "G27", "G29", "G33", "G34", "E6", "E7", "E8", "H3", "H4"], :SchurElement, (t->begin
            function (arg...,)
                return Value(CycPol((chevieget(t, :CycPolSchurElements))[Position(((chevieget(t, :CharInfo))())[:charparams], arg[1])]), -(((arg[2])[1])[1]) // ((arg[2])[1])[2])
            end
        end))
chevieset(["G2", "F4", "G25", "G26", "G32"], :FactorizedSchurElement, (t->begin
            function (arg...,)
                local Y, ci
                Y = Concatenation((arg[2])[chevieget(t, :HyperplaneRepresentatives)])
                ci = (chevieget(t, :SchurData))[Position(((chevieget(t, :CharInfo))())[:charparams], arg[1])]
                return ApplyFunc(VFactorSchurElement, Concatenation([Y, (chevieget(t, :SchurModels))[Symbol(ci[:name])], ci], arg[3:length(arg)]))
            end
        end))
chevieset(["F4", "G25", "G26", "G32"], :SchurElement, (t->begin
            function (arg...,)
                local Y, ci
                Y = Concatenation((arg[2])[chevieget(t, :HyperplaneRepresentatives)])
                ci = (chevieget(t, :SchurData))[Position(((chevieget(t, :CharInfo))())[:charparams], arg[1])]
                return VcycSchurElement(Y, (chevieget(t, :SchurModels))[Symbol(ci[:name])], ci)
            end
        end))
chevieset(["E7", "E8", "F4", "2F4", "G2", "3D4", "H3", "H4", "G24", "G25", "G26", "G27", "G29", "G32", "G33", "G34"], :Ennola, (t->begin
            function (arg...,)
                local uc, res, p, A, b, f
                uc = chevieget(t, :UnipotentCharacters)
                if uc isa Function
                    uc = uc()
                end
                res = uc[:a] * 0
                for f = uc[:families]
                    if haskey(f, :ennola)
                        if IsList(f[:ennola])
                            p = SPerm(f[:ennola])
                        else
                            A = fusion_algebra(f)
                            b = basis(A)
                            if !(haskey(f, :ennola))
                                f[:ennola] = f[:special]
                            end
                            if f[:ennola] > 0
                                p = SPerm(b, b[f[:ennola]] * b)
                            else
                                p = SPerm(b, -(b[-(f[:ennola])]) * b)
                            end
                        end
                    else
                        p = SPerm()
                    end
                    res[f[:charNumbers]] = Permuted(f[:charNumbers], p)
                end
                return SPerm(res)
            end
        end))
