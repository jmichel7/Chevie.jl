using Test
#using Cyclotomics

@testset "Cyclotomics" begin

    @testset "elementary ops" begin
        @test E(5) isa Number
        @test 0+E(5) isa Cyc{Int}
        @test 0+E(5, 0) isa Cyc{Int}
        @test 0+E(5, 6) isa Cyc{Int}

        E₅ = E(5)
        @test Cyc{Float64}(E₅) isa Cyc{Float64}

        E₅fl = Cyc{Float64}(E₅)

        @test zero(E₅) isa Cyc{Int}
        @test zero(E₅fl) isa Cyc{Float64}
        E₅c=E₅+0
        @test one(E₅c) isa Cyc{Int}
        @test one(E₅fl) isa Cyc{Float64}

        @test deepcopy(E₅c).d !== E₅c.d
        @test deepcopy(E₅c).d == E₅c.d
    end

    @testset "io strings" begin
        @test sprint(print, E(5)) == "E(5)"
        @test sprint(show, 2E(5);context=:limit=>true) == "2ζ₅"
        @test sprint(print, -E(5)) == "-E(5)"
        @test sprint(show, -E(5);context=:limit=>true) == "-ζ₅"
        @test sprint(show, E(1);context=:limit=>true) == "1"

        @test sprint(show, -1.0 * E(5);context=:limit=>true) == "-1.0*ζ₅"
        @test sprint(show, 0.0 * E(4);context=:limit=>true) == "0.0"
        @test sprint(print, 0.0 * E(4)) == "0.0"
        @test sprint(print, 1.0 * E(1)) == " 1.0*E(1)^0"

        using Base.Meta
        x = E(5) + 2E(5)^2
        @test sprint(print, x) == "E(5)+2E(5,2)"
        @test sprint(show, x;context=:limit=>true) == "ζ₅+2ζ₅²"
        @test eval(Base.Meta.parse(sprint(print, x))) == x

        x = E(5) - 2E(5)^2
        @test sprint(print, x) == "E(5)-2E(5,2)"
        @test sprint(show, x;context=:limit=>true) == "ζ₅-2ζ₅²"
        @test eval(Base.Meta.parse(sprint(print, x))) == x

        x = -2E(5) + 2E(5)^2
        @test sprint(print, x) == "-2E(5)+2E(5,2)"
        @test sprint(show, x;context=:limit=>true) == "-2ζ₅+2ζ₅²"
        @test eval(Base.Meta.parse(sprint(print, x))) == x
    end

    @testset "indexing and iteration" begin
        x=Cyc(E(5))
        @test x[0] == 0
        @test x[1] == 1
        @test x[-1] == 0
        @test x[-4] == 1
        @test x[6] == 1
        @test collect(pairs(x)) == [1=>1]
        @test coefficients(x) == [0,1,0,0,0]
    end

    @testset "arithmetic: +, -, module: *, //" begin
        x = E(5)
        y = E(5, 2)

        @test 2x isa Cyc{Int}
        @test 2.0x isa Cyc{Float64}
        @test x * 2.0 isa Cyc{Float64}
        @test div(x, 2) isa Cyc{Int}
        @test x // 2 isa Cyc{Rational{Int}}
        @test x / 2.0 isa Cyc{Float64}
        @test x / 2 isa Cyc{Float64}

        @test x + 2y isa Cyc{Int}
        xy = x + 2y
        @test xy[0] == 0
        @test xy[1] == 1
        @test xy[2] == 2

        @test 2.0xy[0] isa Float64
        @test 2.0xy[1] == 2

        @test (xy-2y)[1] == 1

        @test 1 + x isa Cyc{Int}
        @test (1+x)[0] == 0
        @test x + 1 isa Cyc{Int}
        @test 2.0 + x isa Cyc{Float64}
        @test 2.0 - x isa Cyc{Float64}
        @test x + 2.0 isa Cyc{Float64}
        @test (x+2.0)[0] == 2.0

        # Bug: normalform! is needed in div
        x = E(4, 0) - E(4, 2)
        @test isone(div(x, 2))
        x = 5*E(4,0)
        @test isone(div(x - 2, 3))

        # broadcasting on 1.6 is broken

        @test E(3) .* [1, 2] == [E(3), 2E(3)]
        @test eltype(E(3) .* [1.0, 2.0]) <: Cyc{Float64}
        @test eltype(1 // 1 * E(3) .* [1, 2]) <: Cyc{Rational{Int}}
    end

    @testset "*, powering" begin
        x = E(5)+0
        y = E(5, 2)

        w = x * (x + 2y)
        @test x * y isa Cyc{Int}
        @test (x*y)[1] == 0
        @test (x*y)[2] == 0
        @test (x*y)[3] == 1
        w = (1 + x) * (1 + y)
        @test w[0] == 0
        @test w[1] == 0
        @test w[2] == 0
        @test w[3] == 0

        @test (x^2)[1] == 0
        @test (x^2)[2] == 1

        @test ((1+x)^2)[0] == 0
        @test ((1+x)^2)[1] == 1
        @test ((1+x)^2)[2] == 0
        @test ((1+x)^2)[3] == -1

        @test ((1+x^3)^2)[0] == 0
        @test ((1+x^3)^2)[1] == 0
        @test ((1+x^3)^2)[3] == 1
    end

    @testset "normal form" begin
        x = E(45) + E(45, 5)
        x.d
        @test deepcopy(x) !== x

        y = x

        e = E(45)
        w = e + e^2 + e^8 + e^11 + e^17 + e^26 + e^29 + e^38 + e^44

        @test w == y
        @test y.d !== w.d
        @test y.d == w.d

        @test deepcopy(x) == deepcopy(y)
        @test x !== deepcopy(x)
        @test x.d === copy(x).d

        @test hash(deepcopy(x)) == hash(deepcopy(y))
        @test length(
            Set([deepcopy(x), deepcopy(y), deepcopy(x), deepcopy(y)]),
        ) == 1

        @test iszero(1 + x - x - 1)

        @test isone(sum(-E(5)^i for i in 1:4))
        @test isone(E(5)^5)
        x = E(5)^5
        @test x == sum(-E(5)^i for i in 1:4)
    end

    @testset "predicates" begin
        @test isreal(1 + E(5) - E(5))
        @test isreal(E(5, 1) + E(5, 4))
        @test isreal(E(5, 2) + E(5, 3))
        @test !isreal(E(5, 1) + E(5, 2))
        @test !isreal(E(5, 1) + E(5, 3))

        @test isreal(abs2(E(5, 1) + E(5, 2)))

        @test 1 == E(5)^5
        @test E(5)^2 + E(5)^3 ≈ (-1 - sqrt(5)) / 2
        @test E(5)^2 + E(5)^3 != (-1 - sqrt(5)) / 2
        @test float(E(5)^2 + E(5)^3) == (-1 - sqrt(5)) / 2
        @test 2.0(E(5)^2 + E(5)^3) ≈ (-1 - sqrt(5))
        @test 2.0(E(5)^2 + E(5)^3) != (-1 - sqrt(5))

        @test isapprox(1e-17E(5), 0.0, atol = 1e-12)
        @test isapprox(0.0, 1e-17E(5), atol = 1e-12)
    end

    @testset "embedding" begin
        let x = E(45)^5 + E(45)^10
            @test conductor(x) == 9
            y = x
            @test y == E(9)^2 - E(9)^4 - E(9)^7
            @test conductor(x) == 9
        end

        @test conductor(E(6)) == 3
        @test conductor(E(14)) == 7
        @test conductor(E(1234)) == 617
    end

    @testset "tests against GAP" begin
        @test E(9) == -E(9)^4 - E(9)^7
        @test E(9)^3 == E(3)
        @test E(6) == -E(3)^2
        @test E(12) // 3 == -1 // 3 * E(12)^7

        @test E(45)^4 == -E(45)^19 - E(45)^34
        @test E(45)^13 == -E(45)^28 - E(45)^43
        @test E(45)^14 == -E(45)^29 - E(45)^44
        @test E(45)^22 == -E(45)^7 - E(45)^37

        @test E(5) + E(3) ==
              -E(15)^2 - 2 * E(15)^8 - E(15)^11 - E(15)^13 - E(15)^14
        @test (E(5) + E(5)^4)^2 == -2 * E(5) - E(5)^2 - E(5)^3 - 2 * E(5)^4
        @test E(5) / E(3) == E(15)^13
        @test E(5) * E(3) == E(15)^8
    end

    @testset "conjugation and inverse" begin
        function rand1(α::Cyc, u::AbstractRange, k = 5)
            x = zero(eltype(u)) * α
            for (idx, c) in zip(rand(0:conductor(α), k), rand(u, k))
              x+=c*E(conductor(α),idx)
            end
            return x
        end

        for x in [E(45) + E(45)^2, E(45) + E(45)^2 // 1]
            @test galois(x, 1) == x
            y=prod(galois(x,i) for i in 2:conductor(x)if gcd(conductor(x),i)==1)
            @test isreal(x * y)
            @test y ==
                  E(45)^2 + E(45)^3 - E(45)^6 - E(45)^8 + E(45)^11 - E(45)^12 -
                  2 * E(45)^16 +
                  E(45)^17 +
                  E(45)^19 +
                  E(45)^21 - 2 * E(45)^24 - E(45)^26 - E(45)^28 + 2 * E(45)^29 -
                  E(45)^34 + E(45)^37 - 2 * E(45)^42 - E(45)^43 + E(45)^44

            @test galois(x, 1) == x
            @test_throws DomainError galois(x, 5)
        end

        for x in [
            E(45)^5 + E(45)^10,
            E(45) - E(45)^5,
            rand1(E(45)+0, -5:5, 3),
            rand1(E(45)+0, -1:1, 5),
        ]
            iszero(x) || @test isone(x * inv(x // big(1)))
        end

        for x in
            [E(45)^5 + big(1) * E(45)^10, (E(45)^5) // 1 + big(1) * E(45)^10]
            y = x
            @test inv(y) == -E(9)^2 + E(9)^3 - E(9)^4 == inv(x)
            @test inv(y) * x == inv(x) * x == one(x)
        end

        for x in
            [E(45)^5 + big(1) * E(45)^10, (E(45)^5) // 1 + big(1) * E(45)^10]
            y = x
            @test inv(y) == -E(9)^2 + E(9)^3 - E(9)^4 == inv(x)
            @test inv(y) * x == inv(x) * x == one(x)
        end

        for x in [0.5 + 0.75 * E(4), 1 // 2 + 3 // 4 * E(4)]
            @test real(x) == 0.5
            @test imag(x) == 0.75
            @test float(real(x)) isa Float64
            @test float(imag(x)) isa Float64
            @test x == real(x) + im * imag(x)
            @test_throws InexactError float(x)
            @test_throws InexactError Rational{Int}(x)
            @test Rational{Int}(x + conj(x)) isa Rational{Int}
            @test Rational{Int}(x + conj(x)) == x + conj(x)
            if valtype(x) <: Rational
                @test Rational(x + conj(x)) isa Rational{Int}
                @test Rational(x + conj(x)) == 1 // 1
            else
                @test_throws MethodError Rational(x + conj(x))
            end
        end

        let x = E(45)^5 + E(45)^10
            @test isreal(x + conj(x))
            @test float(x + conj(x)) isa Float64
            @test_throws InexactError Float64(x)
            @test_throws InexactError Rational(x)
        end

    end

    @testset "Conversions" begin
        C = typeof(1+E(3))

        @test iszero(C(0))
        @test isone(C(1))
        @test isone(C(1.0))
        @test zeros(typeof(E(3)+0), 2, 2) isa Matrix{<:Cyc}

        v = [Cyc(E(3)^i) for i in 1:3]

        @test (v[1] = 1.0 * E(5)) isa Cyc{Float64}
        @test eltype(v) <: Cyc{Int}
        @test v[1] isa Cyc{Int}
        @test (v[1] = 2.0 * E(5)) isa Cyc{Float64}
        @test v[1] == 2E(5)
        @test_throws InexactError v[1] = 2.5 * E(5)
    end

    @testset "Conversions to julia types" begin
        x = E(3)
        y = x + x^2

        @test Int(y) == -1
        @test y == -1
        @test_throws InexactError Int(x)

        @test float(y) == -1.0
        @test y == -1.0
        @test_throws InexactError float(x)

        @test Float64(y) == -1.0
        # @test y == -1.0
        @test_throws InexactError Float64(x)

        @test float(y // big(1)) isa BigFloat
        @test Float64(y // big(1)) isa Float64

        @test ComplexF64(x) isa ComplexF64
        @test Complex{BigFloat}(x) isa Complex{BigFloat}
        γ = 2 * π / 3
        bγ = 2 * big(π) / 3
        @test Complex{BigFloat}(x) ≈ cos(bγ) + im * sin(bγ)

        @test ComplexF64(x)^2 ≈ ComplexF64(x^2)

        @test complex(x) isa ComplexF64
        @test complex(x) == cos(γ) + im * sin(γ)

        @test complex(big(1) * x) isa Complex{BigFloat}
        @test complex(big(1) * x) ≈ cos(bγ) + im * sin(bγ)

        @test abs(x) ≈ 1.0
        @test abs(big(1) / 3 * x) ≈ big(1) / 3

        @test ComplexF64 == typeof(@inferred ComplexF64(x))
        @test ComplexF64 == typeof(@inferred ComplexF64(y))
        @test ComplexF64 == typeof(@inferred ComplexF64(big(1) * x))
        @test ComplexF64 == typeof(@inferred ComplexF64(big(1) * y))

        @test_logs (
            :error,
            "The cyclotomic is real but it can not be converted to Rational:  1*E(5)^2 + 1*E(5)^3 ≈ -1.618033988749895",
        ) try
            Rational(E(5)^2 + E(5)^3)
        catch ex
            @test ex isa InexactError
        end
    end

    @testset "Complex arithmetic" begin
        @test real(E(4)) == 0
        @test imag(-E(4)) == -1
        @test Cyc(1 + 2im) isa Cyc
        @test Cyc(1.0 - 2.0im) isa Cyc
        @test reim(Cyc(2 - 3im)) == (2, -3)

        @test E(4) * im == -1
        @test E(4) - im == 0

        @test E(9) + im // 3 isa Cyc
    end

end
