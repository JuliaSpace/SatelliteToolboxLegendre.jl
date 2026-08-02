using Test

using OffsetArrays
using SatelliteToolboxLegendre

@testset "Legendre Associated Functions" verbose = true begin
    include("./legendre.jl")
end

@testset "Derivative of the Legendre Associated Functions" verbose = true begin
    include("./dlegendre.jl")
end

@testset "Matrices With Offset Axes" verbose = true begin
    include("./offset_arrays.jl")
end

@testset "Precomputed Coefficients" verbose = true begin
    include("./coefficients.jl")
end
