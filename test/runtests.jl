using Test

using SatelliteToolboxLegendre

@testset "Legendre Associated Functions" verbose = true begin
    include("./legendre.jl")
end

@testset "Derivative of the Legendre Associated Functions" verbose = true begin
    include("./dlegendre.jl")
end
