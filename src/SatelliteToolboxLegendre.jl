module SatelliteToolboxLegendre

############################################################################################
#                                         Includes
############################################################################################

include("./misc.jl")

include("./dlegendre/fully_normalized_dlegendre.jl")
include("./dlegendre/schmidt_quasi_normalized_dlegendre.jl")
include("./dlegendre/unnormalized_dlegendre.jl")
include("./dlegendre/dlegendre.jl")

include("./legendre/fully_normalized_legendre.jl")
include("./legendre/schmidt_quasi_normalized_legendre.jl")
include("./legendre/unnormalized_legendre.jl")
include("./legendre/legendre.jl")


end # module SatelliteToolboxLegendre
