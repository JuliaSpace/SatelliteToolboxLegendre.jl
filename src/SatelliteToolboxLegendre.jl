module SatelliteToolboxLegendre

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./types.jl")

include("./misc.jl")
include("./kernels.jl")
include("./coefficients.jl")

include("./dlegendre/dlegendre.jl")
include("./legendre/legendre.jl")


end # module SatelliteToolboxLegendre
