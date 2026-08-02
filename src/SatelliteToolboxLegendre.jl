module SatelliteToolboxLegendre

using StyledStrings: @styled_str

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./types.jl")

include("./misc.jl")
include("./kernels.jl")
include("./coefficients.jl")

include("./dlegendre.jl")
include("./legendre.jl")

include("./show.jl")

end # module SatelliteToolboxLegendre
