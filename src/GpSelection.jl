module GpSelection
    export fixMatrix, calcSigma, generateZ, generateData, Options, PriorPars

    using Distributions
    using MatrixDepot
    using Combinatorics

    include("utils.jl")
    include("sampledata.jl")
    include("varBayes.jl")
end
