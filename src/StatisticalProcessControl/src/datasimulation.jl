using StatisticalProcessControl

@with_kw struct DataSimulation
    yinit::Vector
    D::Distribution
    maxrl::Int = Int(1e04)
    
    thetaHatVec = zeros(maxrl)
end