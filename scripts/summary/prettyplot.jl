#? Some pretty plots for the control charts which display both the
#? control limit and the cautious limit

using Plots
plot(tmp[:chart_values], legend=:outerright, label="C")
hline!([-tmp[:limit_alarm],tmp[:limit_alarm]], fill=true, fillcolour=:orange, fillalpha=0.15, color=:darkred, style=:dash, label=false)
hline!([-tmp[:limit_cautious],tmp[:limit_cautious]], fill=true, fillcolour=:lightgreen, fillalpha=0.35, color=:darkgreen, style=:dash, label=false)
if !IC
    scatter!([tmp[:t_alarm]], [tmp[:chart_values][end]], label=false, color=:red)    
    vline!([tau], style=:dash, colour=:darkblue, label="tau")
end

plot(tmp[:parameter_updates])
if !IC
    vline!([tau], style=:dash, colour=:darkblue, label="tau")
end
plot(tmp[:di])
if !IC
    vline!([tau], style=:dash, colour=:darkblue, label="tau")
end


function fDist(n, m)
    rand(DiscreteUniform(0,m), n)  
end

function fVec(n,m)
    rand(0:m, n)
end