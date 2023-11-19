module CupaJo

using FCSFiles, Statistics, PointDensityPlots, RecipesBase
import Tables

export Jo, gates, setgate, cleargate, getchannel

"""
```julia
Jo(flowsample)
Jo(table)
```
Create a `Jo` object representing flow cytometry data from a `table` (such as a `DataFrame`)
or a `FCSFiles.FlowSample` (which is created when loading a .fcs file using `FileIO`)

Indexing by `jo[channel]` or `jo[[c1,c2,c3]]` is implemented, and a plot recipe exists for this type to create summary plots.
"""
struct Jo{V}
    channels::Vector{Symbol}
    data :: Dict{Symbol,Vector{V}}
    gates :: Dict{Symbol,Vector{V}}
    function Jo(channels::Vector{Symbol},d::Dict{Symbol,Vector{V}},g::Dict{Symbol,Vector{V}}) where {V<:Number}
        #all channels should be in keys(gates)
        @assert all(c in keys(g) for c in channels)
        new{V}(channels,d,g)
    end

    #Create empty gate Dict if not provided
    function Jo(c::Vector{Symbol},d::Dict{Symbol,Vector{V}}) where {V<:Number}
        g=Dict{Symbol,Vector{V}}(ch => Vector{V}() for ch in c)
        new{V}(c,d,g)
    end
end

function Jo(fs::FlowSample)
    channames = eachindex(fs) |> collect
    data=Dict(Symbol(cn) => collect(fs[cn]) for cn in channames)
    Jo(Symbol.(channames),data)
end

function Jo(table)
    @assert Tables.istable(table)
    columns=Tables.columns(table)
    channels = Tables.columnnames(columns)
    data = Dict(cn => Tables.getcolumn(columns,cn) for cn in channels)
    Jo(channels,data)    
end

Base.names(j::Jo) = j.channels

#Implement the tables.jl interface so we can materialize to DataFrames
Tables.istable(::Type{<:Jo}) = true
Tables.columnaccess(::Type{<:Jo}) = true
Tables.columns(j::Jo) = j
Tables.getcolumn(j::Jo,i::Int) = j.data[j.channels[i]]
Tables.getcolumn(j::Jo,nm::Symbol) = j.data[nm]
Tables.columnnames(j::Jo) = j.channels

#tools for working with gates
"""
```julia
gates(jo)
```
Return all of the gates assigned to a `Jo` struct as a `Dict`
"""
gates(j::Jo) = j.gates

"""
```julia
setgate(jo,channel,value)
```
Add a gate to `jo` at a given `value` for the provided `channel`
"""
function setgate end

function setgate(j::Jo,channel::Symbol,values::Vector)
    @assert channel in names(j) "can only add gates to a channel which exists"
    push!(gates(j)[channel],values...)
end

setgate(j::Jo,channel::Symbol,value::Number) = push!(gates(j)[channel],value)


setgate(j::Jo,channel::String,value) = setgate(j,Symbol(channel),value)

"""
```julia
cleargate(jo [,channel])
```
Clear gates corresponding to a given `channel`. If `channel` is omitted, clear
all gates
"""
function cleargate end

function cleargate(j::Jo{V},channel::Symbol) where {V}
    gates(j)[channel] = Vector{V}()
end

cleargate(j::Jo,c::String) = cleargate(j,Symbol(c))


function cleargate(j::Jo)
    for c in names(j)
        cleargate(j,c)
    end        
end

#methods for pulling channel data
"""
```julia
getchannel(jo,channel)
```
Get the raw values for a single `channel`. Equivalent to `jo[channel]`
"""
function getchannel end

getchannel(j::Jo,channel::Symbol) = j.data[channel]

getchannel(j::Jo,channel::String) = getchannel(j,Symbol(channel))

#make indexing work
Base.getindex(j::Jo,channel::Union{String,Symbol}) = getchannel(j,channel)

function Base.getindex(j::Jo{V},channels::Vector{<:Union{String,Symbol}}) where {V}
    channels = Symbol.(channels)
    d=Dict{Symbol,Vector{V}}()
    g=Dict{Symbol,Vector{V}}()
    for c in channels
        #copy over data
        d[c] = getchannel(j,c)
        g[c] = gates(j)[c]
    end
    Jo(channels,d,g)
end
    
#first make a user recipe for flow point density plots
"""
```julia
flowdensity(x,y;originquantile=0.01)
```
Create a `flowdensity` plot. If negative values exist in the data, these values will
be replaced with values corresponding to the `originquantile` for plotting.
"""
function flowdensity end

@userplot FlowDensity
@recipe function f(fd::FlowDensity)
    xscale --> :log10
    yscale --> :log10
    originquantile --> .01
    seriestype := :pointdensity
    (ox,oy) = fd.args

    #if all values are greater than zero, no editing needed
    if all(ox .> 0) && all(oy .> 0)
        return (ox,oy)
    end
    
    #move very small values to the origin quantile for that axis
    (newx,newy) = map((ox,oy)) do a
        apos=filter(a) do ai
            ai > 0
        end
        aorigin=quantile(apos,plotattributes[:originquantile])
        map(a) do ai
            (ai < aorigin) ? aorigin : ai
        end
    end
    (newx,newy)
end

#build a plot recipe for Jo structs

@recipe function plotjo(j::Jo)
    #build densityplots for each combination of channels, then a histogram for each channel
    nchans = names(j) |> length
    comboindices=collect((i,k) for i in 1:nchans for k in (i+1):nchans)
    chancombos=map(comboindices) do (i,k)
        names(j)[i],names(j)[k]
    end
    layout := @layout [RecipesBase.grid(length(chancombos),1) RecipesBase.grid(nchans,1)]
    plotnum=1
    legend --> false
    for cc in chancombos
        @series begin
            xguide := cc[1]
            yguide := cc[2]
            subplot := plotnum
            FlowDensity((getchannel(j,cc[1]),getchannel(j,cc[2])))
        end

        @series begin
            subplot := plotnum
            xguide := cc[1]
            yguide := cc[2]
            seriestype := :vline
            linecolor --> :black
            gates(j)[cc[1]]
        end

        @series begin
            subplot := plotnum
            xguide := cc[1]
            yguide := cc[2]
            seriestype := :hline
            linecolor --> :black
            gates(j)[cc[2]]
        end
        
        plotnum += 1
    end

    for c in names(j)
        @series begin
            subplot := plotnum
            seriestype := :histogram
            xguide := c
            legend --> false
            getchannel(j,c)
        end

        @series begin
            subplot := plotnum
            seriestype := :vline
            xguide := c
            linecolor --> :black
            gates(j)[c]
        end
        
        plotnum += 1
    end
    ()
end

end # module CupaJo
