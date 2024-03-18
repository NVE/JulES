# Aggregate modelobjects and remove modelobjects not relevant for subsystems
function removeelements!(elements::Vector{DataElement}; aggzone::Dict=Dict(), rm_basebalances::Bool=true) # TODO: Replace with more user settings
    aggzonecopl = Dict()
    for (k,v) in aggzone
        for vv in v
            aggzonecopl["PowerBalance_" * vv] = "PowerBalance_" * k
        end
    end
    
    delix = []
    powerbasebalances = []
    for (i,element) in enumerate(elements)
        if rm_basebalances
            # BaseBalances
            if element.typename == "BaseBalance"
                if element.value["Commodity"] == "Power"
                    push!(delix,i)
                    push!(powerbasebalances, element.instancename)
                end
            end

            # Power balance RHSTerms
            if element.conceptname == "RHSTerm"
                if element.value["Balance"] in powerbasebalances
                    push!(delix,i)
                end
            end
            
            # Residualhints
            if element.typename == "Residualhint"
                push!(delix,i)
            end
        end
        
        # Power balance Arrows
        if element.conceptname == "Arrow"
            if element.value["Balance"] in keys(aggzonecopl)
                value = copy(element.value)
                value["Balance"] = aggzonecopl[element.value["Balance"]]
                elements[i] = DataElement(element.conceptname, element.typename, element.instancename, value)
            end
        end

        if element.typename == "HydroRampingWithout"
            push!(delix,i)
        end
    end
    for deli in sort(delix; rev=true)
        popat!(elements, deli)
    end
    return elements
end