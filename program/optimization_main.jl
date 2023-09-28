using JuMP
import Ipopt

include("data.jl")

# Create the model object
model = Model(Ipopt.Optimizer)

# We set three variables for solving the problem
# p: power generation by generator i
# phi: phase at node j
# v: voltage at node j
@variable(model, 0 <= p[i = 1:n_generators] <= generator_capacity[i])
@variable(model, phi_lb <= phi[1:n_nodes] <= phi_ub)
@variable(model, voltage_lb <= v[1:n_nodes] <= voltage_ub)

# Set the objective function which would minimise the sum of the cost for all generators
@objective(model, Min, sum(generator_cost[i] * p[i] for i in 1:n_generators))

# Loop through the nodes, then calculate and set the individual constraints for all of them
# sum(p_kl(v[i], v[j], phi[i], phi[j], i, j) for j in 1:nodes) is the power flow into a node i
active_power_constraints = []
reactive_power_constraints = []
for i in 1:n_nodes
    # Contrain the flow of (active power + generation power - consumed power) at each node to be 0
    constraint = @NLconstraint(model, 
        - sum(
                (v[i]^2.0 * g_kl[i, j] - v[i] * v[j] * g_kl[i, j]*cos(phi[i] - phi[j]) - v[i] * v[j] * b_kl[i, j] * sin(phi[i] - phi[j]))
                for j in 1:n_nodes
            )  
        +sum(p[j] for j in node_generators[i])  
        -sum(consumer_demand[j] for j in node_consumers[i])
        == 0
    )
    push!(active_power_constraints, constraint)
    # Constrain the net flow of reactive power flowing into each node to be within the allowed limits for the generators at this location 
    constraint = @NLconstraint(model, 
        - node_reactive_power[i]
        <= sum(
                -v[i]^2.0 * b_kl[i, j] + v[i] * v[j] * b_kl[i, j]*cos(phi[i] - phi[j]) - v[i] * v[j] * g_kl[i, j] * sin(phi[i] - phi[j])
                for j in 1:n_nodes
            ) 
        <= node_reactive_power[i]
    )
    push!(reactive_power_constraints, constraint)
end

# Print the optimzation problem in the terminal
# println(model)

# Solve the optimization problem
optimize!(model)

# Printing some of the results for further analysis
println("") # Printing white line after solver output, before printing
println("Termination statue: ", JuMP.termination_status(model))
println("Optimal(?) objective function value: ", JuMP.objective_value(model))
println("\nOptimal(?) power: ", round.(JuMP.value.(p), digits=6))
println("Optimal(?) voltage: ", round.(JuMP.value.(v), digits=6))
println("Optimal(?) phase: ", round.(JuMP.value.(phi), digits=6))
println("\nDual variables/Lagrange multipliers corresponding to some of the constraints: ")
println("active power constraints", round.(JuMP.dual.(active_power_constraints), digits=6))
println("reactive power constraints", round.(JuMP.dual.(reactive_power_constraints), digits=6))
println("power: ", round.(JuMP.dual.(JuMP.UpperBoundRef.(p)), digits=6))
println("voltage: ", round.(JuMP.dual.(JuMP.UpperBoundRef.(v)), digits=6))
println("phase: ", round.(JuMP.dual.(JuMP.UpperBoundRef.(phi)), digits=6))

p = JuMP.value.(p)
v = JuMP.value.(v)
phi = JuMP.value.(phi)

flow_active = zeros(n_nodes, n_nodes)
flow_reactive = zeros(n_nodes, n_nodes)

println()
for (k, l) in pairs
    flow_active[k, l] = p_kl(v[k], v[l], phi[k], phi[l], k, l)
    flow_reactive[k, l] = q_kl(v[k], v[l], phi[k], phi[l], k, l)
    flow_active[l, k] = p_kl(v[l], v[k], phi[l], phi[k], l, k)
    flow_reactive[l, k] = q_kl(v[l], v[k], phi[l], phi[k], l, k)
    println("Active flow from node $k to $l: ", flow_active[k, l])
    println("Active flow from node $l to $k: ", flow_active[l, k])
end

# Calculate the net flow into or from
net_flow = sum(flow_active, dims=2)
println("\n Net flow from (negative values indicate net flow to): ", round.(net_flow, digits=6))

# All flows as matrixes
println("\nFlows as matrixes")
println(round.(flow_active, digits=6))
println(round.(flow_reactive, digits=6))
