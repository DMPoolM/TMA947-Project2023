using JuMP
import Ipopt

#include("project_data.jl")

model = Model(Ipopt.Optimizer)

@variable(model, 0 <= act_power[i = 1:generator_num] <= generator_capacity[i])
@variable(model, theta_lb <= vol_angle[1:num_nodes] <= theta_ub)
@variable(model, voltage_lb <= theta[1:num_nodes] <= voltage_ub)
@variable(model, reactivepower_lb[i] <= react_power[i = 1:generator_num] <= reactivepower_ub[i])

@objective(model, Min, sum(act_power[i] * generator_cost[i] for i in 1:generator_num))

active_power_constraints = []
reactive_power_constraints = []
for i in 1 : num_nodes
  constraint = @NLconstraint(model, 
    - sum(
            (theta[i]^2.0 * g_kl[i, j] - theta[i] * theta[j] * g_kl[i, j]*cos(vol_angle[i] - vol_angle[j]) - theta[i] * theta[j] * b_kl[i, j] * sin(vol_angle[i] - vol_angle[j])) 
            for j in 1:num_nodes
            )
    + sum(
            act_power[j] for j in nodes_with_generators[i]
            )  
    - sum(
            consumerDemand[j] for j in nodes_with_consumers[i]
            )
    == 0)
  push!(active_power_constraints, constraint)
  constraint_1 = @NLconstraint(model, 
      sum(
            react_power[k] for k in nodes_with_generators[i]
            )
    - sum(
            -theta[i]^2.0 * b_kl[i, j] + theta[i] * theta[j] * b_kl[i, j]*cos(vol_angle[i] - vol_angle[j]) - theta[i] * theta[j] * g_kl[i, j] * sin(vol_angle[i] - vol_angle[j]) 
            for j in 1:num_nodes
            ) 
    == 0)
    push!(reactive_power_constraints, constraint_1)
end

optimize!(model)

println()
println("Termination statue: ", JuMP.termination_status(model))
println("Optimal(?) objective function value: ", JuMP.objective_value(model))


println("Power generated: ", round.(JuMP.value.(act_power), digits=6))
println("Max power: ", round.(JuMP.value.(generator_capacity), digits=6))
println("Reactive power generated: ", round.(JuMP.value.(react_power), digits=6))
println("Reactive max +/-: ", round.(JuMP.value.(reactivepower_ub), digits=6))
println("Voltage amplitude: ", round.(JuMP.value.(theta), digits=6))
println("Voltage angle: ", round.(JuMP.value.(vol_angle), digits=6))
println("\nDual variables/Lagrange multipliers corresponding to some of the constraints: ")
println("active power constraints", round.(JuMP.dual.(active_power_constraints), digits=6))
println("reactive power constraints", round.(JuMP.dual.(reactive_power_constraints), digits=6))
println(JuMP.dual.(JuMP.UpperBoundRef.(act_power)))

act_power = JuMP.value.(act_power)
theta = JuMP.value.(theta)
vol_angle = JuMP.value.(vol_angle)
flow_active = zeros(num_nodes, num_nodes)
flow_reactive = zeros(num_nodes, num_nodes)'

println()
for (k, l) in edges
    flow_active[k, l] = p_kl(theta[k], theta[l], vol_angle[k], vol_angle[l], k, l)
    flow_reactive[k, l] = q_kl(theta[k], theta[l], vol_angle[k], vol_angle[l], k, l)
    flow_active[l, k] = p_kl(theta[l], theta[k], vol_angle[l], vol_angle[k], l, k)
    flow_reactive[l, k] = q_kl(theta[l], theta[k], vol_angle[l], vol_angle[k], l, k)
end

println()
for k in 1 : num_nodes
  for l in 1 : num_nodes
    if (flow_active[k,l] > 0)
      println("Flow of active power from $k to $l: ", round.(flow_active[k,l], digits = 6))
    end
  end
end

println()
for k in 1 : num_nodes
  for l in 1 : num_nodes
    if (flow_reactive[k,l] > 0)
      println("Flow of reactive power from $k to $l: ", round.(flow_reactive[k,l], digits = 6))
    end
  end
end
