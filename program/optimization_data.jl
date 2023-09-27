# Contains the data for this project eg. bounds and costs and such

using SparseArrays

n_nodes = 11;

# Which generators are at each node
node_generators = [
    [],         # Node 1
    [1, 2, 3],  # Node 2
    [4],        # Node 3
    [5],        # Node 4
    [6],        # Node 5
    [],         # Node 6
    [7],        # Node 7
    [],         # Node 8
    [8, 9],     # Node 9
    [],         # Node 10
    []          # Node 11
]

# Which consumers are at each node
node_consumers = [
    [1],        # Node 1
    [],         # Node 2
    [],         # Node 3
    [2],        # Node 4
    [],         # Node 5
    [3],        # Node 6
    [],         # Node 7
    [4],        # Node 8
    [5],        # Node 9
    [6],        # Node 10
    [7]         # Node 11
]

# Data for the generators
generator_node = [2, 2, 2, 3, 4, 5, 7, 9, 9]
generator_capacity = [0.02, 0.15, 0.08, 0.07, 0.04, 0.17, 0.17, 0.26, 0.05]
generator_cost = [175, 100, 150, 150, 300, 350, 400, 300, 200]
generator_reactive_ub = 0.03 .* generator_capacity
n_generators = length(generator_node)
node_reactive_power = [length(indexes) > 0 ? sum(generator_reactive_ub[i] for i in indexes) : 0 for indexes in node_generators]

# Data for the consumers
consumer_node = [1, 4, 6, 8, 9, 10, 11]
consumer_demand = [0.10, 0.19, 0.11, 0.09, 0.21, 0.05, 0.04]

# Bounds for the values at the nodes
phi_lb = -pi
phi_ub = pi
voltage_lb = 0.98
voltage_ub = 1.02

# Data for the edges
b_coeff = [20.1, -22.3, -16.8, -17.2, -11.7, -19.4, -10.8, -12.3, -9.2, -13.9, -8.7, -11.3, -14.7, -13.5, -26.7]
g_coeff = [4.12, 5.67, 2.41, 2.78, 1.98, 3.23, 1.59, 1.71, 1.26, 1.11, 1.32, 2.01, 2.41, 2.14, 5.06]
pairs = [(1, 2), (1, 11), (2, 3), (2, 11), (3, 4), (3, 9), (4, 5), (5, 6), (5, 8), (6, 7), (7, 8), (7, 9), (8, 9), (9, 10), (10, 11)]

# Initializing the matrixes with coefficients for the edges
b_kl = sparse(zeros(n_nodes, n_nodes))
g_kl = sparse(zeros(n_nodes, n_nodes))
for (index, (k, l)) in enumerate(pairs)
    b_kl[k, l] = b_coeff[index]
    b_kl[l, k] = b_coeff[index]
    g_kl[k, l] = g_coeff[index]
    g_kl[l, k] = g_coeff[index]
end

function p_kl(v_k, v_l, phi_k, phi_l, k::Int, l::Int)::Float64
    return v_k^2.0 * g_kl[k, l] - v_k * v_l * g_kl[k, l]*cos(phi_k - phi_l) - v_k * v_l * b_kl[k, l] * sin(phi_k - phi_l)
end

function q_kl(v_k, v_l, phi_k, phi_l, k::Int, l::Int)::Float64
    return -v_k^2.0 * b_kl[k, l] + v_k * v_l * b_kl[k, l]*cos(phi_k - phi_l) - v_k * v_l * g_kl[k, l] * sin(phi_k - phi_l)
end
