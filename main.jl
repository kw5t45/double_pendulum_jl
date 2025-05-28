import Pkg; Pkg.add("DifferentialEquations")
import Pkg; Pkg.add("Plots")
import Pkg; Pkg.add("StatsBase")

using DifferentialEquations
using Plots

function double_pendulum!(du, u, p, t)
    # u = [θ₁, ω₁, θ₂, ω₂]
    θ1, ω1, θ2, ω2 = u
    m1, m2, L1, L2, g = p

    Δ = θ2 - θ1

    den1 = (m1 + m2) * L1 - m2 * L1 * cos(Δ)^2
    den2 = (L2 / L1) * den1

    du[1] = ω1

    du[2] = (m2 * L1 * ω1^2 * sin(Δ) * cos(Δ) +
             m2 * g * sin(θ2) * cos(Δ) +
             m2 * L2 * ω2^2 * sin(Δ) -
             (m1 + m2) * g * sin(θ1)) / den1

    du[3] = ω2

    du[4] = (-m2 * L2 * ω2^2 * sin(Δ) * cos(Δ) +
             (m1 + m2) * g * sin(θ1) * cos(Δ) -
             (m1 + m2) * L1 * ω1^2 * sin(Δ) -
             (m1 + m2) * g * sin(θ2)) / den2
end

# Αρχικές συνθήκες: [θ₁, ω₁, θ₂, ω₂]
u₀ = [π / 2, 0.0, π / 2, 0.0]  # Ξεκινάνε από την ίδια θέση
u₀_2 = [π / 2, 0.0, π / 2, 0.000001]  # Ξεκινάνε από την ίδια θέση

# Παράμετροι: m₁, m₂, L₁, L₂, g
p = (10.0, 1.0, 1.0, 1.0, 9.81)

# Χρονικό διάστημα
tspan = (0.0, 10.0)


prob = ODEProblem(double_pendulum!, u₀, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9)

prob_2 = ODEProblem(double_pendulum!, u₀_2, tspan, p)
sol_2 = solve(prob_2, Tsit5(), reltol=1e-9, abstol=1e-9)

# plot(sol.t, sol[1,:], label="θ₁(t)", xlabel="t", ylabel="θ₁ - θ1-2", title="Γωνιακή Θέση θ₁ και θ1-2")
# plot!(sol_2.t, sol_2[1,:], label="θ₁-2(t)", xlabel="t")
# plot(sol.t, sol[1,:], label="θ₁(t)", xlabel="t", ylabel="θ₁", title="Γωνιακή Θέση θ₁ και θ₂")
# plot!(sol.t, sol[3,:], label="θ₂(t)", ylabel="θ (rad)")
# plot(sol.t, sol[1, :], label="θ₁(t)", title="Γωνιακή θέση θ₁(rad)", xlabel="t", ylabel="θ₁")
# # Εξαγωγή λύσης
# θ1 = sol[1, :]
# θ2 = sol[3, :]

# # Μήκη ράβδων
# L1 = p[3]
# L2 = p[4]

# # Υπολογισμός θέσης στο επίπεδο (x₂, y₂)
# x1 = L1 * sin.(θ1)
# y1 = -L1 * cos.(θ1)

# x2 = x1 .+ L2 * sin.(θ2)
# y2 = y1 .- L2 * cos.(θ2)

# # Σχεδίαση τροχιάς του δεύτερου εκκρεμούς
# plot(x2, y2, xlabel="x", ylabel="y", title="Τροχιά του δεύτερου εκκρεμούς στο επίπεδο", legend=false, lw=1.5, aspect_ratio=1)
# Poincaré section: κρατάμε σημεία όπου θ₁ ≈ π/2 και ω₁ > 0
θ2_vals = Float64[]
ω2_vals = Float64[]

for i in 2:length(sol.t)-1
    θ1_prev, ω1_prev = sol[i-1][1], sol[i-1][2]
    θ1,     ω1     = sol[i][1], sol[i][2]
    θ1_next        = sol[i+1][1]

    ϵ = 0.2
    if abs(θ1) < ϵ && ω1 < 0
θ2 = mod(sol[i][3], 2π)
        ω2 = sol[i][4]
        push!(θ2_vals, θ2)
        push!(ω2_vals, ω2)
    end
end

# Πλοκή διαγράμματος Poincaré
scatter(θ2_vals, ω2_vals, markersize=2, legend=false)
xlabel!("θ₂ (rad) modulo 2π")
ylabel!("ω₂ (rad/s)")
title!("Διάγραμμα Poincaré: Τομή για θ₁ ≈ 0 και ω₁ > 0")

using StatsBase  # edw to plots den doylevei

# Δεδομένα: θ₂ και ω₂
x = θ2_vals
y = ω2_vals

# Ορισμός bins (π.χ. 100x100)
xbins = range(0, 2π, length=100)
ybins = range(-15, 15, length=100) 

# Φτιάχνουμε 2D histogram
h = fit(Histogram, (x, y), (xbins, ybins))

# Προβολή με heatmap
heatmap(h.edges[1], h.edges[2], h.weights',
        xlabel="θ₂ (rad) mod π", ylabel="ω₂ (rad/s)",
        xlims=(0, 2π),
        title="Poincaré Heatmap: Τομή για θ₁ ≈ 0 και ω₁ < 0",
        color=:viridis, colorbar=false,
        aspect_ratio=0.1)
# Παράμετροι του συστήματος
p = (1.0, 1.0, 1.0, 1.0, 9.81)
tspan = (0.0, 10.0)

# Πλέγμα αρχικών συνθηκών
θ1_range = range(-π, π, length=300)
θ2_range = range(-π, π, length=300)
output = zeros(length(θ1_range), length(θ2_range))

# Υπολογισμός του χρόνου μέχρι την πρώτη ανατροπή
for (i, θ1₀) in enumerate(θ1_range)
    for (j, θ2₀) in enumerate(θ2_range)
        u₀ = [θ1₀, 0.0, θ2₀, 0.0]
        prob = ODEProblem(double_pendulum!, u₀, tspan, p)
        sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9)

        # Έλεγχος για ανατροπή
        flipped = false
        for u in sol.u
            if abs(u[3]) > π
                output[i, j] = sol.t[findfirst(==(u), sol.u)]
                flipped = true
                break
            end
        end
        if !flipped
            output[i, j] = NaN  # Δεν έγινε ανατροπή
        end
    end
end

# Απεικόνιση του fractal
heatmap(θ2_range, θ1_range, output;
        xlabel="θ₂(0)", ylabel="θ₁(0)",
        title="Χρόνος μέχρι την πρώτη ανατροπή",
        color=:viridis,
        clims=(0, maximum(filter(!isnan, output))),
        colorbar=true)
# Παράμετροι του συστήματος
p = (1.0, 1.0, 1.0, 1.0, 9.81)
tspan = (0.0, 5.0)  # Σταθερός χρόνος προσομοίωσης

# Πλέγμα αρχικών συνθηκών
θ1_range = range(-π, π, length=1000)
θ2_range = range(-π, π, length=1000)
output = zeros(length(θ1_range), length(θ2_range))

# Υπολογισμός της γωνίας θ₂ μετά από 5 δευτερόλεπτα
for (i, θ1₀) in enumerate(θ1_range)
    for (j, θ2₀) in enumerate(θ2_range)
        u₀ = [θ1₀, 0.0, θ2₀, 0.0]
        prob = ODEProblem(double_pendulum!, u₀, tspan, p)
        sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9)

        θ2_final = sol.u[end][3]
        output[i, j] = mod(θ2_final + π, 2π) - π  # Περιορισμός στο [-π, π]
    end
end

# Απεικόνιση του φράκταλ
heatmap(θ2_range, θ1_range, output;
        xlabel="θ₂(0)", ylabel="θ₁(0)",
        title="Γωνία θ₂ μετά από 5 δευτερόλεπτα",
        color=:hsv, clims=(-π, π),
        colorbar=true)