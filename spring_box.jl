using Plots
theme(:rose_pine_dawn)

# ----------------------------------------------------
# Physical constants
# ----------------------------------------------------
const g = 9.81
const L = 1.0

const THETA_0 = π / 3
const THETA_DOT_0 = 0.0
Δt = 0.01
t_max = 10.0

# ----------------------------------------------------
# ODE
# ----------------------------------------------------
function theta_double_dot(theta, theta_dot, μ)
    -μ * theta_dot - (g / L) * sin(theta)
end

# ----------------------------------------------------
# Solver
# ----------------------------------------------------
function solve_pendulum(μ)
    θ = THETA_0
    θ̇ = THETA_DOT_0

    t_vals = 0:Δt:t_max
    θ_vals = Float64[]
    θ̇_vals = Float64[]

    for _ in t_vals
        push!(θ_vals, θ)
        push!(θ̇_vals, θ̇)

        θ̈ = theta_double_dot(θ, θ̇, μ)

        θ += θ̇ * Δt
        θ̇ += θ̈ * Δt
    end

    return t_vals, θ_vals, θ̇_vals
end

# ----------------------------------------------------
# Solve
# ----------------------------------------------------
t, θ1, θ̇1 = solve_pendulum(0.0)
t, θ2, θ̇2 = solve_pendulum(0.1)
t, θ3, θ̇3 = solve_pendulum(0.4)

# ----------------------------------------------------
# Plots
# ----------------------------------------------------
p1 = plot(t, θ1, lw=1.5, label="μ = 0.0",
          xlabel="Time (s)", ylabel="θ (rad)",
          title="Angle vs Time")
plot!(p1, t, θ2, lw=1.5, label="μ = 0.1")
plot!(p1, t, θ3, lw=1.5, label="μ = 0.4")

p2 = plot(t, θ̇1, lw=1.5, label="μ = 0.0",
          xlabel="Time (s)", ylabel="θ̇ (rad/s)",
          title="Angular Velocity vs Time")
plot!(p2, t, θ̇2, lw=1.5, label="μ = 0.1")
plot!(p2, t, θ̇3, lw=1.5, label="μ = 0.4")

p3 = plot(θ1, θ̇1, lw=3, label="μ = 0.0",
          xlabel="θ (rad)", ylabel="θ̇ (rad/s)",
          title="Phase Space")
plot!(p3, θ2, θ̇2, lw=3, label="μ = 0.1")
plot!(p3, θ3, θ̇3, lw=3, label="μ = 0.4")

# ----------------------------------------------------
# Layout
# ----------------------------------------------------
fig = plot(
    p1, p2, p3;
    layout = [1 1;
              1 3],
    size = (700, 700)
)


display(fig)
savefig("C:/Users/jarri/Desktop/pendulum_mu.svg")
