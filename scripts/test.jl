using Oceananigans
using Oceananigans.Units
using OceanBioME
using Oceanostics

include("../src/NOCAggregationModel.jl")

using .NOCAggregationModel

# Parameters
damping_rate = 1/86400 # For the relaxation forcing. Units 1/seconds. Change target functions in forcings section
r₁ = 10*1e-6 # radius in meters?
r₂ = 100*1e-6
r₃ = 200*1e-6
w_sink_speed₁ = 0.1 / 86400 # speed in meters per second?
w_sink_speed₂ = 0.2 / 86400
w_sink_speed₃ = 0.3 / 86400
α = 0.5
D = 2.25
ν₀ = 1e-6
Nx = 2
Ny = 2
Nz = 40
Lx = 2
Ly = 2
Lz = 40

grid = RectilinearGrid(CPU(), size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0))

# Set up biogeochemistry
biogeochemistry = NOCAggregationModelThreeCompartment(; grid, r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, α, D, ν₀)

# Set up relaxation forcing to base profiles
C₁_target(x, y, z, t) = -z
C₁_relaxation = Relaxation(rate=damping_rate, target = C₁_target)
C₂_target(x, y, z, t) = -10z
C₂_relaxation = Relaxation(rate=damping_rate, target = C₂_target)
C₃_target(x, y, z, t) = -100z
C₃_relaxation = Relaxation(rate=damping_rate, target = C₃_target)

model = NonhydrostaticModel(; grid,
                                     biogeochemistry,
                                     forcing = (C₁ = C₁_relaxation, C₂ = C₂_relaxation, C₃ = C₃_relaxation),
                                     closure = ScalarDiffusivity(ν = 0.0001, κ = 0.0001))

# Set initial conditions
# Use this function if only want particles initialised in top half of domain
Cᵢ(x, y, z) = z >= -Lz/2 ? 1e3 : 0 # Concentration in number per meter cubed
set!(model, C₁ = Cᵢ, C₂ = Cᵢ, C₃ = Cᵢ)

# Random noise perturbation to velocity to trigger variation in turbulent_shear for testing
uᵢ(x, y, z) = 1e-3 * randn()
vᵢ(x, y, z) = 1e-3 * randn()
wᵢ(x, y, z) = 1e-3 * randn()
set!(model, u=uᵢ, v=vᵢ, w=wᵢ)

# Set up simulation

simulation = Simulation(model, Δt = 0.1minutes, stop_time = 1day)

# Add callback to compute the turbulent dissipation rate:

turbulent_dissipation_op = KineticEnergyDissipationRate(model)
turbulent_dissipation = model.tracers.turbulent_dissipation

function compute_turbulent_dissipation!(sim)
    turbulent_dissipation .= turbulent_dissipation_op
    return nothing
end

simulation.callbacks[:turbulent_dissipation] = Callback(compute_turbulent_dissipation!)

# Save some output

simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers,
                                                       filename = "scripts/test_tracers.jld2",
                                                       schedule = TimeInterval(10minutes),
                                                       overwrite_existing = true)

# Run

run!(simulation)