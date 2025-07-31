# This script is setup to run a 1D simulation with the NOCAggregationModelThreeCompartment
# We take parameters such as radius and sinking speeds from the Williams et al. paper
# We will use an observed turbulence profile to force the aggregation
# We will relax the smallest size class profiles to the observed profile
# Aim is to reproduce observations, or diagnose why this can't be done
# To begin, the profiles we have used are approximations to the observations in the TS region presented in Williams et al. 

using Oceananigans
using Oceananigans.Units
using OceanBioME
using Oceanostics
using Printf

include("../src/NOCAggregationModel.jl")

using .NOCAggregationModel

# Aggregation parameters
const r₁ = 5e-6 # Small particles <100μm in diameter, take diameter 10μm from optimised parameters in Williams et al.
const r₂ = 87.5e-6 # Medium particles 100-250μm in diameter, take midpoint
const r₃ = 562.5e-6 # Large particles 250-2000μm in diameter, take midpoint
const w_sink_speed₁ = 100 * (2*r₁*1e3)^0.56 / 86400 # speed in meters per second, form taken from Williams et al.
const w_sink_speed₂ = 100 * (2*r₂*1e3)^0.56 / 86400
const w_sink_speed₃ = 100 * (2*r₃*1e3)^0.56 / 86400
const stick = 0.5 # Taken from Williams et al. paper - optimised α and D
const D = 2.25
const ν₀ = 9.86e-7

const C₁_correction = 2 * π * r₁^2 # correction used to calculate the number concentration of small particles from the beam attenuation coefficient (line 191 of Williams et al. manuscript)
const C₂_correction = 4/3 * π * r₂^3 # correction to calculate number conc of medium particles from the volume conc data - volume of one particle
const C₃_correction = 4/3 * π * r₃^3 # as above, but large particles

# Grid and simulation parameters
const Nx = Ny = 1
const Nz = 128
const Lx = Ly = 1
const Lz = 500
const Δt = 10seconds
const stop_time = 7days
const κ = 0.0001 # scalar diffusivity for the concentration fields
const damping_rate = 1/(day/24) # Rate for relaxation forcing for concentration profiles. Units 1/seconds. Will need to test and try varying.

# Observations for turbulence and concentration profiles
function turbulence_observed_profile(x, y, z, t) # Approximation to the observed profile for turbulent dissipation rate
    if z >= -25
        1e-5
    elseif z >= -50
        10^(4*z/25-1)
    elseif z>= -300
        1e-9
    else
        0
    end
end

function C₁_observed_profile(x, y, z, t) # Approximation for observed profiles for number concentration of particles - for small is from beam attenuation data
    if z >= -40
        0.4/C₁_correction
    elseif z >= -100
        (7*z/1200 + 19/30)/C₁_correction
    elseif z >= -300
        (z/4000 + 3/40)/C₁_correction
    else
        0
    end
end

C₁_relaxation_on = true # Toggle whether or not the concentration is relaxed to the observed profile or not

function C₂_observed_profile(x, y, z, t) # Approximation for observed profiles for number concentration of particles - for medium is from volume conc data in ppm (hence 1e-6/C₂_correction)
    if z >= -10
        0.3*1e-6/C₂_correction
    elseif z >= -50
        (-z/40 + 1/20)*1e-6/C₂_correction
    elseif z >= -100
        (z/50 + 23/10)*1e-6/C₂_correction
    elseif z >= -500
        (z/1600 + 29/80)*1e-6/C₂_correction
    else
        0
    end
end

C₂_relaxation_on = false

function C₃_observed_profile(x, y, z, t) # Approximation for observed profiles for number concentration of particles - for large is from volume conc data in ppm (hence 1e-6/C₃_correction)
    if z >= -10 
        1.5*1e-6/C₃_correction
    elseif z >= -70
        (-z/8 + 1/4)*1e-6/C₃_correction
    elseif z >= -500
        (2*z/215 + 415/43)*1e-6/C₃_correction
    else
        0
    end
end

C₃_relaxation_on = false

# Grid
grid = RectilinearGrid(GPU(), size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = (-Lz, 0))

# Biogeochemistry - aggregation
biogeochemistry = NOCAggregationModelThreeCompartment(; grid, r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀)

# Relaxation forcing to observed profiles
C₁_relaxation = C₁_relaxation_on ? Relaxation(rate=damping_rate, target = C₁_observed_profile) : nothing
C₂_relaxation = C₂_relaxation_on ? Relaxation(rate=damping_rate, target = C₂_observed_profile) : nothing
C₃_relaxation = C₃_relaxation_on ? Relaxation(rate=damping_rate, target = C₃_observed_profile) : nothing

# Model
model = NonhydrostaticModel(; grid,
                                biogeochemistry,
                                forcing = (C₁ = C₁_relaxation, C₂ = C₂_relaxation, C₃ = C₃_relaxation),
                                closure = ScalarDiffusivity(ν = 0.0001, κ = (C₁ = κ, C₂ = κ, C₃ = κ, turbulent_dissipation = 0)))

# Initial conditions
# Turbulent dissipation and concentration are set to the observed profiles
turbulence_initial(x, y, z) = turbulence_observed_profile(x, y, z, 0)
C₁_initial(x, y, z) = C₁_observed_profile(x, y, z, 0)
C₂_initial(x, y, z) = C₂_observed_profile(x, y, z, 0)
C₃_initial(x, y, z) = C₃_observed_profile(x, y, z, 0)
set!(model, C₁ = C₁_initial, C₂ = C₂_initial, C₃ = C₃_initial, turbulent_dissipation = turbulence_initial)

# Set up simulation
simulation = Simulation(model, Δt = Δt, stop_time = stop_time)

# Progress printing
io = open("progress_printing.txt", "w")
write(io, "Progress printing\n")
close(io)

function progress_message(sim) 
	io = open("progress_printing.txt", "a");
	@printf(io, "Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time));
	close(io);
end

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(10))

# Output

# Full tracer fields (would be horizontal averages in 3D simulation)
simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers,
                                                       filename = "scripts/1D_simulation_tracers.jld2",
                                                       schedule = TimeInterval(10minutes),
                                                       overwrite_existing = true)

# Growth and loss rates due to aggregation and consumption
z_centers = CenterField(grid)
set!(z_centers, (x, y, z) -> z) # To be used in the consumption calculation - gives the z center coordinates
# Don't output growth for C₁ and loss for C₃ as these are just 0
loss_ts₁ = Field(loss_ts₁_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
loss_ds₁ = Field(loss_ds₁_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
consumption₁ = Field(consumption_func(z_centers, model.tracers.C₁))
growth_ts₂ = Field(growth_ts₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
growth_ds₂ = Field(growth_ds₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
loss_ts₂ = Field(loss_ts₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
loss_ds₂ = Field(loss_ds₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
consumption₂ = Field(consumption_func(z_centers, model.tracers.C₂))
growth_ts₃ = Field(growth_ts₃_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
growth_ds₃ = Field(growth_ds₃_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, model.tracers.C₁, model.tracers.C₂, model.tracers.C₃, model.tracers.turbulent_dissipation))
consumption₃ = Field(consumption_func(z_centers, model.tracers.C₃))

simulation.output_writers[:growth_and_loss] = JLD2Writer(model, (; loss_ts₁, loss_ds₁, consumption₁, growth_ts₂, growth_ds₂, loss_ts₂, loss_ds₂, consumption₂, growth_ts₃, growth_ds₃, consumption₃),
                                                       filename = "scripts/1D_simulation_growth_and_loss.jld2",
                                                       schedule = TimeInterval(10minutes),
                                                       overwrite_existing = true)                                                       

# Run

run!(simulation)