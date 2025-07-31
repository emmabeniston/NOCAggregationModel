# Script to plot the concentration profiles generated from the 1D simulation

using Oceananigans
using CairoMakie
using Printf

# Correction factors - change if change the radii used in 1D_simulation.jl script
const r₁ = 5e-6 # Small particles <100μm in diameter, take diameter 10μm from optimised parameters in Williams et al.
const r₂ = 87.5e-6 # Medium particles 100-250μm in diameter, take midpoint
const r₃ = 562.5e-6 # Large particles 250-2000μm in diameter, take midpoint
const vol_correction₁ = 1e6*4/3*π*r₁^3 # Factor to convert from number conc to ppm
const vol_correction₂ = 1e6*4/3*π*r₂^3
const vol_correction₃ = 1e6*4/3*π*r₃^3

filename = "scripts/1D_simulation_tracers"
filepath = filename * ".jld2"

time_series = (C₁ = FieldTimeSeries(filepath, "C₁"),
               C₂ = FieldTimeSeries(filepath, "C₂"),
               C₃ = FieldTimeSeries(filepath, "C₃"),
               turbulent_dissipation = FieldTimeSeries(filepath, "turbulent_dissipation"))

# Coordinate arrays
xc, yc, zc = nodes(time_series.C₁)
times = time_series.C₁.times
Lx = 1
Ly = 1
Lz = 500

# Animate using observables
n = Observable(1)

 C₁ₙ = @lift interior(time_series.C₁[$n], 1, 1, :)*vol_correction₁ # Convert from number concentration to volume concentration
 C₂ₙ = @lift interior(time_series.C₂[$n], 1, 1, :)*vol_correction₂
 C₃ₙ = @lift interior(time_series.C₃[$n], 1, 1, :)*vol_correction₃
 turbulent_dissipationₙ = @lift interior(time_series.turbulent_dissipation[$n], 1, 1, :)

# Make figure and plots
fig = Figure(size = (750, 750))

axis_kwargs_C₁ = (xlabel="Volume concentration (ppm)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((0, 2), (-Lz, 0)))

axis_kwargs_C₂ = (xlabel="Volume concentration (ppm)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((0, 3), (-Lz, 0)))

axis_kwargs_C₃ = (xlabel="Volume concentration (ppm)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((0, 12), (-Lz, 0)))  

axis_kwargs_turbulent_dissipation = (xlabel="Turbulent dissipation rate (m² s⁻³)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((0, 1.1e-5), (-Lz, 0)))

ax_C₁ = Axis(fig[2, 1]; title = "Small particles, r = $(r₁) m", axis_kwargs_C₁...)
ax_C₂ = Axis(fig[2, 2]; title = "Medium particles, r = $(r₂) m", axis_kwargs_C₂...)
ax_C₃ = Axis(fig[3, 1]; title = "Large particles, r = $(r₃) m", axis_kwargs_C₃...)
ax_turbulent_dissipation = Axis(fig[3, 2]; title = "Turbulent dissipation", axis_kwargs_turbulent_dissipation...)

title = @lift @sprintf("t = %s", prettytime(times[$n]))
fig[1, 1:2] = Label(fig, title, fontsize=24, tellwidth=false)

lines!(ax_C₁, C₁ₙ, zc)

lines!(ax_C₂, C₂ₙ, zc)

lines!(ax_C₃, C₃ₙ, zc)

lines!(ax_turbulent_dissipation, turbulent_dissipationₙ, zc)

# Now record the movie

frames = 1:length(times)

@info "Making a motion picture..."

record(fig, "scripts/1D_volume_conc_profiles.mp4", frames, framerate=8) do i
    n[] = i
end
             