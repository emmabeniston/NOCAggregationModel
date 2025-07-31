# Script to plot the concentration profiles generated from the 1D simulation

using Oceananigans
using CairoMakie
using Printf
using Oceananigans.Units

# Correction factors - change if change the radii used in 1D_simulation.jl script
const r₁ = 5e-6 # Small particles <100μm in diameter, take diameter 10μm from optimised parameters in Williams et al.
const r₂ = 87.5e-6 # Medium particles 100-250μm in diameter, take midpoint
const r₃ = 562.5e-6 # Large particles 250-2000μm in diameter, take midpoint
const vol_day_correction₁ = 1e6*4/3*π*r₁^3*day # Factor to convert from number conc per second to ppm per day
const vol_day_correction₂ = 1e6*4/3*π*r₂^3*day
const vol_day_correction₃ = 1e6*4/3*π*r₃^3*day

filename = "scripts/1D_simulation_growth_and_loss"
filepath = filename * ".jld2"

time_series = (loss_ts₁ = FieldTimeSeries(filepath, "loss_ts₁"),
               loss_ds₁ = FieldTimeSeries(filepath, "loss_ds₁"),
               consumption₁ = FieldTimeSeries(filepath, "consumption₁"),
               growth_ts₂ = FieldTimeSeries(filepath, "growth_ts₂"),
               growth_ds₂ = FieldTimeSeries(filepath, "growth_ds₂"),
               loss_ts₂ = FieldTimeSeries(filepath, "loss_ts₂"),
               loss_ds₂ = FieldTimeSeries(filepath, "loss_ds₂"),
               consumption₂ = FieldTimeSeries(filepath, "consumption₂"),
               growth_ts₃ = FieldTimeSeries(filepath, "growth_ts₃"),
               growth_ds₃ = FieldTimeSeries(filepath, "growth_ds₃"),
               consumption₃ = FieldTimeSeries(filepath, "consumption₃"))

# Coordinate arrays
xc, yc, zc = nodes(time_series.loss_ts₁)
times = time_series.loss_ts₁.times
Lx = 1
Ly = 1
Lz = 500

# Animate using observables
n = Observable(1)

 loss_ts₁ₙ = @lift -interior(time_series.loss_ts₁[$n], 1, 1, :)*vol_day_correction₁ # Convert from number conc per second to ppm per day
 loss_ds₁ₙ = @lift -interior(time_series.loss_ds₁[$n], 1, 1, :)*vol_day_correction₁
 consumption₁ₙ = @lift -interior(time_series.consumption₁[$n], 1, 1, :)*vol_day_correction₁
 total₁ₙ = @lift -interior(time_series.loss_ts₁[$n], 1, 1, :)*vol_day_correction₁-interior(time_series.loss_ds₁[$n], 1, 1, :)*vol_day_correction₁-interior(time_series.consumption₁[$n], 1, 1, :)*vol_day_correction₁
 growth_ts₂ₙ = @lift interior(time_series.growth_ts₂[$n], 1, 1, :)*vol_day_correction₂
 growth_ds₂ₙ = @lift interior(time_series.growth_ds₂[$n], 1, 1, :)*vol_day_correction₂
 loss_ts₂ₙ = @lift -interior(time_series.loss_ts₂[$n], 1, 1, :)*vol_day_correction₂
 loss_ds₂ₙ = @lift -interior(time_series.loss_ds₂[$n], 1, 1, :)*vol_day_correction₂
 consumption₂ₙ = @lift -interior(time_series.consumption₂[$n], 1, 1, :)*vol_day_correction₂
 total₂ₙ = @lift interior(time_series.growth_ts₂[$n], 1, 1, :)*vol_day_correction₂+interior(time_series.growth_ds₂[$n], 1, 1, :)*vol_day_correction₂-interior(time_series.loss_ts₂[$n], 1, 1, :)*vol_day_correction₂-interior(time_series.loss_ds₂[$n], 1, 1, :)*vol_day_correction₂-interior(time_series.consumption₂[$n], 1, 1, :)*vol_day_correction₂
 growth_ts₃ₙ = @lift interior(time_series.growth_ts₃[$n], 1, 1, :)*vol_day_correction₃
 growth_ds₃ₙ = @lift interior(time_series.growth_ds₃[$n], 1, 1, :)*vol_day_correction₃
 consumption₃ₙ = @lift -interior(time_series.consumption₃[$n], 1, 1, :)*vol_day_correction₃
 total₃ₙ = @lift interior(time_series.growth_ts₃[$n], 1, 1, :)*vol_day_correction₃+interior(time_series.growth_ds₃[$n], 1, 1, :)*vol_day_correction₃-interior(time_series.consumption₃[$n], 1, 1, :)*vol_day_correction₃

# Make figure and plots
fig = Figure(size = (1250, 500))

axis_kwargs_C₁ = (xlabel="dV/dt (ppm per day)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((-3, 7), (-Lz, 0)))

axis_kwargs_C₂ = (xlabel="dV/dt (ppm per day)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((-3, 7), (-Lz, 0)))

axis_kwargs_C₃ = (xlabel="dV/dt (ppm per day)",
                    ylabel="z (m)",
                    aspect = AxisAspect(1),
                    limits = ((-10, 30), (-Lz, 0)))

ax_C₁ = Axis(fig[2, 1]; title = "Small particles, r = $(r₁) m", axis_kwargs_C₁...)
ax_C₂ = Axis(fig[2, 2]; title = "Medium particles, r = $(r₂) m", axis_kwargs_C₂...)
ax_C₃ = Axis(fig[2, 3]; title = "Large particles, r = $(r₃) m", axis_kwargs_C₃...)

title = @lift @sprintf("t = %s", prettytime(times[$n]))
fig[1, 1:3] = Label(fig, title, fontsize=24, tellwidth=false)

loss_ts_line = lines!(ax_C₁, loss_ts₁ₙ, zc, color=:red, linestyle=:dot)
loss_ds_line = lines!(ax_C₁, loss_ds₁ₙ, zc, color=:red, linestyle=:solid)
consumption_line = lines!(ax_C₁, consumption₁ₙ, zc, color=:maroon, linestyle=:solid)
total_line = lines!(ax_C₁, total₁ₙ, zc, color=:black, linestyle=:solid)

growth_ts_line = lines!(ax_C₂, growth_ts₂ₙ, zc, color=:blue, linestyle=:dot)
growth_ds_line = lines!(ax_C₂, growth_ds₂ₙ, zc, color=:blue, linestyle=:solid)
lines!(ax_C₂, loss_ts₂ₙ, zc, color=:red, linestyle=:dot)
lines!(ax_C₂, loss_ds₂ₙ, zc, color=:red, linestyle=:solid)
lines!(ax_C₂, consumption₂ₙ, zc, color=:maroon, linestyle=:solid)
lines!(ax_C₂, total₂ₙ, zc, color=:black, linestyle=:solid)

lines!(ax_C₃, growth_ts₃ₙ, zc, color=:blue, linestyle=:dot)
lines!(ax_C₃, growth_ds₃ₙ, zc, color=:blue, linestyle=:solid)
lines!(ax_C₃, consumption₃ₙ, zc, color=:maroon, linestyle=:solid)
lines!(ax_C₃, total₃ₙ, zc, color=:black, linestyle=:solid)

Legend(fig[2, 4],
    [total_line, growth_ts_line, growth_ds_line, loss_ts_line, loss_ds_line, consumption_line],
    ["Net growth", "γᵗˢ (growth due to turbulence)", "γᵈˢ (growth due to differential settling)", "φᵗˢ (loss due to turbulence)", "φᵈˢ (loss due to differential settling)", "Loss due to consumption and remineralisation"])

# Now record the movie

frames = 1:length(times)

@info "Making a motion picture..."

record(fig, "scripts/1D_growth_and_loss_profiles.mp4", frames, framerate=8) do i
    n[] = i
end
             