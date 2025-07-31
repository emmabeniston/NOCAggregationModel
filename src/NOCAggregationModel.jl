# This module generates the structs and equations for the NOC Particle Aggregation model
# This is the model used by Williams et al.
# Currently coded for n=3, but may generalise later on

module NOCAggregationModel

# Export the BGC model
export NOCAggregationModelThreeCompartment

# Export the aggregation and consumption functions that are used in the tendency
export consumption_func
export growth_ts₁_func, growth_ds₁_func, loss_ts₁_func, loss_ds₁_func
export growth_ts₂_func, growth_ds₂_func, loss_ts₂_func, loss_ds₂_func
export growth_ts₃_func, growth_ds₃_func, loss_ts₃_func, loss_ds₃_func

using Oceananigans
using OceanBioME
using Oceananigans.Units
using CUDA
using Adapt

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

using Oceananigans.Fields: ZeroField, ConstantField

import Adapt: adapt_structure
import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

# ## Create NOCAggregationModelThreeCompartment biogeochemical model

# Create biogeochemistry

struct NOCAggregationModelThreeCompartment{A, W} <: AbstractContinuousFormBiogeochemistry
    r₁ :: A # radius of each particle class, in m
    r₂ :: A
    r₃ :: A
    w_sink_speed₁ :: A # sinking speed (positive scalar) of each particle class, in m s⁻¹
    w_sink_speed₂ :: A
    w_sink_speed₃ :: A    
    w_sink_velocity₁ :: W # sinking velocity field of each particle class, in m s⁻¹
    w_sink_velocity₂ :: W
    w_sink_velocity₃ :: W
    stick :: A # stickiness
    D :: A # fractal dimension
    ν₀ :: A # kinematic viscosity of seawater, in m² s⁻¹
end

# Adapt so can run on GPU

adapt_structure(to, bgc::NOCAggregationModelThreeCompartment) = NOCAggregationModelThreeCompartment(adapt(to, bgc.r₁),
                                                                                                    adapt(to, bgc.r₂),
                                                                                                    adapt(to, bgc.r₃),
                                                                                                    adapt(to, bgc.w_sink_speed₁),
                                                                                                    adapt(to, bgc.w_sink_speed₂),
                                                                                                    adapt(to, bgc.w_sink_speed₃),
                                                                                                    adapt(to, bgc.w_sink_velocity₁),
                                                                                                    adapt(to, bgc.w_sink_velocity₂),
                                                                                                    adapt(to, bgc.w_sink_velocity₃),
                                                                                                    adapt(to, bgc.stick),
                                                                                                    adapt(to, bgc.D),
                                                                                                    adapt(to, bgc.ν₀))

# Define function so that OceanBioME and Oceananigans know the tracers

required_biogeochemical_tracers(::NOCAggregationModelThreeCompartment) = (:C₁, :C₂, :C₃, :turbulent_dissipation) 
# Note the Cs are number concentrations (not volume or mass!) and need turbulent dissipation as tracer so can be included as variable in equations

# ## Create functions to specify how our tracers evolve (everything except sinking and relaxation to background state)

# Collision kernels, mass correction, and decay rate

κ_ts(rᵢ, rⱼ, w_sink_speedᵢ, w_sink_speedⱼ, turbulent_dissipation, ν₀) = 1.3 * sqrt(turbulent_dissipation/ν₀) * (rᵢ+rⱼ)^3 # Collision kernel due to turbulent shear
κ_ds(rᵢ, rⱼ, w_sink_speedᵢ, w_sink_speedⱼ, turbulent_dissipation, ν₀) = π * (rᵢ+rⱼ)^2 * abs(w_sink_speedᵢ-w_sink_speedⱼ) # Collision kernel due to differential settling
κ(rᵢ, rⱼ, w_sink_speedᵢ, w_sink_speedⱼ, turbulent_dissipation, ν₀) = κ_ts(rᵢ, rⱼ, w_sink_speedᵢ, w_sink_speedⱼ, turbulent_dissipation, ν₀) + κ_ds(rᵢ, rⱼ, w_sink_speedᵢ, w_sink_speedⱼ, turbulent_dissipation, ν₀) 
mass_correction(r_small, r_large, r_new, D) = 1 / ((r_new/r_small)^D - (r_large/r_small)^D) # Mass correction as need multiple small particles to make a larger one
β(z) = 0.2 * exp(0.17*z/100) / 86400 # Decay rate due to consumption and remineralisation, units of 1/seconds (in paper is 1/days, hence 86400 correction)
consumption_func(z, C) = β(z) * C

# Terms for C₁ - start with functions for each term, which can be exported for output and analysis, and added together for the tendency term - all terms have positive sign so need to add/subtract appropriately

growth_ts₁_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = 0 
growth_ds₁_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = 0
loss_ts₁_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ts(r₁, r₁, w_sink_speed₁, w_sink_speed₁, turbulent_dissipation, ν₀) * C₁ * C₁ * mass_correction(r₁, r₁, r₂, D) + stick * κ_ts(r₁, r₁, w_sink_speed₁, w_sink_speed₁, turbulent_dissipation, ν₀) * C₁ * C₁ + stick * κ_ts(r₂, r₁, w_sink_speed₂, w_sink_speed₁, turbulent_dissipation, ν₀) * C₂ * C₁
loss_ds₁_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ds(r₁, r₁, w_sink_speed₁, w_sink_speed₁, turbulent_dissipation, ν₀) * C₁ * C₁ * mass_correction(r₁, r₁, r₂, D) + stick * κ_ds(r₁, r₁, w_sink_speed₁, w_sink_speed₁, turbulent_dissipation, ν₀) * C₁ * C₁ + stick * κ_ds(r₂, r₁, w_sink_speed₂, w_sink_speed₁, turbulent_dissipation, ν₀) * C₂ * C₁
(bgc::NOCAggregationModelThreeCompartment)(::Val{:C₁}, x, y, z, t, C₁, C₂, C₃, turbulent_dissipation) = growth_ts₁_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) + growth_ds₁_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - loss_ts₁_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - loss_ds₁_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - consumption_func(z, C₁)

# Terms for C₂

growth_ts₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ts(r₁, r₁, w_sink_speed₁, w_sink_speed₁, turbulent_dissipation, ν₀) * C₁ * C₁ * mass_correction(r₁, r₁, r₂, D) 
growth_ds₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ds(r₁, r₁, w_sink_speed₁, w_sink_speed₁, turbulent_dissipation, ν₀) * C₁ * C₁ * mass_correction(r₁, r₁, r₂, D)
loss_ts₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ts(r₁, r₂, w_sink_speed₁, w_sink_speed₂, turbulent_dissipation, ν₀) * C₁ * C₂ * mass_correction(r₁, r₂, r₃, D) + stick * κ_ts(r₂, r₂, w_sink_speed₂, w_sink_speed₂, turbulent_dissipation, ν₀) * C₂ * C₂ * mass_correction(r₂, r₂, r₃, D) + stick * κ_ts(r₂, r₂, w_sink_speed₂, w_sink_speed₂, turbulent_dissipation, ν₀) * C₂ * C₂
loss_ds₂_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ds(r₁, r₂, w_sink_speed₁, w_sink_speed₂, turbulent_dissipation, ν₀) * C₁ * C₂ * mass_correction(r₁, r₂, r₃, D) + stick * κ_ds(r₂, r₂, w_sink_speed₂, w_sink_speed₂, turbulent_dissipation, ν₀) * C₂ * C₂ * mass_correction(r₂, r₂, r₃, D) + stick * κ_ds(r₂, r₂, w_sink_speed₂, w_sink_speed₂, turbulent_dissipation, ν₀) * C₂ * C₂
(bgc::NOCAggregationModelThreeCompartment)(::Val{:C₂}, x, y, z, t, C₁, C₂, C₃, turbulent_dissipation) = growth_ts₂_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) + growth_ds₂_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - loss_ts₂_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - loss_ds₂_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - consumption_func(z, C₂)

# Terms for C₃

growth_ts₃_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ts(r₁, r₂, w_sink_speed₁, w_sink_speed₂, turbulent_dissipation, ν₀) * C₁ * C₂ * mass_correction(r₁, r₂, r₃, D) + stick * κ_ts(r₂, r₂, w_sink_speed₂, w_sink_speed₂, turbulent_dissipation, ν₀) * C₂ * C₂ * mass_correction(r₂, r₂, r₃, D)
growth_ds₃_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = stick * κ_ds(r₁, r₂, w_sink_speed₁, w_sink_speed₂, turbulent_dissipation, ν₀) * C₁ * C₂ * mass_correction(r₁, r₂, r₃, D) + stick * κ_ds(r₂, r₂, w_sink_speed₂, w_sink_speed₂, turbulent_dissipation, ν₀) * C₂ * C₂ * mass_correction(r₂, r₂, r₃, D)
loss_ts₃_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = 0
loss_ds₃_func(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀, C₁, C₂, C₃, turbulent_dissipation) = 0
(bgc::NOCAggregationModelThreeCompartment)(::Val{:C₃}, x, y, z, t, C₁, C₂, C₃, turbulent_dissipation) = growth_ts₃_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) + growth_ds₃_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - loss_ts₃_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - loss_ds₃_func(bgc.r₁, bgc.r₂, bgc.r₃, bgc.w_sink_speed₁, bgc.w_sink_speed₂, bgc.w_sink_speed₃, bgc.stick, bgc.D, bgc.ν₀, C₁, C₂, C₃, turbulent_dissipation) - consumption_func(z, C₃)

# For sinking, add a method to the Oceananigans function biogeochemical_drift_velocity

biogeochemical_drift_velocity(bgc::NOCAggregationModelThreeCompartment, ::Val{:C₁}) = (u = ZeroField(), v = ZeroField(), w = bgc.w_sink_velocity₁)
biogeochemical_drift_velocity(bgc::NOCAggregationModelThreeCompartment, ::Val{:C₂}) = (u = ZeroField(), v = ZeroField(), w = bgc.w_sink_velocity₂)
biogeochemical_drift_velocity(bgc::NOCAggregationModelThreeCompartment, ::Val{:C₃}) = (u = ZeroField(), v = ZeroField(), w = bgc.w_sink_velocity₃)

# Make a function so that into NOCAggregationModelThreeCompartment you give the grid and model parameters
# It will convert a sinking speed scalar into a downwards sinking velocity field
# Can later add defaults for the kwargs if wanted

function NOCAggregationModelThreeCompartment(; grid, r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, stick, D, ν₀) 
    w_sink_velocity₁ = ZFaceField(grid)
    w_sink_velocity_func₁(x, y, z) = z == 0 ? 0 : -w_sink_speed₁
    w_sink_velocity_func₁(z) = z == 0 ? 0 : -w_sink_speed₁ # set two methods in case we have a column model
    set!(w_sink_velocity₁, w_sink_velocity_func₁)

    w_sink_velocity₂ = ZFaceField(grid)
    w_sink_velocity_func₂(x, y, z) = z == 0 ? 0 : -w_sink_speed₂
    w_sink_velocity_func₂(z) = z == 0 ? 0 : -w_sink_speed₂ # set two methods in case we have a column model
    set!(w_sink_velocity₂, w_sink_velocity_func₂)

    w_sink_velocity₃ = ZFaceField(grid)
    w_sink_velocity_func₃(x, y, z) = z == 0 ? 0 : -w_sink_speed₃
    w_sink_velocity_func₃(z) = z == 0 ? 0 : -w_sink_speed₃ # set two methods in case we have a column model
    set!(w_sink_velocity₃, w_sink_velocity_func₃)
    
    return NOCAggregationModelThreeCompartment(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, w_sink_velocity₁, w_sink_velocity₂, w_sink_velocity₃, stick, D, ν₀) 

end

# In main model script, will have to add in the callback to compute the turbulent_dissipation
# The code for that is commented out here:

# # Add callback to compute the turbulent dissipation rate:

# using Oceanostics

# turbulent_dissipation_op = KineticEnergyDissipationRate(model)
# turbulent_dissipation = model.tracers.turbulent_dissipation

# function compute_turbulent_dissipation!(sim)
#     turbulent_dissipation .= turbulent_dissipation_op
#     return nothing
# end

# simulation.callbacks[:turbulent_dissipation] = Callback(compute_turbulent_dissipation!)

end # module NOCAggregationModel
