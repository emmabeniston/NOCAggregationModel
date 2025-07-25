# This module generates the structs and equations for the NOC Particle Aggregation model
# This is the model used by Williams et al.
# Currently coded for n=3, but may generalise later on

module NOCAggregationModel

export NOCAggregationModelThreeCompartment

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
    r₁ :: A # radius of each particle class
    r₂ :: A
    r₃ :: A
    w_sink_speed₁ :: A # sinking speed (positive scalar) of each particle class
    w_sink_speed₂ :: A
    w_sink_speed₃ :: A    
    w_sink_velocity₁ :: W # sinking velocity field of each particle class
    w_sink_velocity₂ :: W
    w_sink_velocity₃ :: W
    α :: A # stickiness
    D :: A # fractal dimension
    ν₀ :: A # kinematic viscosity of seawater
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
                                                                                                    adapt(to, bgc.α),
                                                                                                    adapt(to, bgc.D),
                                                                                                    adapt(to, bgc.ν₀))

# Define function so that OceanBioME and Oceananigans know the tracers

required_biogeochemical_tracers(::NOCAggregationModelThreeCompartment) = (:C₁, :C₂, :C₃, :turbulent_dissipation) 
# Note the Cs are number concentrations (not volume or mass!) and need turbulent dissipation as tracer so can be included as variable in equations

# Create functions to specify how our tracers evolve (everything except sinking and relaxation to background state)

κ(rᵢ, rⱼ, w_sink_speedᵢ, w_sink_speedⱼ, turbulent_dissipation, ν₀) = 1.3 * sqrt(turbulent_dissipation/ν₀) * (rᵢ+rⱼ)^3 + π * (rᵢ+rⱼ)^2 * abs(w_sink_speedᵢ-w_sink_speedⱼ) # General collision kernel
mass_correction(r_small, r_large, r_new, D) = 1 / ((r_new/r_small)^D - (r_large/r_small)^D)
β(z) = 0.2 * exp(0.17*z/100) / 86400 # Units of 1/seconds (in paper is 1/days)

(bgc::NOCAggregationModelThreeCompartment)(::Val{:C₁}, x, y, z, t, C₁, C₂, C₃, turbulent_dissipation) = - bgc.α * κ(bgc.r₁, bgc.r₁, bgc.w_sink_speed₁, bgc.w_sink_speed₁, turbulent_dissipation, bgc.ν₀) * C₁ * C₁ * mass_correction(bgc.r₁, bgc.r₁, bgc.r₂, bgc.D) - bgc.α * κ(bgc.r₁, bgc.r₁, bgc.w_sink_speed₁, bgc.w_sink_speed₁, turbulent_dissipation, bgc.ν₀) * C₁ * C₁ - bgc.α * κ(bgc.r₂, bgc.r₁, bgc.w_sink_speed₂, bgc.w_sink_speed₁, turbulent_dissipation, bgc.ν₀) * C₂ * C₁ - β(z) * C₁
(bgc::NOCAggregationModelThreeCompartment)(::Val{:C₂}, x, y, z, t, C₁, C₂, C₃, turbulent_dissipation) = + bgc.α * κ(bgc.r₁, bgc.r₁, bgc.w_sink_speed₁, bgc.w_sink_speed₁, turbulent_dissipation, bgc.ν₀) * C₁ * C₁ * mass_correction(bgc.r₁, bgc.r₁, bgc.r₂, bgc.D) - bgc.α * κ(bgc.r₁, bgc.r₂, bgc.w_sink_speed₁, bgc.w_sink_speed₂, turbulent_dissipation, bgc.ν₀) * C₁ * C₂ * mass_correction(bgc.r₁, bgc.r₂, bgc.r₃, bgc.D) - bgc.α * κ(bgc.r₂, bgc.r₂, bgc.w_sink_speed₂, bgc.w_sink_speed₂, turbulent_dissipation, bgc.ν₀) * C₂ * C₂ * mass_correction(bgc.r₂, bgc.r₂, bgc.r₃, bgc.D) - bgc.α * κ(bgc.r₂, bgc.r₂, bgc.w_sink_speed₂, bgc.w_sink_speed₂, turbulent_dissipation, bgc.ν₀) * C₂ * C₂ - β(z) * C₂
(bgc::NOCAggregationModelThreeCompartment)(::Val{:C₃}, x, y, z, t, C₁, C₂, C₃, turbulent_dissipation) = + bgc.α * κ(bgc.r₁, bgc.r₂, bgc.w_sink_speed₁, bgc.w_sink_speed₂, turbulent_dissipation, bgc.ν₀) * C₁ * C₂ * mass_correction(bgc.r₁, bgc.r₂, bgc.r₃, bgc.D) + bgc.α * κ(bgc.r₂, bgc.r₂, bgc.w_sink_speed₂, bgc.w_sink_speed₂, turbulent_dissipation, bgc.ν₀) * C₂ * C₂ * mass_correction(bgc.r₂, bgc.r₂, bgc.r₃, bgc.D) - β(z) * C₃

# For sinking, add a method to the Oceananigans function biogeochemical_drift_velocity

biogeochemical_drift_velocity(bgc::NOCAggregationModelThreeCompartment, ::Val{:C₁}) = (u = ZeroField(), v = ZeroField(), w = bgc.w_sink_velocity₁)
biogeochemical_drift_velocity(bgc::NOCAggregationModelThreeCompartment, ::Val{:C₂}) = (u = ZeroField(), v = ZeroField(), w = bgc.w_sink_velocity₂)
biogeochemical_drift_velocity(bgc::NOCAggregationModelThreeCompartment, ::Val{:C₃}) = (u = ZeroField(), v = ZeroField(), w = bgc.w_sink_velocity₃)

# Make a function so that into NOCAggregationModelThreeCompartment you give the grid and model parameters
# It will convert a sinking speed scalar into a downwards sinking velocity field
# Can later add defaults for the kwargs if wanted

function NOCAggregationModelThreeCompartment(; grid, r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, α, D, ν₀) 
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
    
    return NOCAggregationModelThreeCompartment(r₁, r₂, r₃, w_sink_speed₁, w_sink_speed₂, w_sink_speed₃, w_sink_velocity₁, w_sink_velocity₂, w_sink_velocity₃, α, D, ν₀) 

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
