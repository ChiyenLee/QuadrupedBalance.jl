""" 
    RigidBodyModel{T}

A RigidBody model as described by a mechanism from RigidBodyDynamics.jl

# Constructor 
    RigidBodyModel(mech)

where 'mech' is a mechanism created from RigidBodyDynamics.jl. 
""" 

struct RigidBodyModel{T} <: RobotDynamics.AbstractModel 
    mech::Mechanism{Float64}
    statecache::T
    dyncache::DynamicsResultCache{Float64}
    control_indices::Vector{Bool}
end 

function RigidBodyModel(mech::Mechanism, control_indices::Vector{Bool})
    statecache = StateCache(mech)
    dyncache = DynamicsResultCache(mech)
    return RigidBodyModel{typeof(statecache)}(mech, statecache, dyncache, control_indices)
end 

function RigidBodyModel(mech::Mechanism)
    statecache = StateCache(mech)
    dyncache = DynamicsResultCache(mech)
    control_indices = ones(Bool, num_positions(mech))
    return RigidBodyModel{typeof(statecache)}(mech, statecache, dyncache, control_indices)
end 

RobotDynamics.state_dim(model::RigidBodyModel) = num_positions(model.mech) + num_velocities(model.mech)
RobotDynamics.control_dim(model::RigidBodyModel) = sum(model.control_indices)

function RobotDynamics.dynamics(model::RigidBodyModel, x::AbstractVector{T1}, u::AbstractVector{T2}) where {T1, T2} 
    T = promote_type(T1, T2)
    state = model.statecache[T]
    res = model.dyncache[T]

    copyto!(state, x)
    τ = zeros(T, num_velocities(model.mech))
    τ[model.control_indices] .= u 
    dynamics!(res, state, τ)
    q̇ = res.q̇ 
    v̇ = res.v̇ 
    return [q̇;v̇]
end 

function get_mass_matrix(model::RigidBodyModel, x::AbstractVector{T}) where T 
    s = model.statecache[T]
    copyto!(s,x)
    return mass_matrix(s) 
end 

function get_dynamics_bias(model::RigidBodyModel, x::Vector{Float64})
    s = model.statecache[Float64]
    copyto!(s,x)
    return dynamics_bias(s)
end 