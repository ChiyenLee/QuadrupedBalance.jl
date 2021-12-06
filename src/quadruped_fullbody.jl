abstract type AbstractQuadruped <: RobotDynamics.AbstractModel end 

"""
    Full body quadruped model
    x = [q; v]
    q  = [attitude
          position
          hip_angle
          thigh_angle
          calf_angle]
    
    v = [ang_vel
         v_trans (body frame)
         hip_angle
         thigh_angle
         calf_angle]
"""

struct UnitreeA1FullBody <: AbstractQuadruped
    rigidbody::RigidBodyModel 

    function UnitreeA1FullBody(mech::Mechanism) 
        control_indices = Vector{Bool}([zeros(6);ones(12)])
        model = RigidBodyModel(mech, control_indices)
        new(model)
    end 
end 

function dynamics(model::AbstractQuadruped, x::AbstractVector{T1}, u::AbstractVector{T2}) where {T1, T2}
    return RobotDynamics.dynamics(model.rigidbody, x, u)
end 

function get_mass_matrix(model::AbstractQuadruped, x) 
    return get_mass_matrix(model.rigidbody, x)
end 

function get_dynamics_bias(model::AbstractQuadruped, x)
    return get_dynamics_bias(model.rigidbody, x)
end 

RobotDynamics.state_dim(model::AbstractQuadruped) = RobotDynamics.state_dim(model.rigidbody)
RobotDynamics.control_dim(model::AbstractQuadruped) = RobotDynamics.control_dim(model.rigidbody)

function pinned_dynamics(A1::QuadrupedBalance.AbstractQuadruped, x::AbstractVector{T1}, u::AbstractVector{T2}, λ::AbstractVector{T3}, foot_indices) where {T1, T2, T3}
    T = promote_type(typeof(x), typeof(u), typeof(λ))
    x = convert(T, x)
    u = convert(T, u)
    λ = convert(T, λ)
    T = promote_type(T1, T2, T3)
    # attitude_error_jacobian = blockdiag(sparse(QuadrupedBalance.quaternion_differential(x[1:4])), sparse(I(33)) )
    attitude_error_jacobian = blockdiag(sparse(0.5*QuadrupedBalance.quaternion_differential(x[1:4])),
                                        sparse(Rotations.UnitQuaternion(x[1:4])),
                                        sparse(I(30)) )
    J = QuadrupedBalance.dfk_world(x)[foot_indices,:] * attitude_error_jacobian
    J = J[:, 1:18]
    # J = contacts_jacobian(A1, x)[foot_indices,:]
    M = QuadrupedBalance.get_mass_matrix(A1, x)
    ẋ = QuadrupedBalance.dynamics(A1, x, u)    
    ẋ[20:end] .= ẋ[20:end] .+ inv(M) * J' * λ
    return ẋ
end 