struct UnitreeA1FlyWheel <: AbstractQuadruped
    rigidbody::RigidBodyModel

    function UnitreeA1FlyWheel(mech::Mechanism)
        control_indices = Vector{Bool}([zeros(6);ones(14)])
        model = RigidBodyModel(mech, control_indices)
        new(model)
    end 
end 

function RobotDynamics.dynamics(model::UnitreeA1FlyWheel, x::AbstractVector{T1}, u::AbstractVector{T2}) where {T1, T2}
    return RobotDynamics.dynamics(model.rigidbody, x, u)
end 

# RobotDynamics.state_dim(model::AbstractQuadruped) = RobotDynamics.state_dim(model.rigidbody)
# RobotDynamics.control_dim(model::AbstractQuadruped) = RobotDynamics.control_dim(model.rigidbody)
