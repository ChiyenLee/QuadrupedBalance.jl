module QuadrupedBalance

using RigidBodyDynamics
using RobotDynamics
using ForwardDiff
using Rotations
using StaticArrays
using LinearAlgebra
using SparseArrays
include("rigidbodymodel.jl")
include("quadruped_fullbody.jl")
include("quadruped_flywheel.jl")
include("centroidal_model.jl")
include("forward_kinematics.jl")
include("centroidal_pendulum.jl")
include("simple_sim.jl")
include("utils.jl")
end # module
