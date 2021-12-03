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
include("forward_kinematics.jl")
include("utils.jl")
end # module
