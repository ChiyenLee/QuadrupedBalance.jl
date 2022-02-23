using Pkg; Pkg.activate(".")
using Revise
using QuadrupedBalance
using TOML
using RobotDynamics
using LinearAlgebra
using OSQP
using SparseArrays
using RigidBodyDynamics
using Rotations
using PyPlot
using MeshCat
using MeshCatMechanisms
using StaticArrays
using DelimitedFiles
const QB = QuadrupedBalance
const RD = RobotDynamics

mutable struct JointController 
    joint_des::Vector{Float64}
    contact_indices::Vector{Int64}
    Kp::Float64 
    Kd::Float64 
end 

function joint_pd_control(q, q_des, q_v, Kp, Kd)
    q_diff = q - q_des 
    u = -Kp .* q_diff - Kd .* q_v 
end 

function bisection_search(sim, controller, x_guess_high, x_guess_low, tf_sim ; joint_noise = 0.0)
    x_foot_high = x_guess_high   # x position for end point 1 (front foot, back foot)
    x_foot_low = x_guess_low # x position for end point 2
    x_foot_now = (x_foot_low + x_foot_high) / 2
    y_foot = 0.1308     # y position. We keep this cnstant 

    # Guesses for inverse kinematics 
    q_guess = [ones(4)*0.0; ones(4)*0.9; ones(4)*-1.5]
    q_stand =QB.inv_kin([x_foot_now[1], -y_foot, -0.2,
                         x_foot_now[1],   y_foot, -0.2, 
                         x_foot_now[2], -y_foot, -0.2,
                         x_foot_now[2],  y_foot, -0.2], q_guess)
    q_stand = q_stand + randn(12) * joint_noise
    x0 = zeros(37); x0[1] = 1; x0[7] = 0.2; x0[8:19] = q_stand
    controller.joint_des = q_stand
    sim.foot_pos_constraints = QB.fk_world(x0)
    x_pos_final = balance_rollout(sim, controller, x0, tf_sim, dt)

    max_iter = 20
    counter = 0
    while(norm(x_pos_final) > 1e-5 && max_iter > counter)
        counter += 1
        x0[8:19] = q_stand
        controller.joint_des = q_stand 
        sim.foot_pos_constraints = QB.fk_world(x0)
        x_pos_final = balance_rollout(sim, controller, x0, tf_sim, dt)

        # if the robot leans forward, that means CoM is in front of the line. vice versa 
        # we shift the line in the direction of the offset 
        q_guess = copy(q_stand)
        if x_pos_final > 0
            x_foot_low = copy(x_foot_now)
            x_foot_now = (x_foot_now + x_foot_high)/2
        else 
            x_foot_high = copy(x_foot_now) 
            x_foot_now = (x_foot_now + x_foot_low)/2
        end 
        q_stand = QB.inv_kin([x_foot_now[1], -y_foot, -0.2,
                              x_foot_now[1],  y_foot, -0.2, 
                              x_foot_now[2], -y_foot, -0.2,
                              x_foot_now[2],  y_foot, -0.2], q_guess)
        q_stand = q_stand + randn(12) * joint_noise
        # println(sim.contact_indices)
        # println(x_foot_now)
    end 

    return x_foot_now
end 

function find_com(x_guess_high, x_guess_low, dt, tf_sim; joint_noise=0.0)
    # Line 1 bisectional search 
    sim = QB.ConstrainedSim(zeros(12), [true, false, false, true], 0.1, dt)
    controller = JointController(zeros(12), sim.contact_indices, 300, 3)
    x_final1 = bisection_search(sim, controller, x_guess_high, x_guess_low, tf_sim; joint_noise=joint_noise)

    # Line 2 bisectional search 
    sim = QB.ConstrainedSim(zeros(12), [false, true, true, false], 0.1, dt)
    controller = JointController(zeros(12), sim.contact_indices, 300, 3)
    x_final2 = bisection_search(sim, controller, x_guess_high, x_guess_low, tf_sim; joint_noise = joint_noise)

    y = 0.1308 
    p_FR = [x_final1[1], -y]
    p_RL = [x_final1[2],  y]

    p_FL = [x_final2[1], y]
    p_RR = [x_final2[2], -y]

    p_com = find_intersection(p_FR, p_RL, p_FL, p_RR)
    return p_com
end 

"""
    Given four points find the intersection on the line segment 
"""
function find_intersection(p1, p2, q1, q2) 
    k1 = (p2[2] - p1[2]) / (p2[1] - p1[1])
    k2 = (q2[2] - q1[2]) / (q2[1] - q1[1])

    A = [-k1 1; 
         -k2 1] 
    b = [p1[2] - k1 * p1[1]; 
         q1[2] - k2 * q1[1]]
    return A \ b 
end 

"""
    Do a simple balance rollout and return the results 
"""
function balance_rollout(sim, joint_controller, x0, tf, dt, record=false)
    times = 0:dt:tf    
    xs = zeros(length(times), length(x0))
    xs[1,:] = copy(x0)
    u = zeros(12)
    for i in 1:length(times)-1
        joint_positions = xs[i,8:19]
        joint_v = xs[i,26:end]
        u_j = joint_pd_control(joint_positions, joint_controller.joint_des, joint_v, 
                               joint_controller.Kp, joint_controller.Kd)
        u[:]= u_j
        xs[i+1,:], _ = QB.semi_implicit_euler(sim, A1, xs[i,:], u)
    end 

    if record
        return xs 
    else
        x_pos_final = xs[end,5]
        return x_pos_final
    end
end 

urdfpath = joinpath(@__DIR__,"..", "src","a1","urdf","a1.urdf")
A1mech = parse_urdf(urdfpath, floating=true, remove_fixed_tree_joints=false)
A1 = QuadrupedBalance.UnitreeA1FullBody(A1mech);

dt = 0.001
tf = 0.3
sim = QB.ConstrainedSim(ϕ_cons, [true, false, false, true], 0.1, dt)
controller = JointController(q_stand, sim.contact_indices, 300, 3)

x_foot_now = [0.130, -0.130] 
x_guess_high = x_foot_now .+ 0.05
x_guess_low = x_foot_now .- 0.05

com_list = []
for i in 1:200
    com_now = find_com(x_guess_high, x_guess_low, dt, tf; joint_noise = 0.01)
    push!(com_list, com_now)
end 
com_mat = hcat(com_list...)
println("standard deviation ", std(com_mat, dims=2))
println("mean ", mean(com_mat, dims=2))
# x_final = bisection_search(sim, controller, x_guess_high, x_guess_low, tf)


# q_guess = [ones(4)*0.0; ones(4)*0.9; ones(4)*-1.5]
# x_foot_now = x_final
# y_foot = 0.1308
# q_stand =QB.inv_kin([x_foot_now[1], -y_foot, -0.2,
#                      x_foot_now[1],   y_foot, -0.2, 
#                      x_foot_now[2], -y_foot, -0.2,
#                      x_foot_now[2],  y_foot, -0.2], q_guess)

# x0 = zeros(37); x0[1] = 1; x0[7] = 0.2; x0[8:19] = q_stand; 
# ϕ_cons = QB.fk_world(x0)
# sim = QB.ConstrainedSim(ϕ_cons, [false, true, true, false], 0.1, dt)
# x_out = balance_rollout(sim, controller, x0, tf, sim.h, true);




# Setting up Visualization 
vis = Visualizer()
cur_path=pwd() 
cd(joinpath(@__DIR__,"..","src", "a1", "urdf"))
mvis = MechanismVisualizer(A1mech, URDFVisuals(urdfpath), vis)
cd(cur_path)

open(vis)

# Running vis 
times = 0:sim.h:tf
q_anim = [x_out[i,1:19] for i in 1:length(times)-1]
animation = Animation(mvis, times[1:50:end-1], q_anim[1:50:end]);
setanimation!(mvis, animation);
