using StaticArrays 
using OSQP
using ForwardDiff 
using Dojo
using MeshCat 
using QuadrupedBalance
const QB = QuadrupedBalance 
include("dojo_utils.jl")
mutable struct BalanceControl 
    x_des::Vector # desired state 
    f_des::Vector # desired control vector 
    A::Matrix     # Linearized state transition jacobian
    B::Matrix     # Linearized control jacobian
    Q::Matrix     # LQR gain
    R::Matrix     # LQR regularization
    P::Matrix
    h::Float64    # time step size     
end 

""" 
    Given two foot points on the ground. Calcualte reaction in an equilibrium point 
"""
function calc_eq_forces(model, x)    
    ## Matrices that map grnd reaction to body frame
    rot = QB.H'QB.L(x[7:10])*QB.R(x[7:10])'QB.H #Rotation matrix
    r̂1 = rot'*QB.hat((model.p1-x[1:3])) # rotate to body frame
    r̂2 = rot'*QB.hat((model.p2-x[1:3]))
    r̂3 = rot'*QB.hat((model.p3-x[1:3]))
    r̂4 = rot'*QB.hat((model.p4-x[1:3]))
    M = [I(3) I(3) I(3) I(3);
        r̂1   r̂2    r̂3   r̂4] 
    g = [0, 0, -9.81]
    b = [model.m * g; zeros(3)] # Desired translation + rotational acceleration
    
    ## QP matrix and vector weightd 
    W = Diagonal(ones(6))
    P = M' * W * M + 1e-5*I(12)
    q = b' * W * M
    
    problem = OSQP.Model()
    OSQP.setup!(problem, P=P, q=q', A=sparse(consMatrix), l=lb, u=ub, verbose=0)
    results = OSQP.solve!(problem)
    return results.x
end 

"""
    Backward ricatti 
"""
function ricatti(A, B, Q, R)
    P = copy(Q)
    P_prev = zero(P)
    K = zeros(size(B'))
    res = norm(P-P_prev)
    while res > 1e-7
        K = (R + B'*P_prev*B)\ B'*P*A
        P = Q + A'*P_prev*(A-B*K)
        res = norm(P-P_prev)
        P_prev = copy(P)
    end
    return P, K
end 

"""
    Dojo controller interface 
"""
function balance_controller(mechanism, k)
    x = get_minimal_state(mechanism)
    ind_flipped = SVector{12,Int64}([10:12;7:9;4:6;1:3]) # for some reason the indices are flipped 
    ind_foot_flipped = SVector{4,Int64}([4,3,2,1])
    encs = x[13:end][1:2:end][ind_flipped]
    pos = x[1:3] + [-0.01, 0.0, 0.0]

    # controller.x_des[7:10] = QB.ρ([0., 0., y_des[k]])

    # error calculation 
    δx = zeros(12)
    δx[1:3] = pos - controller.x_des[1:3] # position
    δx[4:6] = x[7:9]                    # velocity 

    # δx[7:9] = x[4:6]                    # attitude 
    rot_vec = RotationVec(x[4:6]...)
    diff = Rotations.rotation_error(rot_vec, UnitQuaternion(controller.x_des[7:10]...),
                                    Rotations.CayleyMap())
    δx[7:9] = -diff               
              


    # println(δx[1:3])

    A = controller.A
    B = controller.B
    Q = controller.Q
    R = controller.R
    P_uctg = controller.P

    W = sparse(R + B' * P_uctg * B )
    q = (B' * P_uctg * A * δx)'
    OSQP.setup!(problem, P=W, q=q', A=sparse(consMatrix), l=lb, u=ub, verbose=0, eps_abs=1e-6, eps_rel=1e-6)
    # OSQP.setup!(problem, P=W, q=q', verbose=0, eps_abs=1e-6, eps_rel=1e-6)
    results = OSQP.solve!(problem)
    f = results.x  
    
    rot = RotationVec(x[4:6]...)
    rot = kron(I(4), rot')
    u = -QB.dfk(dojo_to_qb(encs))[:,:]' * rot * (f + controller.f_des) 
    u = qb_to_dojo(u)
    q_stand_dojo = qb_to_dojo(q_stand)

    # extract joints
    leg_joints = mechanism.joints[2:end]
    for i in 1:12
        if contacts[ind_foot_flipped][(i-1) ÷3 + 1] 
            set_input!(leg_joints[i], [u[i]]*mechanism.timestep)
        else
            ui = -200 * (encs[i] - q_stand_dojo[i])
            set_input!(leg_joints[i], [ui*mechanism.timestep])
        end 
    end 
    
    println(δx)
end 
 
## Foot placements 
robot_height = 0.2 
p_FR = [0.1208, -0.1308, 0.0]
p_FL = [0.1208,  0.1308, 0.0]
p_RR = [-0.1408, -0.1308, 0.0]
p_RL = [-0.1408,  0.1308, 0.0]
q_guess = [ones(4)*0.0; ones(4)*0.9; ones(4)*-1.5]
q_stand = QB.inv_kin([p_FR[1],  p_FR[2], -robot_height,
                      p_FL[1],  p_FL[2], -robot_height, 
                      p_RR[1],  p_RR[2], -robot_height,
                      p_RL[1],  p_RL[2], -robot_height], q_guess)


## Initialize model
I_inertia = Diagonal([0.1,0.25,0.3])
model = QB.CentroidalModel(I_inertia, m, p_FR, p_FL, p_RR, p_RL)

## Initialize state and calculate equilibrium pose 
## In this example we will balance on the p_FR + p_RL feet 
x_des = zeros(13)
x_des[7:10] = [1.0, 0.0, 0.0, 0.0]
x_des[1:3] = (p_FR + p_RL) / 2
x_des[3] = robot_height + 0.02 

## Set up the constraint matrix and compute the equilibrium foot forces 
F_MIN = 5
F_MAX = 200
num_legs = 4
μ = 0.7
contacts = [true, false, false, true]
consMatrix = zeros(num_legs + 4*num_legs, 3*num_legs) # constraint matrix
lb = zeros(num_legs + 4*num_legs)
ub = zeros(num_legs + 4*num_legs) 
for i in 1:4
    # establish lowerbound on z 
    consMatrix[i,3+(i-1)*3] = 1 
    c_flag = contacts[i] ? 1.0 : 0.0
    lb[i] = c_flag * F_MIN 
    ub[i] = c_flag * F_MAX

    consMatrix[(num_legs+1) + (i-1)*4, 1+(i-1)*3] = 1
    consMatrix[(num_legs+1) + (i-1)*4, 3+(i-1)*3] = -μ
    lb[(num_legs+1) + (i-1)*4] = -OSQP.OSQP_INFTY

    consMatrix[(num_legs+1) + (i-1)*4 + 1, 1+(i-1)*3] = -1
    consMatrix[(num_legs+1) + (i-1)*4 + 1, 3+(i-1)*3] = -μ
    lb[(num_legs+1) + (i-1)*4+1] = -OSQP.OSQP_INFTY

    consMatrix[(num_legs+1) + (i-1)*4 + 2, 2+(i-1)*3] = 1
    consMatrix[(num_legs+1) + (i-1)*4 + 2, 3+(i-1)*3] = -μ
    lb[(num_legs+1) + (i-1)*4+2] = -OSQP.OSQP_INFTY

    consMatrix[(num_legs+1) + (i-1)*4 + 3, 2+(i-1)*3] = -1
    consMatrix[(num_legs+1) + (i-1)*4 + 3, 3+(i-1)*3] = -μ
    lb[(num_legs+1) + (i-1)*4+3] = -OSQP.OSQP_INFTY
end 
f_des = zeros(12) # since we have twelve grf forces we assigned only to the FR and RL leg 
f_des = calc_eq_forces(model, x_des)

## Linearize about the equilibrium pose and apply the attitude error jacobian for quaternion operation
h = 0.001 # discretization step
A = ForwardDiff.jacobian(t->dynamics(model,t, f_des), x_des)
B = ForwardDiff.jacobian(t->dynamics(model,x_des, t), f_des)
Ad = A*h + I # super coarse approximate discretizatoin step 
Bd = B*h 
G = blockdiag(sparse(1.0I(6)), 
              sparse(QB.quaternion_differential(x_des[7:10])), 
              sparse(1.0I(3)))
Ad = G' * Ad * G 
Bd = G' * Bd

## Update the lowerbound and upperbound to include equilibrium values 
for i in 1:4
    c_flag = contacts[i] ? 1.0 : 0.0
    lb[i] = c_flag * (F_MIN - f_des[(i-1)*3+3]) 
    ub[i] = c_flag * (F_MAX - f_des[(i-1)*3+3])

    lb[(num_legs+1) + (i-1)*4] = -OSQP.OSQP_INFTY
    ub[(num_legs+1) + (i-1)*4] = μ*f_des[(i-1)*3+3] - f_des[(i-1)*3+1]

    lb[(num_legs+1) + (i-1)*4 + 1] = -OSQP.OSQP_INFTY
    ub[(num_legs+1) + (i-1)*4 + 1] = μ*f_des[(i-1)*3+3] + f_des[(i-1)*3+1]

    lb[(num_legs+1) + (i-1)*4 + 2] = -OSQP.OSQP_INFTY
    ub[(num_legs+1) + (i-1)*4 + 2] = μ*f_des[(i-1)*3+3] - f_des[(i-1)*3+2]

    lb[(num_legs+1) + (i-1)*4 + 3] = -OSQP.OSQP_INFTY
    ub[(num_legs+1) + (i-1)*4 + 3] = μ*f_des[(i-1)*3+3] + f_des[(i-1)*3+2]

end 

## Initialize the controller 
Q = Matrix(Diagonal([20, 20, 10000., 1.0, 1.0, 10., 2600, 2800, 1500., 10., 5, 1.])) *1e-2
R = Matrix(Diagonal(kron(ones(4), [0.5, 0.5, .5]))) * 1e-4 
P_uctg, K = ricatti(Ad,Bd,Q,R) # unconstrained cost to go 
controller = BalanceControl(x_des, f_des, Ad, Bd, Q, R, P_uctg, h)

## initialize simulator 
problem = OSQP.Model()
q_stand[10] -= 0.3
q_stand[11] -= 0.3

vis = Visualizer()
open(vis)

dt = 0.01 
tf = 1.5
times = 0:dt:tf

# generate time sinusoidal trajectory 
y_des = sin.(times ./ 0.2) * 0.15
controller.x_des[3] = robot_height + 0.02
controller.x_des[7:10] = [1.0,0.0,0.0,0.0]
mechanism = Dojo.get_mechanism(:quadruped, damper = 1.5,friction_coefficient=2.0, timestep=dt, gravity=-9.81);
set_state(mechanism, q_stand, [-0.01,0,0.2+0.02], [-0.0, 0.0, 0.0])
storage = simulate!(mechanism, tf, balance_controller, record=true, verbose=false);

Dojo.visualize(mechanism, storage, vis=vis);

