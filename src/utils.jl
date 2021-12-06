######################## Interface Utilities ###############
struct MotorIDs 
    FR_Hip::Int; FR_Thigh::Int; FR_Calf::Int
    FL_Hip::Int; FL_Thigh::Int; FL_Calf::Int 
    RR_Hip::Int; RR_Thigh::Int; RR_Calf::Int
    RL_Hip::Int; RL_Thigh::Int; RL_Calf::Int
end   

const MotorIDs_c = MotorIDs(1,2,3,4,5,6,7,8,9,10,11,12)
const MotorIDs_rgb = MotorIDs(1,5,9,2,6,10,3,7,11,4,8,12)

"""Given a motor command/array, map from IDs1 to IDs2"""
function mapMotorArrays(cmd1::AbstractVector{Float64}, ids1::MotorIDs, ids2::MotorIDs)
    cmd2 = @MVector zeros(length(cmd1))
    for motor in fieldnames(MotorIDs) 
        id1 = getfield(ids1,motor)
        id2 = getfield(ids2,motor)
        cmd2[id2] = cmd1[id1]
    end 
    cmd2 = SVector(cmd2)
    return cmd2
end 

########################### Math Utilities ###################
"""Returns the cross product matrix """
function cross_mat(v)
    return [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] -v[1] 0]
end 

"""Given quaternion returns G matrix"""
function get_G(quat)
    H = [zeros(3) I]'
    return L(quat) * H
end 

function quaternion_differential(quat)
    H = [zeros(3) I(3)]' 
    return L(quat) * H
end 

"""Calculate the adjoint matrix of homogenous transform g"""
function adj(g)
    R = g[1:3,1:3]
    p = g[1:3,end]
    return [R zeros(3,3); cross_mat(p)*R R]
    
end 

"""Calculate the inverse adjoint matrix """
function adj_inv(g)
    R = g[1:3,1:3]
    p = g[1:3,end]
    return [R' -R'*cross_mat(p); I R']
end

function L(q)
    s = q[1]
    v = q[2:4]
    L = [s    -v';
         v  s*I+hat(v)]
    return L
end

function hat(v)
    return [0 -v[3] v[2];
            v[3] 0 -v[1];
            -v[2] v[1] 0]
end

"""Quaternion error to axis angle"""
function ζ(v)
    v_norm = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    if v_norm != 0
        vec = sin(0.5*v_norm) * v/v_norm 
        scalar = cos(0.5*v_norm)
    else 
        scalar = 1 
        vec = zeros(3)
    end
    return [scalar;vec]

end

function q_inv(q)
    return [q[1];-q[2:end]]
end



function dynamics_rk4(model::RobotDynamics.AbstractModel, x,u,h)
    #RK4 integration with zero-order hold on u
    f1 = dynamics(model, x, u)
    f2 = dynamics(model, x + 0.5*h*f1, u)
    f3 = dynamics(model, x + 0.5*h*f2, u)
    f4 = dynamics(model, x + h*f3, u)
    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
end

"""Get linearized A, B, C dynamics matrix"""
function dynamics_jacobians_discrete(model::RobotDynamics.AbstractModel, x, u, h)
    A = ForwardDiff.jacobian(t->dynamics_rk4(model,t,u,h), x)
    B = ForwardDiff.jacobian(t->dynamics_rk4(model,x,t,h), u)
    return A, B
end 

"""Get linearized A, B, C dynamics matrix"""
function dynamics_jacobians(model::RobotDynamics.AbstractModel, x, u)
    A = ForwardDiff.jacobian(t->dynamics(model,t,u), x)
    B = ForwardDiff.jacobian(t->dynamics(model,x,t), u)
    return A, B
end 

function dynamics_jacobians(A1, xf, uf, λf)
    A = ForwardDiff.jacobian(t->pinned_dynamics(A1,t,uf, λf), xf)
    B = ForwardDiff.jacobian(t->pinned_dynamics(A1,xf,t, λf), uf)
    C = ForwardDiff.jacobian(t->pinned_dynamics(A1,xf,uf, t), λf)
    return A, B, C
end 

function dynamics_jacobians(A1, xf, uf, λf, foot_indices)
    A = ForwardDiff.jacobian(t->pinned_dynamics(A1,t,uf, λf, foot_indices), xf)
    B = ForwardDiff.jacobian(t->pinned_dynamics(A1,xf,t, λf, foot_indices), uf)
    C = ForwardDiff.jacobian(t->pinned_dynamics(A1,xf,uf, t, foot_indices), λf)
    return A, B, C
end 

    