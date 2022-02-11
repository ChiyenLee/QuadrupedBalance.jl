struct CentroidalPendulum <: RobotDynamics.AbstractModel
    J::Matrix
    m::Real 
    foot_pos_world::Vector{Float64} # FR, FL, RR, RL
    foot_pos_body::Vector{Float64}  # FR, FL, RR, RL
    contacts::Vector{Bool}
    function CentroidalPendulum(J::Matrix, m::Real, foot_pos_world::Vector{Float64}, foot_pos_body::Vector{Float64}, contacts::Vector{Bool})
        
        new(J, m, foot_pos_world, foot_pos_body, contacts)
    end
end 

function RobotDynamics.dynamics(model::CentroidalPendulum, x, u)
    g = [0.0, 0.0, 9.81]
    p = x[1:3]
    ṗ = x[4:6]
    quat = x[7:10]
    ω = x[11:13]
    ρ = x[14:15]   # angular momentum of the fly wheel 

    q = [p;quat]
    q̇ = [ṗ ; ω]

    u_f = [1.0 0.0;
           0.0 1.0;  
           0.0 0.0] * u 
    
    ρ_all = [1.0 0.0;
             0.0 1.0;
             0.0 0.0] * ρ

    # mass matrix 
    M = [model.m*I(3) zeros(3,3);
         zeros(3,3) model.J]
    
    # dynamic bias (coriolis + gravity) 
    cor = [g..., (hat(ω) * (model.J * ω + ρ_all))...]
    
    # constraint dynamics
    J = dcdq(model, q)
    d = ddcdq(model, q,q̇) * q̇
    A = [M J'; 
         J I(6)*1e-8]
    res = [-cor - [0,0,0, u_f...]; 
           -d]
    out = A \ res 
    q̈ = out[1:6]

    p̈ = q̈[1:3]
    ω_dot = q̈[4:6]

    # kinematics
    quat_dot = 0.5*L(quat)*H*ω
    ρ_dot = u
    

    return [ṗ; p̈; quat_dot; ω_dot; ρ_dot]

end

function constraints(model::CentroidalPendulum, q::AbstractVector{T}) where T
    c = Array{T,1}()
    p = q[1:3]
    quat = q[4:7]
    for i in 1:length(model.contacts)
        if model.contacts[i]
            r_w = model.foot_pos_world[(i-1)*3+1:(i-1)*3+3]
            r_b = model.foot_pos_body[(i-1)*3+1:(i-1)*3+3]
            append!(c, UnitQuaternion(quat) * r_b + p - r_w)
        end 
    end 

    return c 
end 

function dcdq(model::CentroidalPendulum, q)
    attitude_error_jacoian = SparseArrays.blockdiag(sparse(I(3)), sparse(0.5*quaternion_differential(q[4:7])))
    return ForwardDiff.jacobian(t->constraints(model,t), q) * attitude_error_jacoian
end 

"""d(J(q)*q̇)"""
function ddcdq(model::CentroidalPendulum, q, q̇)
    attitude_error_jacoian = SparseArrays.blockdiag(sparse(I(3)), sparse(0.5*quaternion_differential(q[4:7])))
    f(model, q) = dcdq(model, q) * q̇
    return ForwardDiff.jacobian(t->f(model,t), q) * attitude_error_jacoian
end 

function linearize(model, x, u)
    A = ForwardDiff.jacobian(t->RobotDynamics.dynamics(model, t, u), x)
    B = ForwardDiff.jacobian(t->RobotDynamics.dynamics(model, x, t), u)

    return A, B
end 