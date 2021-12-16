struct CentroidalModel <: RobotDynamics.AbstractModel
    J::Matrix
    m::Real 
    p1::Vector # FR 
    p2::Vector # FL 
    p3::Vector # RR 
    p4::Vector # RL
end 

function RobotDynamics.dynamics(model::CentroidalModel, x, u)
    g = [0.0, 0.0, -9.81]
    p = x[1:3]
    ṗ = x[4:6]
    q = x[7:10]
    ω = x[11:13]

    Q = H'L(q)*R(q)'H #Rotation matrix
    r1 = (model.p1 - p)
    r2 = (model.p2 - p)
    r3 = (model.p3 - p)
    r4 = (model.p4 - p)
    M = [I(3) I(3) I(3) I(3);
         hat(r1) hat(r2) hat(r3) hat(r4)]
    u_out = M * u
    f = u_out[1:3]
    τ = u_out[4:6]

    # linear dynamics
    p̈ = 1/model.m * f + g
    
    # attitude dynamics
    damping = 0.0  
    q̇ = 0.5*L(q)*H*ω
    ω_dot =  model.J \ (Q'*τ - damping *ω - hat(ω) * model.J * ω)

    return [ṗ; p̈; q̇; ω_dot]

end

function linearize(model, x, u)
    A = ForwardDiff.jacobian(t->RobotDynamics.dynamics(model, t, u), x)
    B = ForwardDiff.jacobian(t->RobotDynamics.dynamics(model, x, t), u)

    return A, B
end 