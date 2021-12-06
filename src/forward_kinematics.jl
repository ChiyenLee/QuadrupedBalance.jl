const x_hip = 0.183 # hip dispalcement from the trunk frame 
const y_hip = 0.047 
const Δy_thigh = 0.08505 # thigh displacement from the hip frame 
const l_limb = 0.2 

"""
    fk(q::AbstractVector{T})
        Given joint angles compute the position of feet relative to obody rigin  
        p = [x_FR, y_FR, z_FR, x_FL, y_FL, z_FL, ...]
        feet order: FR, FL, RR, RL 
"""
function fk(q::AbstractVector)
    q_FR = q[[1,5,9]]
    q_FL = q[[2,6,10]]
    q_RR = q[[3,7,11]]
    q_RL = q[[4,8,12]]


    fk(x_mir, y_mir, θ) = [-l_limb * sin(θ[2] + θ[3]) - l_limb * sin(θ[2]) + x_hip * x_mir; 
                           y_hip * y_mir + Δy_thigh* y_mir * cos(θ[1]) + l_limb*sin(θ[1])*cos(θ[2]+θ[3]) + l_limb*sin(θ[1])*cos(θ[2]); 
                           Δy_thigh * y_mir *  sin(θ[1]) - l_limb * cos(θ[1]) * cos(θ[2] + θ[3]) - l_limb * cos(θ[1])*cos(θ[2]) ] 

    p = [fk(1, -1, q_FR)...;
         fk(1, 1, q_FL)...; 
         fk(-1,-1, q_RR)...; 
         fk(-1,1, q_RL)...]                   
    return p
end

"""
    fk_world(q::AbstractVector{T})
        Given full state vectors (dim=37) compute the position of feet relative to obody rigin  
        p = [x_FR, y_FR, z_FR, x_FL, y_FL, z_FL, ...]
        feet order: FR, FL, RR, RL 
"""
function fk_world(x::AbstractVector)
    q = x[8:19]
    quat = Rotations.UnitQuaternion(x[1:4])
    pos = x[5:7]

    p_body = fk(q)
    p_world = copy(p_body)

    for i in 1:4 
        p_world[(i-1)*3 .+ (1:3)] = quat * p_body[(i-1)*3 .+ (1:3)] + pos 
    end 
    return p_world 
end 

"""
    dfk(q::AbstractVector{T})
        Given joint angles (dim=12) compute the (12 x 12) kinematic jacobian
        p = [x_FR, y_FR, z_FR, x_FL, y_FL, z_FL, ...]
        feet order: FR, FL, RR, RL 
"""
function dfk(q::AbstractVector)
    return ForwardDiff.jacobian(t->fk(t), q)
end 

"""
    fk_world(q::AbstractVector{T})
        Given full state vectors (dim=37) compute (12 x 37) the kinematic jacobian
        p = [x_FR, y_FR, z_FR, x_FL, y_FL, z_FL, ...]
        feet order: FR, FL, RR, RL 
"""
function dfk_world(x::AbstractVector)
    return ForwardDiff.jacobian(t->fk_world(t), x)
end 