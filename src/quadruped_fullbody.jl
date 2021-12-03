abstract type AbstractQuadruped <: RobotDynamics.AbstractModel end 

struct UnitreeA1FullBody <: AbstractQuadruped
    rigidbody::RigidBodyModel 
    contacts::Vector{Point3D}

    function UnitreeA1FullBody(mech::Mechanism) 
        contacts = build_contact_points(mech)
        control_indices = Vector{Bool}([zeros(6);ones(12)])
        model = RigidBodyModel(mech, control_indices)
        new(model, contacts)
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

RobotDynamics.state_dim(model::AbstractQuadruped) = state_dim(model.rigidbody)
RobotDynamics.control_dim(model::AbstractQuadruped) = control_dim(model.rigidbody)


"""
    build_contact_points(model::UnitreeA1)
    Establish a list of contact points (Point3D). Returns them in a list  
    First four are RL, RR, FL, FR foot point 
    Last 8 are the eight corners that bounds the torso 
"""
function build_contact_points(mech::Mechanism)
    RL_foot = findbody(mech, "RL_foot")
    RR_foot = findbody(mech, "RR_foot")
    FL_foot = findbody(mech, "FL_foot")
    FR_foot = findbody(mech, "FR_foot")
    body = findbody(mech, "trunk")
    
    # Foot contact points 
    RL_point = Point3D(default_frame(RL_foot), 0,0,0)
    RR_point = Point3D(default_frame(RR_foot), 0,0,0)
    FL_point = Point3D(default_frame(FL_foot), 0,0,0)
    FR_point = Point3D(default_frame(FR_foot), 0,0,0)
    
    # Body contact points 
    w_b = 0.267/2 # width
    l_b = 0.194/2 # length
    h_b = 0.114/2   # height
    bp1 = Point3D(default_frame(body), -w_b,-l_b, -h_b)
    bp2 = Point3D(default_frame(body), -w_b,-l_b,  h_b)
    bp3 = Point3D(default_frame(body), -w_b, l_b, -h_b)
    bp4 = Point3D(default_frame(body), -w_b, l_b,  h_b)
    bp5 = Point3D(default_frame(body), w_b, -l_b, -h_b)
    bp6 = Point3D(default_frame(body), w_b, -l_b,  h_b)
    bp7 = Point3D(default_frame(body), w_b,  l_b, -h_b)
    bp8 = Point3D(default_frame(body), w_b,  l_b,  h_b)

    return [RL_point, RR_point, FL_point, FR_point, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8]
end



""" 
    fk_world(model::UnitreeA1, q::AbstractVector{T}, contacts)
    Given model, model configuration, an array of Point3D, returns the forward kinematics 
    of the list of contact points in the global frame.  

    fk = [x_RL, y_RL, z_RL, x_RR, y_RL, z_RL, ...]
"""
function fk_world(model::AbstractQuadruped, x::AbstractVector{T}) where T
    contacts = model.contacts
    world = findbody(model.rigidbody.mech, "world")
    fk = zeros(T, length(contacts)*3)
    state = model.rigidbody.statecache[T]
    # set_configuration!(state, q)
    copyto!(state, x)
    
    for i in 1:length(contacts)
        g = relative_transform(state, contacts[i].frame, default_frame(world)).mat
        R = @view g[1:3, 1:3]
        p = @view g[1:3,end]        
        fk[(i-1)*3+1:i*3] .= p + R * contacts[i].v
    end 
    
    return fk
end

""" 
    fk_world(model::UnitreeA1, q::AbstractVector{T}, contacts)
    Given model, model configuration, an array of Point3D, returns the forward kinematics 
    of the list of contact points in the body frame. 
    fk = [x_RL, y_RL, z_RL, x_RR, y_RL, z_RL, ...]
"""
function fk_body(model::AbstractQuadruped, x::AbstractVector{T}) where T 
    contacts = model.contacts
    body = findbody(model.rigidbody.mech, "base")
    fk = zeros(T, length(contacts)*3)
    state = model.rigidbody.statecache[T]
    copyto!(state,x)

    for i in 1:length(contacts)
        g = relative_transform(state, contacts[i].frame, default_frame(body)).mat
        p = @view g[1:3,end]
        R = @view g[1:3,end]
        fk[(i-1)*3+1:i*3] .= p
    end
    return fk
end 

""" 
    Kinematic jacobian of the contact points in body frame 
"""
function body_jacobian(model::AbstractQuadruped, x::AbstractVector{T}) where T
    return ForwardDiff.jacobian(t->(fk_body(model, t)), x)
end 

""" 
    Kinematic jacobian of the contact points in spatial frame 
"""
function spatial_jacobian(model::AbstractQuadruped, x::AbstractVector{T}) where T 
    return ForwardDiff.jacobian(t->(fk_world(model, t)), x)
end 


###### Experimental below ########
function contacts_body_jacobian(model::AbstractQuadruped, x::AbstractVector{T}) where T 
    # body jacobian 
    contacts = model.contacts
    imu_link = findbody(model.rigidbody.mech, "imu_link")
    J = zeros(T, length(contacts)*3, num_velocities(model.rigidbody.mech))
    state = model.rigidbody.statecache[T]
    # set_configuration!(state, q)
    copyto!(state, x)
    
    for i in 1:length(contacts)
        body = body_fixed_frame_to_body(model.rigidbody.mech, contacts[i].frame) # return body fixed to the frame 
        imu2body = path(model.rigidbody.mech, imu_link, body)
        g = relative_transform(state, contacts[i].frame, default_frame(imu_link)).mat
        # g = relative_transform(state, contacts[i].frame, default_frame(imu_link)).mat
        J_p = g[1:3,1:3] * point_jacobian(state, imu2body, contacts[i]).J # the original J is in body frame
        # the g[1:3, 1:3] rotation matrix rotate the velocity vector back to world frame 
        J[(i-1)*3+1:i*3, :] .= J_p
    end 
    return J
end 

""" 
    contacts_jacobian(model::UnitreeA1FullBody, q::AbstractVector{T}, contacts)
Given model, model configuration, an array of Point3D in contacts, return a stacked matrix 
of the joint Jacobian for each contact points.  
"""
# point_jacobian: https://juliarobotics.org/RigidBodyDynamics.jl/stable/algorithms/#RigidBodyDynamics.point_jacobian-Union{Tuple{C},%20Tuple{M},%20Tuple{X},%20Tuple{MechanismState{X,M,C,JointCollection}%20where%20JointCollection,RigidBodyDynamics.Graphs.TreePath{RigidBody{M},Joint{M,JT}%20where%20JT%3C:JointType{M}},Point3D}}%20where%20C%20where%20M%20where%20X
function contacts_jacobian(model::AbstractQuadruped, x::AbstractVector{T}) where T 
    # body jacobian 
    contacts = model.contacts
    world = findbody(model.rigidbody.mech, "world")
    imu_link = findbody(model.rigidbody.mech, "imu_link")
    J = zeros(T, length(contacts)*3, num_velocities(model.rigidbody.mech))
    state = model.rigidbody.statecache[T]
    # set_configuration!(state, q)
    copyto!(state, x)
    
    for i in 1:length(contacts)
        body = body_fixed_frame_to_body(model.rigidbody.mech, contacts[i].frame) # return body fixed to the frame 
        world2body = path(model.rigidbody.mech, world, body)
        g = relative_transform(state, contacts[i].frame, default_frame(world)).mat
        # g = relative_transform(state, contacts[i].frame, default_frame(imu_link)).mat
        J_p = g[1:3,1:3] * point_jacobian(state, world2body, contacts[i]).J # the original J is in body frame
        # the g[1:3, 1:3] rotation matrix rotate the velocity vector back to world frame 
        J[(i-1)*3+1:i*3, :] .= J_p
    end 
    return J
end 

function contacts_spatial_jacobian(model::AbstractQuadruped, x::AbstractVector{T}) where T 
    contacts = model.contacts
    world = findbody(model.rigidbody.mech, "world")
    imu_link = findbody(model.rigidbody.mech, "imu_link")
    J = zeros(T, length(contacts)*3, num_velocities(model.rigidbody.mech))
    state = model.rigidbody.statecache[T]
    # set_configuration!(state, q)
    copyto!(state, x)
    
    for i in 1:length(contacts)
        body = body_fixed_frame_to_body(model.rigidbody.mech, contacts[i].frame) # return body fixed to the frame 
        world2body = path(model.rigidbody.mech, world, body)
        g_body2world = relative_transform(state, contacts[i].frame, default_frame(world)).mat
        g_body2imu = relative_transform(state, contacts[i].frame, default_frame(imu_link)).mat

        J_p = point_jacobian(state, world2body, contacts[i]).J # the original J is in body frame
        J_p[:,1:3] = g_body2world[1:3,1:3] * J_p[:,1:3]
        J_p[:,4:6] = g_body2imu[1:3,1:3] * J_p[:,4:6] # this part confuses me the most. Why is this to body frame?
        J_p[:,7:end] = g_body2world[1:3,1:3] * J_p[:,7:end]
        J[(i-1)*3+1:i*3, :] .= J_p
    end 
    return J
end 