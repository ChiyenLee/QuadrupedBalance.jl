
""" 
    mechanism: Dojo.mechanism
    joint_states: 
""" 
function set_state(mechanism, joint_states::Vector{Float64}, tran::Vector{Float64}, rot::Vector{Float64})
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :FR_hip_joint), [joint_states[1]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :FL_hip_joint), [joint_states[2]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :RR_hip_joint), [joint_states[3]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :RL_hip_joint), [joint_states[4]])

    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :FR_thigh_joint), [joint_states[5]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :FL_thigh_joint), [joint_states[6]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :RR_thigh_joint), [joint_states[7]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :RL_thigh_joint), [joint_states[8]])

    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :FR_calf_joint), [joint_states[9]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :FL_calf_joint), [joint_states[10]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :RR_calf_joint), [joint_states[11]])
    Dojo.set_position!(mechanism, get_joint_constraint(mechanism, :RL_calf_joint), [joint_states[12]])

    set_position!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [tran; rot])

    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [zeros(3); zeros(3)])

    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :FR_hip_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :FR_thigh_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :FR_calf_joint), [0.])

    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :FL_hip_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :FL_thigh_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :FL_calf_joint), [0.])

    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :RR_hip_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :RR_thigh_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :RR_calf_joint), [0.])

    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :RL_hip_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :RL_thigh_joint), [0.])
    Dojo.set_velocity!(mechanism, get_joint_constraint(mechanism, :RL_calf_joint), [0.])

end 


"""
    Remap from Dojo joint indices to QB indices 
"""
function dojo_to_qb(q1)
    inds_dojo = @SVector [1,4,7,10,2,5,8,11,3,6,9,12]
    q2 = q1[inds_dojo]
    return q2
end 

"""
    FR_hip, FR_thigh, FR_calf, FL_hip, FL_thigh, FL_calf 
"""
function qb_to_dojo(q1)
    inds_qb = @SVector [1,5,9,2,6,10,3,7,11,4,8,12]
    q2 = q1[inds_qb]
    return q2
end 