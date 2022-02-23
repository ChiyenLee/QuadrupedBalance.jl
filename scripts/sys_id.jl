using Pkg; Pkg.activate(".")
using Revise 
using QuadrupedBalance 
using RobotDynamics 
using LinearAlgebra 
using SparseArrays 
using PyPlot 
using ForwardDiff
using FFTW
using Rotations 
const QB = QuadrupedBalance
pygui(true)

## Initialize a model 
m = 12.454 # kg 
p_FR = [0.1308, -0.1308, 0.0]
p_FL = [0.1308,  0.1308, 0.0]
p_RR = [-0.1308, -0.1308, 0.0]
p_RL = [-0.1308,  0.1308, 0.0]
# I_inertia = Diagonal([0.1,0.25,0.3])
I_intertia = [[0.1, 0.0, 0.0], 
              [0.0, 0.25,0.],
              [0.0, 0., 0.3]]
model = QB.CentroidalModel(I_inertia, m, p_FR, p_FL, p_RR, p_RL)

## Come up with a simple LQR controller to control headings 
x_eq = zeros(13); x_eq[7] = 1.0; 
u_eq = zeros(12); u_eq[3:3:12] .= m * 9.81 / 4

h = 0.001 # discretization step
A = ForwardDiff.jacobian(t->dynamics(model,t, u_eq), x_eq)
B = ForwardDiff.jacobian(t->dynamics(model, x_eq, t), u_eq)
Ad = A*h + I # super coarse approximate discretizatoin step 
Bd = B*h 
G = blockdiag(sparse(1.0I(6)), 
              sparse(QB.quaternion_differential(x_eq[7:10])), 
              sparse(1.0I(3)))
Ad = G' * Ad * G 
Bd = G' * Bd

Q = Matrix(Diagonal([20, 20, 10000., 1.0, 1.0, 10., 1600, 1800, 1500., 10., 5, 1.])) *1e-2
R = Matrix(Diagonal(kron(ones(4), [0.5, 0.5, .5]))) * 1e-2 
P, K = QB.ricatti(Ad, Bd, Q, R)

## Simulattion variables and params  
dt = 0.001 
tf = 5.0 
times = 0:dt:tf 
xs = zeros(length(times), 13); xs[1,7] = 1.0
ω₀ = 2 * π * 1.0
ang_des_list = sin.(times * ω₀) * 0.05
x_des = zeros(13); x_des[7] = 1.0
δx = zeros(12)

# roll, pitch, yaw velocity data to be collected 
ω1s = zeros(length(times), 3) 
ω2s = zeros(length(times), 3) 
ω3s = zeros(length(times), 3) 
u1s = zeros(length(times), 12)
u2s = zeros(length(times), 12)
u3s = zeros(length(times), 12) 
t1s = copy(times)
t2s = copy(times)
t3s = copy(times)

## Oscillation for yaw 
x3s = zeros(length(times), 13); x3s[1,7] = 1.0
for i in 1:length(times)-1
    quat = x3s[i,7:10]
    ang_des = QB.ρ([0,0,ang_des_list[i]])
    quat_err = QB.L(quat) * ang_des
    ang_err = quat_err[2:end]

    δx[1:6] = x3s[i,1:6] - x_des[1:6]
    δx[7:9] = ang_err 
    δx[10:12] = x3s[i,11:13] - x_des[11:13]
    u = -K * δx 
    u3s[i,:] = u

    x3s[i+1, :] = QB.dynamics_rk4(model, x3s[i,:], u, dt)
end 
ω3s = x3s[:,13]

## Oscillation for pitch 
x2s = zeros(length(times), 13); x2s[1,7] = 1.0 
for i in 1:length(times)-1
    quat = x2s[i,7:10]
    ang_des = QB.ρ([0,ang_des_list[i],0.])
    quat_err = QB.L(quat) * ang_des
    ang_err = quat_err[2:end]

    δx[1:6] = x2s[i,1:6] - x_des[1:6]
    δx[7:9] = ang_err 
    δx[10:12] = x2s[i,11:13] - x_des[11:13]
    u = -K * δx 
    u2s[i,:] = u

    x2s[i+1, :] = QB.dynamics_rk4(model, x2s[i,:], u, dt)
end 
ω2s = x2s[:,12]

## Oscillation for roll 
x1s = zeros(length(times), 13); x1s[1,7] = 1.0 
for i in 1:length(times)-1
    quat = x1s[i,7:10]
    ang_des = QB.ρ([ang_des_list[i],0., 0.])
    quat_err = QB.L(quat) * ang_des
    ang_err = quat_err[2:end]

    δx[1:6] = x1s[i,1:6] - x_des[1:6]
    δx[7:9] = ang_err 
    δx[10:12] = x1s[i,11:13] - x_des[11:13]
    u = -K * δx 
    u1s[i,:] = u

    x1s[i+1, :] = QB.dynamics_rk4(model, x1s[i,:], u, dt)
end 
ω1s = x1s[:,11] #.+ randn(length(times),1) * 0.1

## Run FFTs to get the exact sinusoids for the omegas 
# assume we hit steady state after t = 2.0 
function fit_cosine(data, f_sample)
    N = size(data,1)

    fs = collect(((0: (N÷2)-1))/N .* f_sample) # get the list of frequency
    as = fft(data,1)[1:N÷2,:] # fourier coefficient 
    a_off1 = as[2:end, :]
    inds_max = argmax(abs.(a_off1), dims=1) # find the max 
    magnitudes = abs.(a_off1[inds_max]) * 2/N
    phases = angle.(a_off1[inds_max])
    f = fs[[i[1] + 1 for i in inds_max ]]
    offset = abs.(as[1,:]) * 2/N .* sign.(angle.(as[1,:]))
    return magnitudes, phases, f, offset 
end 

# Calculate the interpolated ωs. Take their derivative to get ω̇
t_cutoff = 1.0 
ind_cutoff = convert(Int64, t_cutoff ÷ dt)
a1, ϕ1, f1, off1 = fit_cosine(ω1s[ind_cutoff:end], 1/dt)
a2, ϕ2, f2, off2 = fit_cosine(ω2s[ind_cutoff:end], 1/dt)
a3, ϕ3, f3, off3 = fit_cosine(ω3s[ind_cutoff:end], 1/dt)
cos_wave(a, ϕ, f, off, t) = a .* cos.(2*π*f .* t .+ ϕ) .+ off'
sin_wave(a, ϕ, f, off, t) = a .* sin.(2*π*f .* t .+ ϕ) .+ off'

## Fit control inputs. However, it's not exactly sinusoidal. 
## need to filter out the higher frequency noise and keep the rest 
function grf_to_torques(model, x, u)
    quat = x[7:10]
    p = x[1:3]
    Q = UnitQuaternion(quat)
    r1 = (model.p1 - p)
    r2 = (model.p2 - p)
    r3 = (model.p3 - p)
    r4 = (model.p4 - p)
    M = [I(3) I(3) I(3) I(3);
         QB.hat(r1) QB.hat(r2) QB.hat(r3) QB.hat(r4)]
    u_out = M * u 

    τ_body = Q' * u_out[4:6]
    return τ_body
end 

τ1s = zeros(length(times), 3)
τ2s = zeros(length(times), 3)
τ3s = zeros(length(times), 3)
for i in 1:size(u1s,1)
    τ1s[i,:] = grf_to_torques(model, x1s[i,:], u1s[i,:])
    τ2s[i,:] = grf_to_torques(model, x2s[i,:], u2s[i,:])
    τ3s[i,:] = grf_to_torques(model, x3s[i,:], u3s[i,:])
end 

at1, ϕt1, ft1, offt1 = fit_cosine(τ1s[ind_cutoff:end, :], 1/dt)
at2, ϕt2, ft2, offt2 = fit_cosine(τ2s[ind_cutoff:end, :], 1/dt)
at3, ϕt3, ft3, offt3 = fit_cosine(τ3s[ind_cutoff:end, :], 1/dt)

τ1_interp = cos_wave(at1, ϕt1, ft1, offt1, t1s)
τ2_interp = cos_wave(at2, ϕt2, ft2, offt2, t2s)
τ3_interp = cos_wave(at3, ϕt3, ft3, offt3, t3s)
ω1_interp = cos_wave(a1, ϕ1, f1, off1, t1s)
ω2_interp = cos_wave(a2, ϕ2, f2, off2, t2s)
ω3_interp = cos_wave(a3, ϕ3, f3, off3, t3s)


## Experiments with noise 
noise_level = 0.1 
ω1_noise = ω1s + randn(size(ω1s)) * noise_level
ω2_noise = ω2s + randn(size(ω2s)) * noise_level 
ω3_noise = ω3s + randn(size(ω3s)) * noise_level 
τ1_noise = τ1s + randn(size(τ1s)) * noise_level 
τ2_noise = τ2s + randn(size(τ2s)) * noise_level 
τ3_noise = τ3s + randn(size(τ3s)) * noise_level 

at1, ϕt1, ft1, offt1 = fit_cosine(τ1_noise[ind_cutoff:end, :], 1/dt)
at2, ϕt2, ft2, offt2 = fit_cosine(τ2_noise[ind_cutoff:end, :], 1/dt)
at3, ϕt3, ft3, offt3 = fit_cosine(τ3_noise[ind_cutoff:end, :], 1/dt)
a1, ϕ1, f1, off1 = fit_cosine(ω1_noise[ind_cutoff:end], 1/dt)
a2, ϕ2, f2, off2 = fit_cosine(ω2_noise[ind_cutoff:end], 1/dt)
a3, ϕ3, f3, off3 = fit_cosine(ω3_noise[ind_cutoff:end], 1/dt)

times_fit = 0:dt:tf
τ1_interp = cos_wave(at1, ϕt1, ft1, offt1, times_fit)
τ2_interp = cos_wave(at2, ϕt2, ft2, offt2, times_fit)
τ3_interp = cos_wave(at3, ϕt3, ft3, offt3, times_fit)
ω1_interp = cos_wave(a1, ϕ1, f1, off1, times_fit)
ω2_interp = cos_wave(a2, ϕ2, f2, off2, times_fit)
ω3_interp = cos_wave(a3, ϕ3, f3, off3, times_fit)
ω1_dot = sin_wave(-a1 * 2*π*f1, ϕ1, f1, 0.0, times_fit)
ω2_dot = sin_wave(-a2 * 2*π*f2, ϕ2, f2, 0.0, times_fit)
ω3_dot = sin_wave(-a3 * 2*π*f3, ϕ3, f3, 0.0, times_fit)
b = []
H = []

for i in 1:length(times_fit)
    A1x = [ω1_interp[i] 0 0 0 0 0; 
           0 ω1_interp[i] 0 0 0 0;
           0 0 ω1_interp[i] 0 0 0]
    A1y = [ 0 ω2_interp[i] 0 0 0 0; 
            0 0 0 ω2_interp[i] 0 0;
            0 0 0 0 ω2_interp[i] 0]
    A1z = [ 0 0 ω3_interp[i] 0 0 0; 
            0 0 0 0 ω3_interp[i] 0;
            0 0 0 0 0 ω3_interp[i]]
    A1 = vcat([A1x, A1y, A1z])

    A2x = [ω1_dot[i] 0 0 0 0 0; 
            0 ω1_dot[i] 0 0 0 0;
            0 0 ω1_dot[i] 0 0 0]
    A2y = [ 0 ω2_dot[i] 0 0 0 0; 
             0 0 0 ω2_dot[i] 0 0;
             0 0 0 0 ω2_dot[i] 0]
    A2z = [ 0 0 ω3_dot[i] 0 0 0; 
             0 0 0 0 ω3_dot[i] 0;
             0 0 0 0 0 ω3_dot[i]]
    A2 = vcat([A2x, A2y, A2z])

    append!(H, A1 + A2)    
    append!(b,[τ1_interp[i,:]...,τ2_interp[i,:]..., τ3_interp[i,:]...] )

end 
H = vcat(H...)
I_est = pinv(H) * b

## Plotting the three exicitation routines
subplot(3,1,1)
plt.plot(t1s, ω1_noise, label="roll velocity")
plt.plot(t1s, ω1_interp, label="interpolated roll vel")
subplot(3,1,2)
plt.plot(t2s, ω2_noise, label="pitch velocity")
plt.plot(t2s, ω2_interp, label="interpolated roll vel")
subplot(3,1,3)
plt.plot(t3s, ω3_noise, label="yaw velocity")
plt.plot(t3s, ω3_interp, label="interpolated roll vel")

plt.legend() 

figure() 
subplot(3,1,1)
plot(t1s, τ1_noise)
plot(t1s, τ1_interp)
subplot(3,1,2)
plot(t2s, τ2_noise)
plot(t2s, τ2_interp)
subplot(3,1,3)
plot(t3s, τ3_noise)
plot(t3s, τ3_interp)
