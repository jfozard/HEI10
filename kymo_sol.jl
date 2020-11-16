


using ArgParse

using DifferentialEquations
using DiffEqCallbacks
using LinearAlgebra
using Random
using SparseArrays
using StatsBase
using Distributions


# Piecewise linear functions

# Piecewise linear interpolation at multiple positions (specialized for Float64)
function interp1d(x::Array{Float64,1}, x0::Array{Float64,1}, y0::Array{Float64,1})
    # x - position to evaluate function
    # x0 - positions of discontinuities in function
    # y0 - values of function on each region: y = y0[i] for x0[i] < x < x0[i+1]


    y = zeros(Float64, size(x))
    jj = 1
    for ii in 1:length(y)
        if x0[1]<=x[ii]<x0[end]
            jj = findnext(val->val>x[ii],x0,jj)-1
            y[ii] = y0[jj] + (y0[jj+1]-y0[jj])/(x0[jj+1]-x0[jj])*(x[ii]-x0[jj])
        end
    end
    return y
end

# Piecewise linear interpolation at a single position
function interp1d(x::Number, x0::Array, y0::Array)
    # x - position to evaluate function
    # x0 - positions of discontinuities in function slope
    # y0 - values of function at each discontinuity. Function takes the value 0 outside the range [x0[1], x0[length(x0)]]

    jj = 1
    y = 0.0
    if x0[1]<=x<x0[end]
        jj = findnext(val->val>x,x0,jj)-1
        y = y0[jj] + (y0[jj+1]-y0[jj])/(x0[jj+1]-x0[jj])*(x-x0[jj])
    end
    return y
end


# Piecewise constant interpolation 
function piecewise(x::Number, x0::Array, y0::Array)
    # x - position to evaluate function
    # x0 - positions of discontinuities in function
    # y0 - values of function on each region: y = y0[i] for x0[i] < x < x0[i+1]

    jj = 1
    y = 0.0
    if x0[1]<=x<x0[end]
        jj = findnext(val->val>x,x0,jj)-1
        y = y0[jj]
    end
    return y
end


# Inverse of a piecewise constant density function
function piecewise_density_inv(x0::Array{T,1}, y0::Array{T,1}) where {T <: Number}
    csum = zeros(T, size(x0))
    for i in 2:length(x0)
        csum[i] = (x0[i] - x0[i-1])*y0[i-1] + csum[i-1]
    end
    csum = csum / csum[end]
    csum, x0
end


# Inverse of a piecewise linear density function
function linear_density_inverse(x0::Array{T,1}, y0::Array{T,1}) where {T <: Number}
    csum = zeros(T, size(x0))
    for	i in 2:length(x0)
    	csum[i] = (x0[i] - x0[i-1])*(y0[i-1] + y0[i])/2	+ csum[i-1]
    end
    csum_t = csum[end]
    function inverse(x::T)
        z = x*csum_t
        jj = findnext(val->val>z, csum, 1)-1
        z = z - csum[jj]
        y = y0[jj]
        dy = (y0[jj+1] - y0[jj])/(x0[jj+1]-x0[jj])
        if abs(dy)>1e-12
            return x0[jj]+((sqrt(y*y+2*z*dy)-y)/dy)
        elseif abs(y)>1e-12
            return x0[jj]+z/y
        else
            return x0[jj]
        end
    end
    inverse
end

# L^2 / sum-squared vector norm
function norm(x)
    sqrt(sum(abs2,x))
end

# RI HEI10 escape rate function (\beta(C)*C in equation (2) )
function betaC(C, K, n_c)
    C ./ (1 .+ (C / K) .^ n_c)
end

# RI HEI10 escape rate function \beta(C)*C (evaluated in place)
function betaC!(dC, C, K, n_c)
    @. dC = C ./ (1 .+ (C ./ K) .^ n_c)
end

# Derivative of (\beta(C)*C) wrt C 
function d_betaC(C::Float64, K::Float64, n_c::Float64)::Float64
    1 ./ (1.0 .+ (C ./ K) .^ n_c) .-
    C ./ (1.0 .+ (C ./ K) .^ n_c) .^ 2 .* (n_c ./ K) .* (C ./ K) .^ (n_c - 1)
end

# Derivative of (\beta(C)*C) wrt C (evaluated in place)
function d_betaC!(dbC, C, K, n_c)
    @. dbC = 1 ./ (1.0 .+ (C ./ K) .^ n_c) .-
    C ./ (1.0 .+ (C ./ K) .^ n_c) .^ 2 .* (n_c ./ K) .* (C ./ K) .^ (n_c - 1)
end



# Function to calculate the RHS of the ODE system
# Result returned in-place as first argument
function basic_ode!(dyy::Array{Float64,1}, yy::Array{Float64,1}, p, t::Float64)
    m,           # Number of SC spatial compartments for discretization
    N,           # Number of RIs
    K,           # Hill threshold
    nodes,       # Indices of SC spatial compartments containing RIs
    alpha,       # Rate of HEI10 production on SC - here set to be zero
    beta,        # Rate of HEI10 escape from SC - here set to be zero
    D,           # RI Diffusion coefficient on SC (D_c)
    dL,          # Length of each spatial compartment
    betaC!,      # beta(C)*C function (passed to ODE for efficiency)
    d_betaC,     # derivative of beta(C)*C function wrt C
    n_c,         # Hill coefficient (\gamma)
    a_nodes,     # RI HEI10 absorption rate from SC (\alpha)
    b_nodes,     # RI HEI10 escape rate back onto SC (\beta)
    b2_nodes,    # RI HEI10 escape rate into nucleoplasm  (not used)
    b3_nodes,    # RI HEI10 rate  (not used)
    bC,          # Placeholder to re-use vector
    A = p

    # Extract variables from ODE State vector
    @inbounds y = @view yy[1:m] # HEI10 concentration on the SC (c_i, 1<=i<=m)
    @inbounds dy = @view dyy[1:m] 
    @inbounds C = @view yy[m+1:m+N] # HEI10 amounts at the RH (C_j 1<=j<=N)
    @inbounds dC = @view dyy[m+1:m+N]

    h = D/dL/dL

    # Diffusion (production and decay - not used) of HEI10 on the SC.
    # Standard second-order finite difference approximation of 1D laplacian
    @inbounds dyy[1] = alpha - beta * y[1] + h * (y[2] - y[1])
    @inbounds @simd for i = 2:(m-1)
        dyy[i] = alpha - beta * y[i] +
                h * (y[i+1] - 2 * y[i] + y[i-1])
    end

    @inbounds dyy[m] = alpha - beta * y[m] + h * (y[m-1] - y[m])

    # Calculate escape of HEI10 back onto SC.
    betaC!(bC, max.(C,0.0), K, n_c)

    # Calculate net fluxes from SC to RIs, and update ODE RHS for appropriate SC compartment
    @inbounds for j = 1:N
        A[j] = a_nodes * y[nodes[j]] - b_nodes * bC[j]
        dyy[nodes[j]] = dyy[nodes[j]] - A[j] / dL
    end

    # Rate of change of HEI10 amounts at the RIs
    @. dC = A .- b2_nodes .* bC .- b3_nodes .* C
end


# Function to calculate the Jacobian of the ODE system RHS - helps implicit integration method be more efficient 
# Result returned in-place as first argument
function jac!(J, yy, p, t)
    m, # For parameter interpretation please see basic ODE
    N,
    K,
    nodes,
    alpha,
    beta,
    D,
    dL,
    betaC!,
    d_betaC,
    n_c,
    a_nodes,
    b_nodes,
    b2_nodes,
    b3_nodes,
    bC,
    A = p
    r = D / dL / dL
    C = @view yy[m+1:m+N]
    #@show typeof(J)
    J[1, 1] = -beta - r
    J[1, 2] = r
    @inbounds for j = 2:m-1
        J[j, j-1] = r
        J[j, j] = -beta - 2 * r
        J[j, j+1] = r
    end
    @inbounds J[m, m-1] = r
    @inbounds J[m, m] = -beta - r
    dbC::Float64 = 0.0
    idx::Int64 = 0
    @inbounds for j = 1:N
        dbC = d_betaC(max(0, C[j]), K, n_c)
        J[m+j, m+j] = -(b_nodes + b2_nodes) * dbC - b3_nodes
        J[nodes[j], m+j] = b_nodes * dbC / dL
        J[m+j, nodes[j]] = a_nodes
        J[nodes[j], nodes[j]] -= a_nodes / dL
    end
end



function simulate_video(idx, args)

    mt = MersenneTwister(idx)


    L::Float64 = Float64(args["L"])+Float64(args["Ls"])*clamp(randn(mt), -3, 3)
    N::Int64 = Int64(round(L*args["density"]))
    m::Int64 = Int64(args["m"])

    x::Array{Float64,1} = L * sort(rand(mt, N))

    dL::Float64 = L / m
    nodes = map(a -> ceil(Int, a), x / dL)


    xm = (0.5:1.0:m) * dL

    n_c::Float64 = Float64(args["n_c"])
    K::Float64 = Float64(args["K"])
    alpha = 0.0
    beta = 0.0
    a_nodes::Float64 = Float64(args["a_nodes"])
    b_nodes::Float64 = Float64(args["b_nodes"])
    b2_nodes = 0.0
    b3_nodes = 0.0
    D::Float64 = Float64(args["D"])

    bC = zeros(Float64, N)
    A = zeros(Float64, N)

    u0::Array{Float64,1} = fill(Float64(args["u0"]), m)
    C0_t = 1.0 
    C0_t2::Float64 = Float64(args["t_C0_ratio2"])


    d = Truncated(Normal(Float64(args["C0"]), Float64(args["C0_noise"])), max(0, Float64(args["C0"])-3*Float64(args["C0_noise"])), Float64(args["C0"])+3*Float64(args["C0_noise"]))


    C0::Array{Float64,1} = rand(mt, d, N)
    C0[C0.<0] .= 0.0

    t_L::Float64 = Float64(args["t_L"])


    # Additional HEI10 loading for RI near ends of SCs (telomeres)
    if !args["t_exp"]
        # Piecewise linear function ($f$)
        for i in 1:N
    	    p = x[i]/L
            if p < t_L
               r = p/t_L
            C0[i] *= C0_t*(r) + C0_t2*(1-r)
            elseif p> 1-t_L
               r = (1-p)/t_L
               C0[i] *= C0_t*(r) + C0_t2*(1-r)
       	     end
	 end
    else
        # Exponentially decaying RI end weights
	for i in 1:N
    	    p = x[i]/L
	    C0[i] *= (C0_t2-C0_t)*(exp(-2*p/t_L) + exp(-2*(1-p)/t_L))+ C0_t
	 end
    end

    # Construct initial ODE state vector
    y0::Array{Float64,1} = vcat(u0, C0)

    # ODE system parameters
    p = (
        m,
        N,
        K,
        nodes,
        alpha,
        beta,
        D,
        dL,
        betaC!,
        d_betaC,
	n_c,
        a_nodes,
        b_nodes,
        b2_nodes,
        b3_nodes,
        bC,
        A,
    )

    # Initialize ODE problem with jacobian
    j = spzeros(m + N, m + N)
    g2 = ODEFunction(basic_ode!; jac=jac!, jac_prototype = j)
    prob_g = ODEProblem(g2, y0, (0.0, Float64(args["n_ts"]) * Float64(args["dt"])), p)
    cb = PositiveDomain(g2) # Additional constraint to ensure non-negativity of concentration values during integration

    # Solve ODE system using implicit Rodas5 solver, with tight relative error tolerance
    sol = solve(prob_g,Rodas5(autodiff=:false), cb=cb, reltol=1e-12)
    # Return simulation results at all timepoints
    (sol, x, xm, L, N)
end



function parse_commandline()

# Parse command line options and parameters
# Note that the driver file (run_simulations.py) overrides most of these parameter values

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--filename"
            default="kymo.png"
            help="base filename"
	"--timecourse"
        action = :store_true
	    help="timecourse instead"
	"--start"
           arg_type=Int
           default=1
           help="first seed"
        "--n_ts"
            arg_type=Int
            default=100
            help="number of timesteps"
        "--L"
	    arg_type=Float64
	    default=45.3 # (um)
            help="length of chromosome" # Mean length of this SC
        "--Ls"
            arg_type=Float64
	    default=5.52 # (um)
            help="std dev in length of chromosome"  # Variation in length of this SC
        "--dt"
     	    arg_type=Float64
	    default=360.0 # (s)
            help="time interval" # Length of each integration period (timestep is shorter, controlled automatically by error estimate)
        "--m"
	    arg_type=Int
	    default=2000 # (no unit)
            help="number of elements" # (M) Number of intervals for spatial discretization of SC 
        "--K"
     	    arg_type=Float64
            default=1.0 # (a.u.)
            help="threshold for RI unbinding" # (K_C) Hill function thresold for escape of HEI10 from RI 
        "--n_c"
    	    arg_type=Float64
    	    default=1.25 # (no units)
            help="Hill coefficient for RI unbinding" # (\gamma) Hill coefficient for escape of HEI10 from RI
        "--a_nodes"
            arg_type=Float64
            default=2.1 # (um s^{-1})
            help="binding rate at RIs" # (\alpha) RI HEI10 absorption rate from SC
        "--b_nodes"
            arg_type=Float64 
            default=0.5 # (s^{-1})
            help="unbinding rate at RIs" # (\beta) RI HEI10 escape rate back onto SC
        "--u0"
            arg_type=Float64
            default=1.2 # (a.u.)
            help="initial HEI10 concentration" # (c_0) HEI10 initial loading on SC
        "--C0"
            arg_type=Float64
            default=6.8 # (a.u.)
            help="initial RI HEI10" # (C_0) HEI10 initial RI loading
        "--C0_noise"
            arg_type=Float64
            default=2.2 # *(a.u.)
            help="initial RI HEI10 noise" # (\sigma) Noise in HEI10 initial RI loading 
        "--D"
            arg_type=Float64
            default=1.1 # (um^2 s^{-1}) 
            help="HEI10 diffusion constant on SC"  # (D_c) RI Diffusion coefficient on SC 
        "--t_C0_ratio2"
            arg_type=Float64
            default=2.0 # (no units)
	    help="initial hei10 ratio for RIs in telomeres"  # (f_e) Increased RI loading at SC telomeres
	"--t_exp"
	    help="exponential functions for telomeres"
	    action = :store_true
        "--t_L"
            arg_type=Float64
            default=0.1 # (no units)
            help="telomere rel length" # (x_e) Telomere length as fraction of bivalent 
        "--density"
            arg_type=Float64
            default=0.5         # (um^{-1})
            help="RI density"  # (\rho) SC RI density
        end
    return parse_args(s)
end


function kymo_sol()
    n_args = parse_commandline()
    println("Parsed args:")

    for (arg,val) in n_args
        println("$arg=$val")
    end

    sol, x, xm, L, N = simulate_video(n_args["start"], n_args)

    t = range(0, stop=sol.t[end], length=100)
    u0 = hcat(sol.u...)
    u = hcat(sol(t)...)

    @show(u[2001:end, end])
    @show(sum(u[2001:end, end].>1))


    open(n_args["filename"],"w") do io
       print(io, string(L) * ", " *string(sol.t[end]) *'\n'*join(map(string, x), ',')*'\n')
       for i=1:100
           v=sol(t[i])[2001:end]
           print(io, join(map(string, v), ',')*'\n') 
       end
    end

end

kymo_sol()
