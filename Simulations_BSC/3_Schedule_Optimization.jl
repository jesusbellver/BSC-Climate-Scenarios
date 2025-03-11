import Pkg

#Packages for cluster
#Pkg.add("DifferentialEquations")
#Pkg.add("DataFrames")
#Pkg.add("CSV")
#Pkg.add("Interpolations")
#Pkg.add("Dates")
#Pkg.add("JLD2")
#Pkg.add("ODE")
#Pkg.add("LSODA")
#Pkg.add("XLSX")
#Pkg.add("LatinHypercubeSampling")
#Pkg.add("Random")
#Pkg.add("DelimitedFiles")
#Pkg.add("QuadGK")
#Pkg.add("StatsBase")


using DifferentialEquations
using DataFrames
using CSV
using Interpolations
using Dates
using JLD2
using ODE
using LSODA
using XLSX
using Base.Threads
using Random
using DelimitedFiles
using QuadGK
using StatsBase
using LinearAlgebra

#println("Packages Installed")
#flush(stdout)


#PARA PODER LANZARLO EN EL CLUSTER CON year COMO ARGUMENTO. COMENTAR SI NO!!!
# args = ARGS
# const year = parse(Int, args[1])
# const BETA_SIM = args[2]
# const idx_model = args[3]

#codigo .sh
#!/bin/bash

# export JULIA_NUM_THREADS=10

# read -p "Enter the simulation year: " year
# read -p "Enter the value for BETA_SIM: " beta_sim
# read -p "Enter the Model: " model

# /home/usuaris/jbellver/.juliaup/bin/julia BCN-Model/Simulations_BCN/6_Schedule_Optimization.jl "$year" "$beta_sim" "$model" > output_OptSch${year}_${beta_sim}_${model}.txt


#FUERA DEL CLUSTER
const year = 2019
const BETA_SIM = 0.1
const idx_model = 1


##########
#Specify model
list_models = ["CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "MIROC6", "MPI-ESM1-2-HR", "NorESM2-MM"]
model = list_models[idx_model]


println("Year: ", year, " Beta: ", BETA_SIM, " Model: ", model)
flush(stdout)

# Define the parameters (From Param_fit_R0)
const alpha = 1.1
const w_min = 3.5
const w_flush = 0.2
const K1 = 4
const K2 = 0.15


const iters_annealing = 8000
#C_T = 0.99985 #Para iters=20000
#C_T = 0.9997 #Para iters=10000
#C_T = 0.9998 #Para iters=20000


# C_T = 0.9996 #Para iters=10000
# Init_T = 0.072


#betas=1,0.8,0.6,0.4,0.2,0
# const C_T = 0.9996 #Para iters=10000
# const Init_T = 0.05
const C_T = 0.99955 #Para iters=10000
const Init_T = 0.0425
const Max_zones_day=3
const Max_days_week=4


#Functions to generate M, c and W
#Compute M
function generate_M(alpha)
    M=zeros(n,n)
    for i in 1:n
        for j in 1:n
            if i==j
                M[i, j] = 0  #Influence on themselves
            else
                M[i,j]=exp(-alpha*D[i,j])  #Influence on the rest of the zone
            end
        end
    end

    return M
end


#Compute Flushings
function generate_flushings(w_flush)
    Flushings=Dict{Int, Vector{Float64}}()
    for i in 1:n
        for j in 1:length(t_interp)
            if Rainfall(X_emb[i],Y_emb[i],t_interp[j]) > Prof_Sorrers[i] * w_flush
                if haskey(Flushings, i)
                    push!(Flushings[i], t_interp[j])
                else
                    Flushings[i] = [t_interp[j]]
                end
            end
        end
    end

    return Flushings
end

#Compute water dynamics without flushing
function generate_W_no_flush(w_min)
    W_no_flush=zeros(n,length(t_interp))
    for i in 1:n
        for j in 1:length(t_interp)
            if t_interp[j] <= t_reformas[i] && Rainfall_21(X_emb[i],Y_emb[i],t_interp[j]) >= w_min # Suitability of bs after they become directes = 0.
                W_no_flush[i, j] = 1            
            end
        end
    end

    return W_no_flush
end

#Compute c_W (matrix incorporating control, water dynamics and flushings)
function generate_W()

    # Initialize W
    W0 = copy(W_no_flush)
    for i in 1:n
        if haskey(Flushings, i)
            for flush_time in Flushings[i]
                # Find the closest index in t_interp that is greater than or equal to flush_time
                flush_index = findfirst(t -> t >= flush_time, t_interp)
                if flush_index !== nothing
                    W0[i, flush_index] = 0
                    W0[i, flush_index+1] = 0
                end
            end
        end
    end

    #Interpolate
    W = []
    # Perform interpolation for each coordinate separately
    for i in 1:n
        push!(W, LinearInterpolation(t_interp, W0[i,:], extrapolation_bc=Flat()))
    end

    return W
end

function generate_c(Interventions)

    #Compute control without flushing
    c_no_flush=zeros(n,length(t_interp))
    for i in 1:n
        if haskey(Interventions, i)
            for intervention_time in Interventions[i]
                if W[i](intervention_time) >= 0.0001 # No tratamos si no hay agua!!!!
                    t_index = findfirst(x -> x >= intervention_time, t_interp)
                    for j in t_index:length(t_interp)
                        if t_interp[j] <= intervention_time + 49
                            c_no_flush[i,j]=1
                        end
                    end
                end
            end
        end
    end

    # Add flushing to control
    c0 = copy(c_no_flush)
    for i in 1:n
        if haskey(Flushings, i)
            for flush_time in Flushings[i]
                # Find the closest index in t_interp that is greater than or equal to flush_time
                flush_index = findfirst(t -> t >= flush_time, t_interp)
                if flush_index !== nothing

                    # Find the intervention times greater than flush_time
                    if haskey(Interventions, i)
                        intervention_times = Interventions[i]
                        upcoming_intervention = findfirst(t -> t > flush_time, intervention_times)

                        if upcoming_intervention !== nothing
                            intervention_time = intervention_times[upcoming_intervention]
                            for j in flush_index:length(t_interp)
                                if t_interp[j] >= flush_time && t_interp[j] < intervention_time
                                    c0[i, j] = 0
                                end
                            end
                        else
                            for j in flush_index:length(t_interp)
                                if t_interp[j] >= flush_time
                                    c0[i, j] = 0
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    #Interpolate
    c = []
    # Perform interpolation for each coordinate separately
    for i in 1:n
        push!(c, LinearInterpolation(t_interp, c0[i,:], extrapolation_bc=Flat()))
    end

    return c
end



#LOAD DATA FOR SIMULATIONS
#Load Rainfall Data
const Rainfall=load_object(joinpath(@__DIR__, "Parameters/Rainfall_$(year)_$(model)"))
const Rainfall_21=load_object(joinpath(@__DIR__, "Parameters/Rainfall_21_$(year)_$(model)"))

#Gestionar reformas de embornales
const Unique_BS = CSV.read(joinpath(@__DIR__,"Dades_netes/Unique_Items.csv"), DataFrame) 

#Load Zones data
const N=length(unique(Unique_BS.nom_zr)) #number of risk zones
const Zone_Names = unique(Unique_BS.nom_zr)
const Zone_Indices = Dict(zone => findall(row -> row.nom_zr == zone, eachrow(Unique_BS)) for zone in Zone_Names)


#Load embornal data
D = load_object(joinpath(@__DIR__, "Parameters/D_Matrix_$(year)"))
const n = length(D[1,:])
const X_emb = collect(values(load_object(joinpath(@__DIR__, "Parameters/X_coords"))))
const Y_emb = collect(values(load_object(joinpath(@__DIR__, "Parameters/Y_coords"))))
const Prof_Sorrers = collect(values(load_object(joinpath(@__DIR__, "Parameters/Profunditat_Sorrers"))))


#Load R0
const R0 = load_object(joinpath(@__DIR__, "Parameters/R0_culex_$(year)_$(model)"))



#=         #Select one zone
zone_name="Jardins de Vil·la Amèlia"
Revisions=Revisions[Revisions.nom_zr.==zone_name,:]
n=length(unique(Revisions.id_item))
D=D[Unique_BS.nom_zr.==zone_name,Unique_BS.nom_zr.==zone_name]
c=c[Unique_BS.nom_zr.==zone_name]
X_emb=X_emb[Unique_BS.nom_zr.==zone_name]
Y_emb=Y_emb[Unique_BS.nom_zr.==zone_name]
Prof_Sorrers=Prof_Sorrers[Unique_BS.nom_zr.==zone_name]
Unique_BS=Unique_BS[Unique_BS.nom_zr.==zone_name,:] =#

#FUNCTIONS FOR THE SIMULATIONS

function generate_schedule(init_times, periodicity)
    Schedule = Dict{String, Vector{Float64}}()

    for i in 1:N
        Schedule[Zone_Names[i]]=[init_times[i] + k * periodicity[i] for k in 0:(div(end_mosq_season - init_times[i], periodicity[i]))]
    end

    return Schedule
end

function check_schedule(Schedule)

    #Check Max zones per day is satisfied
    check_times =sort(vcat(collect(values(Schedule))...))
    day_counts = vcat(collect(values(countmap(check_times)))...)
  
    #check max days per week is satisfied
    max_count = 0
    count = 0
    for i in 1:length(check_times)
        for j in (i+1):length(check_times)
            if check_times[j] - check_times[i] < 7
                count += 1
            else
                max_count = max(max_count,count) #store the maximum amount of days in a week so far
                count = 0 # Restart the counter and stop checking if the difference is already 7 or more
                break
            end
        end
    end

    if any(x -> x > Max_zones_day, day_counts) == false && max_count <= Max_days_week
        return true, Schedule
    else
        return false
    end
end

function random_Schedule()
    while true
        #generate times and periodicities
        init_times = round.(init_mosq_season .+ 5 .+ 51 .* rand(N)) .+ 0.5 #Los tratamientos van siempre a las 12:00 del mediodia. Generar trats iniciales entre 1 Abril y 31 Mayo pero evitando extremos.
        periodicity = round.(35 .+ (end_mosq_season-init_mosq_season + 1 -35) .* rand(N)) #los periodos tienen que ser días enteros. Generar periodios iniciales entre 25 y end_mosq_season-init_mosq_season -25 días (límites 1 y (end_mosq_season-init_mosq_season)+1)
        Schedule = generate_schedule(init_times,periodicity)

        #Check Interventions respect constraints
        Schedule = check_schedule(Schedule) #Schedule[1]=True,False, si True Interventions[2] contiene intervenciones

        if Schedule[1] == true
            return init_times, periodicity, Schedule[2]
        end
    end
end


function update_Schedule(init_times_0, periodicity_0, eps_t, eps_p)

    while true
        #alter times and periodicities
        init_times_1 = clamp.(init_times_0 .+ eps_t .* (2 .* round.(rand(N)) .- 1), init_mosq_season .+ 0.5, init_mosq_season .+ 61 .+ 0.5) #2*rand-1 produces a rand between -1 and 1. clamp entre 1 Abril y 31 Mayo
        periodicity_1 = clamp.(periodicity_0 .+ eps_p .* (2 .* round.(rand(N)) .- 1), 1, end_mosq_season-init_mosq_season + 1)
        Schedule = generate_schedule(init_times_1,periodicity_1)

        #Check Interventions respect constraints
        Schedule = check_schedule(Schedule) #Interventions[1]=True,False, si True Interventions[2] contiene intervenciones

        if Schedule[1] == true
            return init_times_1, periodicity_1, Schedule[2]
        end
    end
end

function generate_Interventions(Schedule)
    Interventions = Dict{Float64, Vector{Float64}}()
    for i in Zone_Names
        for j in Zone_Indices[i]
            Interventions[j] = Schedule[i]
        end
    end
    return Interventions
end

#Functions for integration
function fun!(du, u, p, t)
    tmp = similar(u)
    mul!(tmp, M, u)  # In-place multiplication: tmp = M * u
    @inbounds for i in eachindex(u)
        # Combine element-wise operations in a single loop to avoid temporary arrays.
        du[i] = tmp[i] * (1.0 - u[i]) * R0(X_emb[i], Y_emb[i], t) * K1 - K2 * u[i]
    end
end

function condition(u, t, integrator)
    #Define condition for discontinuity
    c = integrator.p[2]

    for i in 1:n
        if W[i](t) < 0.0001
            return true
        end
    end

    for i in 1:n
        if c[i](t) > 0.9999
            return true
        end
    end

    return false
end

function affect!(integrator)
    # Modify the solution when discontinuity occurs
# Inicializamos las condiciones como arreglos booleanos
    condition_1 = Vector{Bool}(undef, n)
    for i in 1:n
        condition_1[i] = W[i](integrator.t) < 0.0001
    end
    condition_2 = fill(false, n)
    Interventions = integrator.p[1]
    c = integrator.p[2]


    # Iteramos sobre cada zona y sus índices
    for zone in Zone_Names
        for i in Zone_Indices[zone]
            if haskey(Interventions, i)
                # Buscar la última intervención que se hizo <= integrator.t
                last_intervention_index = nothing
                for idx in reverse(eachindex(Interventions[i]))
                    if Interventions[i][idx] <= integrator.t
                        last_intervention_index = idx
                        break 
                    end
                end

                if last_intervention_index !== nothing #Si se hizo intervención
                    #Ver si hay larvas en la zona (>1)
                    s = 0.0
                    t_last = Interventions[i][last_intervention_index]
                    sol_last = integrator.sol(t_last)
                    for j in Zone_Indices[zone]
                        s += sol_last[j]
                    end

                    if s >= 1 #Si había larvas mantener el tratamiento hasta que c<1, si no eliminar tratamiento
                        condition_2[i] = c[i](integrator.t) > 0.9999
                    end
                end
            end
        end
    end


    for i in 1:n
        if condition_1[i] || condition_2[i]
            integrator.u[i] = 0.0001
        end
    end

    return integrator.u
end

# Función para integración trapezoidal
function trapz(t, y)
    s = 0.0
    for i in 1:length(t)-1
        s += (t[i+1] - t[i]) * (y[i] + y[i+1]) / 2
    end
    return s
end


#Cost function
function Cost_function(sol, periodicity, beta, norm)
    t_points = sol.t
    values = [ sum(sol(t)) for t in t_points ]
    INT= trapz(t_points,values) / norm
    PER=sum(periodicity) / (7*(end_mosq_season-init_mosq_season+1)) #Max possible is end_mosq_season-init_mosq_season days * 8 zones
    println("Int: ", INT)
    println("Per: ", 1-PER)
    flush(stdout)
    return beta * INT + (1-beta) * (1-PER), INT, 1-PER
end

#Set simulation Parameters
# alg=RK4()
# atol=1e-17
# max_iters=1e6

const alg = Tsit5()
const tol = 1e-14
const max_iters = 1e8

#SET SIMULATION PARAMETERS
const d0 = Dates.DateTime("$(year)-03-01", "yyyy-mm-dd")
const t0 = 1
const df = Dates.DateTime("$(year)-12-31", "yyyy-mm-dd") #último día
const tf = 1 + Dates.value(df - d0)/(24*60*60*1000) + 1 #El último + 1 es porque simulamos hasta el final del día
const tspan = (t0, tf)
const t_vect=t0:tf
const t_interp = t_vect[1:(end-1)] #points where there's data for interpolation

const dims = Dates.DateTime("$(year)-04-01", "yyyy-mm-dd") # 1 de Abril 
const dems = Dates.DateTime("$(year)-11-30", "yyyy-mm-dd") #30 de noviembre
const init_mosq_season = 1 + Dates.value(dims - d0)/(24*60*60*1000)
const end_mosq_season = 1 + Dates.value(dems - d0)/(24*60*60*1000)
const t_mosq_season = init_mosq_season:end_mosq_season

const u0 = [0.0001 for i in 1:n]  # Initial conditions

const M = generate_M(alpha)
const W_no_flush = generate_W_no_flush(w_min)
const Flushings = generate_flushings(w_flush)
const W = generate_W()

#Simulated Annealing
function Simulated_Annealing(beta, beta_str, Temp, init_times_0, periodicity_0, J0, c, prob, p, cb, Results, norm, h_tag)
    for i in 1:iters_annealing

        println("Beta ", beta, " Iteration ", i)
        flush(stdout)

        Temp = Temp * C_T

        eps_t = 1 + round(4 * rand() * (1 - i / iters_annealing))  # genera entero entre 1 y 5, sesgado hacia 1 conforme avanza la simulación
        eps_p = 1 + round(2 * rand() * (1 - i / iters_annealing))  # genera entero entre 1 y 3, sesgado hacia 1 conforme avanza la simulación
        init_times_1, periodicity_1, Schedule_1 = update_Schedule(init_times_0, periodicity_0, eps_t, eps_p)
        Interventions_1 = generate_Interventions(Schedule_1)
        c = generate_c(Interventions_1)

        # Bandera para indicar si se acepta la actualización o no
        valid_update = true

        try
            p = (Interventions_1, c)
            prob = ODEProblem(fun!, u0, tspan, p)
            cb = DiscreteCallback(condition, affect!, save_positions=(true, true))  # Define el callback
           
           @time begin
                sol = solve(prob, alg, callback=cb, abstol=tol, maxiters=max_iters)  # Soluciona el problema ODE con el callback
                J1 = Cost_function(sol, periodicity_1, beta, norm)
            end

            # Verificar si se obtuvo un resultado problemático
            if J1[2] < 0.01 || sol.retcode != :Success
                println("Inestabilidad numérica")
                valid_update = false  # Se descarta la actualización; se mantienen los valores _0
            else
                # Evaluar la diferencia y aplicar el criterio de aceptación
                E = J1[1] - J0[1]
                if E < 0 || rand() < exp(-E / Temp)
                    J0 = J1
                    init_times_0 = init_times_1
                    periodicity_0 = periodicity_1
                end
            end

        catch e
            println("ERROR/WARNING: ", e)
            flush(stdout)
            valid_update = false
        end

        Results[i] = [J0, init_times_0, periodicity_0, Temp, valid_update]

        if i % 100 == 0
            println("beta=", beta, " iteration=", i, " invalid updates: ", count(x -> x[5] == false, values(Results)))
            save_object(joinpath(@__DIR__, "Results/Simulated_Annealing_$(year)_b$(beta_str)_$(model)_$(h_tag)"), Results)
            if count(x -> x[5] == false, values(Results)) >= 1000
                break
            end
        end
    end

    return Results
end



#SIMULATION
function Run_Sim(BETA_SIM)
    @threads for beta in [BETA_SIM for z in 1:nthreads()]
        h_tag = "$(Dates.format(now(), "mm_dd"))_$(Threads.threadid())"
        #beta=BETA_SIM
        #h_tag = "$(Dates.format(now(), "mm_dd"))"
        if beta == 1 || beta==0
            beta = Int(beta)  #For some reason I need this for correct saving name
        end
        beta_str = string(beta) #convert beta to string for convenience in saving the
        beta_str = replace(beta_str, "." => "")

        #Initialization. Hay que inicializar estas variables porque si las inicilaizo en el loop luego se pierden (cosas de Julia)
        Temp = Init_T
        Results = Dict()
        init_times_0 = nothing
        periodicity_0 = nothing
        Schedule_0 = nothing
        Interventions_0 = nothing
        c = nothing
        p = nothing
        prob = nothing
        cb = nothing
        J0 = nothing

        #Compute the INT value of the uncontrolled system
        c  = generate_c(Dict())  #Interventions=Dict() -> No interventions
        p = (Dict(), c) 
        prob = ODEProblem(fun!, u0, tspan, p)
        cb = DiscreteCallback(condition, affect!, save_positions=(true,true))   # Define the callback
        sol = solve(prob, alg, callback=cb, abstol=tol, maxiters=max_iters)   # Solve the ODE problem with the callback
        norm = trapz(sol.t,[ sum(sol(t)) for t in sol.t ])


        while true # Generar una condición inicial no muy mala
            init_times_0, periodicity_0, Schedule_0 = random_Schedule()
            Interventions_0 = generate_Interventions(Schedule_0)

            c  = generate_c(Interventions_0)  #no entiendo del todo el global, pero si no no funciona. Creo que es porque c, W no se pasa como argumento a la funcion "condition"
            p = (Interventions_0, c) #la p tiene que estar aunque no se utilize
            prob = ODEProblem(fun!, u0, tspan, p)
            cb = DiscreteCallback(condition, affect!, save_positions=(true,true))   # Define the callback
            sol = solve(prob, alg, callback=cb, abstol=tol, maxiters=max_iters)   # Solve the ODE problem with the callback

            J0 = Cost_function(sol, periodicity_0, beta, norm)
            if J0[1]<0.6
                Results[0] = [J0, init_times_0, periodicity_0, Temp, true]
                break
            end
        end

        Results = Simulated_Annealing(beta, beta_str, Temp, init_times_0, periodicity_0, J0, c, prob, p, cb, Results, norm, h_tag)
        save_object(joinpath(@__DIR__, "Results/Simulated_Annealing_$(year)_b$(beta_str)_$(model)_$(h_tag)"), Results)

    end
end

Run_Sim(BETA_SIM)

