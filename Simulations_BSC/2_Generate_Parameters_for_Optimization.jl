import Pkg

#Pkg.add("DifferentialEquations")
#Pkg.add("GeoStats")

using DifferentialEquations
using Plots
using CSV
using Interpolations
using Dates
using JLD2
using Statistics
using DataFrames
using XLSX
using DelimitedFiles


#SOME FUNCTIONS
function smooth_data(data, window_size)
    smoothed_data = similar(data, length(data) - window_size + 1)
    for i in 1:length(smoothed_data)
        smoothed_data[i] = mean(@view data[i:i+window_size-1])
    end
    return smoothed_data
end

#INTERPOLATE DATA SPATIALLY AND TEMPORALLY

#Create a function to interpolate spatially using Inverse Distance Squares
#The function needs a dataframe containig df.longitude, df.latitude,
#df.time with time in float vlaues and df.values with values being the column 
#containig the precipitation, temperature or whatever
function IDW_Interp(long, lat, t, df)
    # Find indices where time matches t
    time_indices = df.time .== t

    # Calculate distances
    distances = sqrt.((long .- df.longitude).^2 .+ (lat .- df.latitude).^2)
    
    # Filter out NaN values and corresponding distances
    valid_indices = .!ismissing.(df.values) .& time_indices
    distances = distances[valid_indices]
    values = df.values[valid_indices]

    # Check if any distance is 0
    if any(distances .== 0)
        return values[distances .== 0][1]  # Return the corresponding value
    else
        weights = 1 ./ distances
        weighted_values = weights .* values
        interpolated_values = sum(weighted_values) ./ sum(weights)
        return interpolated_values
    end
end


#Create a function to interpolate Data
#The function needs a dataframe containig df.longitude, df.latitude,
#df.time with time in float vlaues and df.values with values being the column 
#containig the precipitation, temperature or whatever
function Data_Interp(long, lat, t, df)
    # Check if t is in the dataframe
    if t in df.time
        return IDW_Interp(long, lat, t, df)
    else
        # Find the two closest time values in df.time such that t is in the middle
        sorted_times = sort(df.time)
        lower_idx = searchsortedlast(sorted_times, t)
        upper_idx = searchsortedfirst(sorted_times, t)

        lower_time = sorted_times[lower_idx]
        upper_time = sorted_times[upper_idx]
        
        # Compute the weights for linear interpolation
        lower_weight = (upper_time - t) / (upper_time - lower_time)
        upper_weight = 1 - lower_weight

        # Compute IDW_Interp at these two time values and perform linear interpolation
        lower_val = IDW_Interp(long, lat, lower_time, df)
        upper_val = IDW_Interp(long, lat, upper_time, df)
        interpolated_value = lower_weight * lower_val + upper_weight * upper_val
        
        return interpolated_value
    end
end


#MOSQUITO FUNCTIONS
## Thermal responses Aedes Albopictus from Mordecai 2017:
#Basic Functions
function Briere_func(cte, tmin, tmax, temp)
    if tmax >= temp && tmin<=temp
        outp = temp * cte * (temp - tmin) * sqrt(tmax - temp)
    else
        return 0.00001
    end
    return outp
end

function Quad_func(cte, tmin, tmax, temp)
    outp = cte * (temp - tmin) * (tmax-temp)
    if temp<=tmin || temp>=tmax
        return 0.00001
    end
    return outp
end

function QuadN_func(cte, c1, c2, temp)
    outp = cte * temp^2 + c1 * temp + c2
    if outp < 0.00001
        outp = 0.00001
    end
    return outp
end

function Lin_func(m, z, temp)
    outp = -m*temp + z
    if outp < 0.00001
        outp = 0.00001
    end
    return outp
end


models=["CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "MIROC6", "MPI-ESM1-2-HR", "NorESM2-MM"]

for model in models

    println(model)

    #LOAD DATA 
    # Open the Excel file for reading
    Meteo = CSV.read(joinpath(@__DIR__, "Dades_netes/SSP245/Meteo_SSP245_31_60_$(model).csv"), DataFrame);
    datestamp = map(DateTime, Meteo.Fecha)
    Meteo.year = Dates.year.(Meteo.Fecha)


    #PREPARE PARAMETERS FOR EACH YEAR
    for year in 2031:2060

        # Filter Meteo for the current year
        Meteo_year = filter(row -> row.year == year, Meteo)
        datestamp_year = filter(x -> Dates.year(x) == year, datestamp)


        d0=Dates.DateTime("$(year)-03-01", "yyyy-mm-dd")
        t = map(x -> 1 + Dates.value(x - d0)/(24*60*60*1000) , datestamp_year)  #convert dates to time vector (in days)
        Meteo_year.time = t


        #Create interpolated objects for weather
        X_long = sort(Float64.(unique(Meteo_year.longitude)))
        Y_lat = sort(Float64.(unique(Meteo_year.latitude)))
        t_interp=unique(t)

        #Rainfall
        df_precip = select(Meteo_year, [:longitude, :latitude, :time, :Rainfall]) #create df for precipitation
        rename!(df_precip, :Rainfall => :values)
        precip_data = [Data_Interp(i, j, k, df_precip) for i in X_long, j in Y_lat, k in t_interp]
        Rainfall = LinearInterpolation((X_long,Y_lat,t_interp), precip_data, extrapolation_bc=Flat())  #Out of bounds, extrapolate with constant closest value


        #21Days Cumulative Rainfall
        df_precip_21 = select(Meteo_year, [:longitude, :latitude, :time, :Rainfall_21]) #create df for precipitation
        rename!(df_precip_21, :Rainfall_21 => :values)
        precip_data_21 = [Data_Interp(i, j, k, df_precip_21) for i in X_long, j in Y_lat, k in t_interp]
        Rainfall_21 = LinearInterpolation((X_long,Y_lat,t_interp), precip_data_21, extrapolation_bc=Flat())  #Out of bounds, extrapolate with constant closest value


        #Mean Temperature
        df_mean_temp = select(Meteo_year, [:longitude, :latitude, :time, :Mean_T]) #create df for precipitation
        rename!(df_mean_temp, :Mean_T => :values)
        mean_temp_data = [Data_Interp(i, j, k, df_mean_temp) for i in X_long, j in Y_lat, k in t_interp]
        Mean_Temp =  LinearInterpolation((X_long,Y_lat,t_interp), mean_temp_data, extrapolation_bc=Flat())
        #print(df_mean_temp.values[5601:5700])


        #Settings for smoothing interpolations
        window=7
        t_smooth=LinRange(0,last(t_interp),length(t_interp)-window+1)

        mean_temp_data_smooth=Array{Float64}(undef,length(X_long),length(Y_lat),length(t_smooth))
        for i in 1:length(X_long)
            for j in 1:length(Y_lat)
                mean_temp_data_smooth[i,j,:]=smooth_data(Mean_Temp(X_long[i],Y_lat[j],t_interp),window)
            end
        end
        Mean_Temp_smooth =LinearInterpolation((X_long,Y_lat,t_smooth), mean_temp_data_smooth , extrapolation_bc=Flat())


        save_object(joinpath(@__DIR__, "Parameters/Rainfall_$(year)_$(model)"), Rainfall)
        save_object(joinpath(@__DIR__, "Parameters/Rainfall_21_$(year)_$(model)"), Rainfall_21)
        save_object(joinpath(@__DIR__, "Parameters/Mean_Temperature_$(year)_$(model)"), Mean_Temp_smooth)



        #2019 Mordecai
        bit_rate_alb(temp) = Briere_func(0.00019, 10.4, 38.1, temp)  # Biting rate
        fecundity_alb(temp) = Briere_func(0.0477, 7.9, 35.6, temp)  # Fecundity (eggs/female/gonotrophic cycle)
        lf_alb(temp) = Quad_func(1.39,13.5,31.4, temp) #lifespan, 1/death rate
        pEA_alb(temp) = Quad_func(0.00356, 9.1, 36.2, temp) #egg to adult survival



        #Basic parameters Culex
        bit_rate_cul(temp) = Briere_func(0.00017, 9.4, 39.6, temp)  # Biting rate
        fecundity_cul(temp) = Quad_func(0.598, 5.3, 38.9, temp)  # Fecundity (eggs/female/gonotrophic cycle)
        lf_cul(temp) = Lin_func(4.86, 169.8, temp)  # Adult life span
        pLA_cul(temp) = Quad_func(0.0036, 7.8, 38.4, temp)
        pEL_cul(temp) = Quad_func(0.00211, 3.2, 42.6, temp)  # Egg development rate


        #R0 Cul
        function R0_cul(temp)

            a = bit_rate_cul(temp)
            f = (1 / 2) * fecundity_cul(temp)
            deltaa = lf_cul(temp)
            probel = pEL_cul(temp)
            probla = pLA_cul(temp)

            R0 = ((f * a * deltaa) * probla * probel)^(1 / 3)

            return R0

        end



        #Interpolate Final parameters R0_culex
        R0_Cul=LinearInterpolation((X_long,Y_lat,t_smooth), R0_cul.(Mean_Temp_smooth(X_long,Y_lat,t_smooth)), extrapolation_bc=Flat())

        #Save
        save_object(joinpath(@__DIR__, "Parameters/R0_culex_$(year)_$(model)"), R0_Cul) # Culex pipiens

        #SAVE D MATRIX AS PARAMETER
        D = readdlm(joinpath(@__DIR__, "Parameters/ZR_Distance_Matrix.txt"), '\t')

        save_object(joinpath(@__DIR__, "Parameters/D_Matrix"), D)


        #SAVE COORDINATES AND PROFUNDITAT AS PARAMETERS
        Unique_BS = CSV.read(joinpath(@__DIR__, "Dades_netes/Unique_Items.csv"), DataFrame) 
        XY = hcat(Unique_BS[!, :x_lon], Unique_BS[!, :y_lat]) # Create a matrix with XY info
        X_emb=Dict()
        for (i,value) in enumerate(Unique_BS.x_lon)
            X_emb[i]=value
        end
        Y_emb=Dict()
        for (i,value) in enumerate(Unique_BS.y_lat)
            Y_emb[i]=value
        end

        #Store Prof_Sorers as Dictionary
        Prof_Sorrers=Dict()
        for (i, value) in enumerate(Unique_BS.prof_sorrer)
            Prof_Sorrers[i] = Float64(value)
        end

        save_object(joinpath(@__DIR__, "Parameters/X_coords"), X_emb)
        save_object(joinpath(@__DIR__, "Parameters/Y_coords"), Y_emb)
        save_object(joinpath(@__DIR__, "Parameters/Profunditat_Sorrers"), Prof_Sorrers)


    end
end
