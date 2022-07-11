# activate the environment
using Pkg; Pkg.activate(".")

using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using AnalyticMHDTestSolutions
using StatsBase
using Statistics
using LaTeXStrings
using Printf
using DelimitedFiles
using LinearAlgebra

"""
    get_ideal_solution_cluster_shock(t)

Calculates the analytic solution at time `t`.
"""
function get_ideal_solution_cluster_shock(t)

    fi = "data/tests/injection/NLDSA/snap_IC"
    rho = read_block(fi, "RHO", parttype = 0)
    u = read_block(fi, "U", parttype = 0)
    rhol = Float64(maximum(rho))
    rhor = Float64(minimum(rho))
    Ul = Float64(maximum(u))
    Ur = Float64(minimum(u))
    x_contact = 140.0

    dsa_model = 2
    xs_first_guess = 3.3
    Mach = 5.541278364096982
    par = RiemannParameters(; rhol, rhor, Ul, Ur, t, x_contact, dsa_model, xs_first_guess, Mach)

    x_ideal = collect(LinRange(0.0, 280.0, 2))
    sol = AnalyticMHDTestSolutions.solve(x_ideal, par)

    return sol, par
end


"""
    get_relative_slope_histograms(fi, q_ideal)

Constructs histograms of slopes in the spectrum, relative to the ideal slope `q_ideal`.
"""
function get_relative_slope_histograms(fi, q_ideal)

    h  = read_header(fi)
    sol, par = get_ideal_solution_cluster_shock(h.time)

    pos_contact =  sol.v34 * par.t + par.x_contact
    pos_shock   =  sol.vs  * par.t + par.x_contact
    pos_downstream = (pos_contact + pos_shock) * 0.5
    width_downstream = (pos_shock - pos_contact )

    x_range = [ pos_downstream - 0.5width_downstream, pos_downstream + 0.5width_downstream]

    x = read_snap(fi, "POS", 0)[1,:]

    k = findall( x_range[1] .<= x .< x_range[2] )

    crp_slope = read_snap(fi, "CRpS", 0)[:,k]
    slopes = reduce(vcat, crp_slope)

    println("mean slopes: $(mean(slopes)) ± $(std(slopes))")

    slopes ./= q_ideal

    bins = LinRange(0.9, 1.1, 40)

    h = fit(Histogram, slopes, bins)
    counts = h.weights ./ (size(k, 1) * size(crp_slope, 1))

    bin_centers = Vector{Float64}(undef, size(bins,1)-1)
    for i = 1:size(bins,1)-1
        bin_centers[i] = 0.5*( bins[i] + bins[i+1])
    end

    return bin_centers, counts
end



"""
    read_to_dict(fi, field_names, fields, id)

Helper function to write the quantities of a specific particle to a `Dict`
"""
function read_to_dict(fi, field_names, fields, id)

    snap_info = read_info(fi)

    d = Dict{String, Union{Vector{T}, T} where T}()

    for i = 1:size(fields, 1)
        dim   = snap_info[getfield.(snap_info, :block_name) .== fields[i]][1].n_dim
        if dim == 1
            d[field_names[i]] = read_block(fi, fields[i], parttype=0)[id]
        else
            d[field_names[i]] = read_block(fi, fields[i], parttype=0)[:,id]
        end
    end

    return d
end

"""
    read_up_down(fi, fields, field_names, 
                 fields_shocked, field_names_shocked)

Read up- and downstream quantities of the shock.
"""
function read_up_down(fi, fields, field_names, 
                          fields_shocked, field_names_shocked)

    h  = read_header(fi)
    sol, par = get_ideal_solution_cluster_shock(h.time)

    GU = GadgetPhysical()

    pos_contact =  sol.v34 * par.t + par.x_contact
    pos_shock   =  sol.vs  * par.t + par.x_contact
    pos_upstream   =  1.2 * sol.vs  * par.t + par.x_contact
    pos_downstream = (pos_contact + pos_shock) * 0.5

    x       = read_block(fi, "POS", parttype=0)[1,:] 
    down_id = findfirst( x .≈ pos_downstream)
    up_id   = findfirst( x .≈ pos_upstream)

    println(pos_downstream, " ", pos_upstream)
    k       = findall(pos_downstream .< x .< pos_upstream )
    mach    = read_block(fi, "MACH", parttype=0)[k]
    max_id  = findmax(mach)[2]

    shock_down       = read_to_dict(fi, field_names, fields, down_id)
    shock_down["vA"] =  √( sum( @views shock_down["bfld"][i]^2 for i = 1:3) / ( 4π * shock_down["rho"][1] * GU.rho_cgs) ) / GU.v_cgs

    shock_up       = read_to_dict(fi, field_names, fields, up_id)
    shock_up["vA"] = √( sum(shock_up["bfld"][i]^2 for i = 1:3) / ( 4π * shock_up["rho"][1] * GU.rho_cgs) ) / GU.v_cgs

    # read properties from shockfinder
    shock_max  = read_to_dict(fi, field_names_shocked, fields_shocked, k[max_id])
    
    return shock_down, shock_up, shock_max
end

"""
    analytic_q(shock_up, shock_down)

Analytic solutions for DSA and NLDSA slopes.
"""
function analytic_q(shock_up, shock_down)

    xs =  shock_down["rho"] / shock_up["rho"]
    α  = shock_down["vA"] / norm(shock_down["v"])

    q_dsa   = 3xs / ( xs - 1.0 )
    q_nldsa = 3xs / ( xs - 1.0 - α )

    return q_dsa, q_nldsa

end

"""
    plot_slope_histograms_relative(files, q_dsa, q_nldsa, labels, colors, plot_name)

Main plot function.
"""
function plot_slope_histograms_relative(files, q_dsa, q_nldsa, labels, colors, plot_name)

    lw = 2.5

    x_pixels = 700
    fig = get_figure(1.0; x_pixels)
    plot_styling!(x_pixels)

    ax = gca()

    axis_ticks_styling!(ax)

    xlabel("Rel. Difference of Spectral Slope  " * L"q/q_{\mathrm{ideal}}")
    ylabel("Probability  " * L"N/N_{\mathrm{tot}}")

    ax.set_xlim([0.95, 1.05])
    ax.set_ylim([0.0, 1.0])

    for i = 1:size(files, 1)
        @info "$(labels[i])"
        if i ∈ [1]
            q_id = q_dsa
        else
            q_id = q_nldsa
        end
        bins, counts_dsa = get_relative_slope_histograms(files[i], q_id)
        step(bins, counts_dsa,
            color = colors[i], label = labels[i], lw = lw)
    end

    legend(frameon = false, loc = "upper right")

    savefig(plot_name, bbox_inches = "tight")

    close(fig)

end

folders = [ "dsa",
            "dsa_recalc",
            "nldsa",
            "nldsa_recalc"
            ]

colors = ["#0207ab", "teal", "red", "orange", "#e31607", "#a6009a", "#5e0858"]

fi = "data/tests/injection/NLDSA/" .* folders .* "/snap_000"

labels = ["DSA", "DSA+recalc", 
          "NLDSA", "NLDSA+recalc"
          ]

plot_name = "Plots/Fig08.png"


fields      = [ "RHO", "U", "BFLD", "VEL" ]
field_names = [ "rho", "u", "bfld", "v"]

fields_shocked      = [ "RHO", "U", "BFLD", "VEL", "SHVD", "SHVU", "MALF", "SHAD", "SHAU", 
                        "SHCP", "MACH", "SHRH", "SHPD", "SHPU", "SHSP" ]
field_names_shocked = [ "rho", "u", "bfld", "v", 
                        "v_down", "v_up", "alfven_mach", "alfven_down", "alfven_up", 
                        "shock_compress", "mach", "rho_up", "p_down", "p_up", "v_sh"]

# calculate the analytic solutions
shock_down, shock_up, shock_max = read_up_down(fi[1], fields, field_names, fields_shocked, field_names_shocked)
q_dsa, q_nldsa = analytic_q(shock_up, shock_down)

plot_slope_histograms_relative(fi, q_dsa, q_nldsa, labels, colors, plot_name)