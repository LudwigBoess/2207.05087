# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO
using PyPlot, PyPlotUtility
using GadgetUnits
using AnalyticMHDTestSolutions
using StatsBase
using Statistics
using Printf
using DelimitedFiles
using LinearAlgebra


"""
    read_eta(fi, mach_ideal, dsa_model)

Reads pre-computed efficiencies. Solving the Riemann problem for every Mach number and efficiency model is inpractical.
"""
function read_eta(fi, mach_ideal, dsa_model)

    ideal_data = readdlm("data/tests/injection/eta_M/etaM_analytic.txt")

    mach_arr = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
        10.0, 15.0,
        20.0, 30.0, 40.0, 50.0,
        60.0, 70.0, 80.0, 90.0, 100.0]

    x_lim = [85.0, 87.0]

    x = read_block(fi, "POS", parttype = 0)[1, :]

    k = findall(x_lim[1] .< x .<= x_lim[2])

    rho = read_block(fi, "RHO", parttype = 0)[k]
    u = read_block(fi, "U", parttype = 0)[k]

    γ_cr = 4 / 3

    crpp = read_block(fi, "CRpP", parttype = 0)[k]
    CRpE = crpp ./ ((γ_cr - 1.0) .* rho)

    η_ideal = ideal_data[:, dsa_model+1][mach_arr.==mach_ideal][1]
    U_E_ratio = mean(CRpE ./ u)

    return η_ideal, U_E_ratio
end


"""
    plot_mach_efficiencies_compare(mach_arr, sim_paths, dsa_models, model_names, acc_colors, plot_name)

Main plot routine
"""
function plot_mach_efficiencies_compare(mach_arr, sim_paths, dsa_models, model_names, acc_colors, plot_name)

    lw = 2.5
    x_pixels = 600

    fig = get_figure(1.1; x_pixels)

    #fig = plt.Figure()

    plot_styling!(x_pixels)

    gs = plt.GridSpec(4, 1, figure = fig)

    subplot(get_gs(gs, 0:3, 0))

    ax1 = gca()
    ax1.set_xlim([1.0, 150.0])
    ax1.set_ylim([1.e-3, 2.0])
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    axis_ticks_styling!(ax1)

    ax1.set_xticklabels([])

    ylabel("Downstr. Energy Ratio  " * L"E_{cr,2}/E_{th,2}")

    subplot(get_gs(gs, 3, 0))

    ax2 = gca()
    ax2.set_xlim([1.0, 150.0])
    ax2.set_ylim([1.e-3, 2.0])
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    axis_ticks_styling!(ax2)

    xlabel("Sonic Mach Number  " * L"\mathcal{M}_s")

    ylabel("Rel. Error")

    for m ∈ 1:length(dsa_models)

        @info "DSA model: $(model_names[m])"
        # read data
        η_ideal = Vector{Float64}(undef, size(mach_arr, 1))
        U_E_ratio_S = Vector{Float64}(undef, size(mach_arr, 1))
        U_E_ratio_vs = Vector{Float64}(undef, size(mach_arr, 1))

        for i = 1:size(mach_arr, 1)
            @info "\tMach $(mach_arr[i])"
            fi = sim_paths[m] * "Mach_$(@sprintf("%i", mach_arr[i]))_0/snap_000"
            η_ideal[i], U_E_ratio_S[i] = read_eta(fi, mach_arr[i], dsa_models[m])

        end

        L1_errors_S = @. abs((U_E_ratio_S - η_ideal) / η_ideal)
        L1_errors_vs = @. abs((U_E_ratio_vs - η_ideal) / η_ideal)

        ax1.plot(mach_arr, U_E_ratio_S, color = acc_colors[m]; lw)#, label = model_names[m]; lw)

        ax1.plot(mach_arr, η_ideal, color = "white", alpha = 0.5, lw = 2lw)
        ax1.plot(mach_arr, η_ideal, color = "k", ls = ":", alpha = 0.8; lw)

        ax2.plot(mach_arr, L1_errors_S, color = acc_colors[m]; lw)

    end


    for m ∈ 1:length(dsa_models)
        ax1.plot([NaN], [NaN], color = acc_colors[m], label = model_names[m]; lw)
    end
    ax1.legend(frameon = false, loc = "lower right")

    subplots_adjust(hspace = 0, wspace = 0)

    savefig(plot_name, bbox_inches = "tight")
    close(fig)
end



mach_arr = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
    10.0, 15.0,
    20.0, 30.0, 40.0, 50.0,
    60.0, 70.0, 80.0, 90.0, 100.0]

acc_colors = ["#000c8f", "#a700c9", "#c90000", "#f28500"]
model_names = ["Kang+ 2007", "Kang&Ryu 2013", "Caprioli&Spitkovsky 2014", "Ryu+ 2019"]
dsa_models = [0, 1, 3, 2]

general_path = "data/tests/injection/eta_M/"
dsa_model_names = [ "Kang07/", "KR13/",
                    "CS14/",   "Ryu_19/"]

sim_paths = general_path .* dsa_model_names

plot_name = "Plots/Fig05.pdf"
plot_mach_efficiencies_compare(mach_arr, sim_paths, dsa_models, model_names, acc_colors, plot_name)

