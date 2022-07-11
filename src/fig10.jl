# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using DSAModels
using PyPlot, PyPlotUtility
using Statistics
using Unitful
using Distributions
using StatsBase

# edit this to where you stored the simulation snapshots
const simulation_path = "/path/to/snapshots"

"""
    filter_shocked_particles(folder::String, relic::String)

Reads all particles currently identified as shocked.
"""
function filter_shocked_particles(folder::String, relic::String)

    sim_path = "$simulation_path/$folder/"

    if relic == "NR"
        fi = sim_path * "snap_023"
        h = read_header(fi)
        center_pos = h.boxsize .* [0.3, 0.475, 0.5]
        xyz_size = 0.3h.boxsize
    else
        fi = sim_path * "snap_020"
        h = read_header(fi)
        center_pos = h.boxsize .* [0.5, 0.475, 0.5]
        xyz_size = 0.2h.boxsize
    end

    shock =  GadgetCube(center_pos - [0.5xyz_size, 0.5xyz_size, 0.5xyz_size], 
                        center_pos + [0.5xyz_size, 0.5xyz_size, 0.5xyz_size])

    pos = read_block(fi, "POS", parttype = 0)

    sel = GadgetIO.get_geometry_mask(shock, pos)

    mach = read_block(fi, "MACH", parttype = 0)[sel]

    sel_mach = findall(mach .> 1.0)
    mach = mach[sel_mach]

    Bfld = read_block(fi, "BFLD", parttype = 0)[:,sel[sel_mach]]
    B = @. √(Bfld[1, :]^2 + Bfld[2, :]^2 + Bfld[3, :]^2)

    return Dict("M" => mach, "B" => B)
end

"""

Computes the histograms for `M_s` and `B` of all shocked particles
"""
function get_histograms(folder)

    # store relevant particles
    NR = filter_shocked_particles(folder, "NR")
    SR = filter_shocked_particles(folder, "SR")

    # read Mach data

    mach_bins = LinRange(1.0, 7.0, 26)

    # mach histograms
    h = fit(Histogram, NR["M"], mach_bins)
    NM_NR = h.weights 

    h = fit(Histogram, SR["M"], mach_bins)
    NM_SR = h.weights 

    
    # B histograms
    B_bins = 10.0 .^ LinRange(-10, -5, 51)

    h = fit(Histogram, NR["B"], B_bins)
    NB_NR = h.weights 

    h = fit(Histogram, SR["B"], B_bins)
    NB_SR = h.weights 


    return Dict("Ms" => Dict("bins" => mach_bins,
            "NR" => NM_NR,
            "SR" => NM_SR),
        "B" => Dict("bins" => 1.e6 .* B_bins,
            "NR" => NB_NR,
            "SR" => NB_SR))

end

"""
    set_up_axis!(ax)

Helper function to get plot layout.
"""
function set_up_axis!(ax)

    # upper left
    subplot(2, 2, 1)
    ax[1] = gca()
    axis_ticks_styling!(ax[1])

    ax[1].set_yscale("log")
    ax[1].set_xlim([1.0, 7.0])
    ax[1].set_ylim([100, 2.e5])
    ax[1].set_ylabel("SPH particles")
    ax[1].set_xticklabels([])

    # upper right
    subplot(2, 2, 2)
    ax[2] = gca()
    axis_ticks_styling!(ax[2])

    ax[2].set_yscale("log")
    ax[2].set_xscale("log")
    ax[2].set_xlim([1.0e-3, 10.0])
    ax[2].set_ylim([10, 2.e5])
    ax[2].set_ylabel("SPH particles")
    ax[2].set_xticklabels([])

    # lower left
    subplot(2, 2, 3)
    ax[3] = gca()
    axis_ticks_styling!(ax[3])

    ax[3].set_yscale("log")
    ax[3].set_xlim([1.0, 7.0])
    ax[3].set_ylim([100, 2.e5])
    ax[3].set_ylabel("SPH particles")
    ax[3].set_xlabel("Sonic Mach Number " * L"\mathcal{M}_s")

    # lower right
    subplot(2, 2, 4)
    ax[4] = gca()
    axis_ticks_styling!(ax[4])

    ax[4].set_yscale("log")
    ax[4].set_xscale("log")
    ax[4].set_xlim([1.0e-3, 10.0])
    ax[4].set_ylim([10, 2.e5])
    ax[4].set_ylabel("SPH particles")
    ax[4].set_xlabel("Magnetic Field Strength B [μG]")

    ax
end


"""
    plot_histograms(folders, colors, labels, plot_name)

Main plot routine.
"""
function plot_histograms(folders, colors, labels, plot_name)

    Nfiles = length(folders)

    fig = get_figure(2.0)
    plot_styling!()

    ax = Vector{Any}(undef, 4)
    set_up_axis!(ax)

    for i = 1:Nfiles
        
        data = get_histograms(folders[i])

        # Mach PDF left
        ax[1].step(data["Ms"]["bins"][1:end-1], data["Ms"]["NR"], where = "post", color = colors[i], label = labels[i])
        # B PDF left
        ax[2].step(data["B"]["bins"][1:end-1], data["B"]["NR"], where = "post", color = colors[i], label = labels[i])

        # Mach PDF right 
        ax[3].step(data["Ms"]["bins"][1:end-1], data["Ms"]["SR"], where = "post", color = colors[i], label = labels[i])
        # B PDF right 
        ax[4].step(data["B"]["bins"][1:end-1], data["B"]["SR"], where = "post", color = colors[i], label = labels[i])

    end
    ax[3].legend(frameon = false, fontsize=11, loc="lower right")
    ax[1].text(6.0, 7.e4, "NR", fontsize = 20)
    ax[3].text(6.0, 7.e4, "SR", fontsize = 20)
    subplots_adjust( hspace= 0)

    savefig(plot_name, bbox_inches = "tight")
    close(fig)
end

folders = [ "KR13d_thetaB", 
            "KR13t", 
            "KR13t_thetaB",
            "Ryu19t",
            "Ryu19t_thetaB",
            "Ryu19t_thetaB_pinj", 
            "Ryu19t_thetaB_pinj_q"
            ]

colors = ["#fc17f1", "purple", "darkblue", "dodgerblue", "teal", "#82f70c", "#dbc202"]

labels = [ "KR13" * L"d \theta_B",
           "KR13" * L"t", 
           "KR13"  * L"t \theta_B",
           "Ryu19" * L"t",
           "Ryu19" * L"t \theta_B",
           "Ryu19" * L"t \theta_B p_\mathrm{inj}",
           "Ryu19" * L"t \theta_B p_\mathrm{inj} q_\alpha"]

plot_name = "Plots/Fig10.png"

plot_histograms(folders, colors, labels, plot_name)
