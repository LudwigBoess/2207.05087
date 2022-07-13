# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using ProgressMeter
using SpectralCRsUtility
using Base.Threads
using Printf
using Cosmology

# edit this to where you stored the simulation snapshots
const simulation_path = "/path/to/snapshots/"

"""
    get_synch_power(data, h; 
                    pmin = 1.e-1, pmax = 1.e5)

Compute the synchrotron power within the cutout frame for B, 10B, 100B and 1muG.
"""
function get_synch_power(data, h; 
                        pmin = 1.e-1, pmax = 1.e5)

    @info "setting up spectra"

    GU = GadgetPhysical(h)

    # number of datapoint for spectrum
    Nnu = 20
    # set observational frequencies
    ν_arr = 10.0.^(LinRange(8, 11, Nnu))

    # cr setup 
    Nbins = size(data["CReN"],1)
    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)

    Npart = size(data["CReN"],2)

    # allocate arrays
    # synchrotron power
    P   = Matrix{Float64}(undef, 4, Npart)
    # synchrotron emissivity
    j_ν = Matrix{Float64}(undef, 4, Npart)

    @info "calculating synch emissivity"

    # convert to cgs
    m   = data["MASS"] .* GU.m_cgs
    rho = data["RHO"]  .* GU.rho_cgs

    # loop over all required frequencies
    @showprogress for iNu = 1:Nnu

        ν = ν_arr[iNu]

        # loop over all particles
        @sync begin
            @inbounds for i = 1:Npart
                Base.Threads.@spawn begin

                    # convert to physical units
                    norm   = GU.CR_norm .* 10.0.^data["CReN"][:,i]
                    slope  = Float64.(data["CReS"][:,i])
                    cut    = Float64(data["CReC"][i])

                    # compute magnetic field
                    B      = √( data["BFLD"][1,i]^2 + 
                                data["BFLD"][2,i]^2 + 
                                data["BFLD"][3,i]^2 )

                    # compute emissivities
                    j_ν[1, i] = synchrotron_emission(norm, slope, cut, B, par, ν0 = ν, 
                                            reduce_spectrum = true,
                                            integrate_pitch_angle = true)

                    j_ν[2, i] = synchrotron_emission(norm, slope, cut, 10B, par, ν0 = ν, 
                                            reduce_spectrum = true,
                                            integrate_pitch_angle = true)

                    j_ν[3, i] = synchrotron_emission(norm, slope, cut, 100B, par, ν0 = ν, 
                                            reduce_spectrum = true,
                                            integrate_pitch_angle = true)

                    j_ν[4, i] = synchrotron_emission(norm, slope, cut, 1.e-6, par, ν0 = ν, 
                                            reduce_spectrum = true,
                                            integrate_pitch_angle = true)
                end
            end
        end

        # compute synchrotron power in [W / Hz]
        P[1, iNu] = sum( @. j_ν[1,:] * m / rho ) * 1.e-7
        P[2, iNu] = sum( @. j_ν[2,:] * m / rho ) * 1.e-7
        P[3, iNu] = sum( @. j_ν[3,:] * m / rho ) * 1.e-7
        P[4, iNu] = sum( @. j_ν[3,:] * m / rho ) * 1.e-7

    end # Nu

    @info "done"

    return ν_arr, P
end

"""
    write_synch_spectrum(filename, ν_arr, P)

Write the synchrotron spectrum to a binary file.
"""
function write_synch_spectrum(filename, ν_arr, P)

    f = open(filename, "w")
    
    write(f, Int64(length(ν_arr)))
    write(f, ν_arr)
    for i = 1:length(ν_arr)
        write(f, P[1, i])
    end
    for i = 1:length(ν_arr)
        write(f, P[2, i])
    end
    for i = 1:length(ν_arr)
        write(f, P[3, i])
    end
    for i = 1:length(ν_arr)
        write(f, P[4, i])
    end

    close(f)
end


"""
    compute_synch_power(sim_path, center_pos, xyz_size, snap, P_filename)

Cuts out a cube around the `center_pos` of the relic and computes the synchrotron power contained in that cube.
"""
function compute_synch_power(sim_path, center_pos, xyz_size, snap, P_filename)

    filename = sim_path * "snap_$(@sprintf("%03i", snap))"

    @info "reading data"
    h = read_header(filename)

    blocks = ["POS", "HSML", "RHO", "U", "BFLD", "CReN", "CReS", "CReC"]

    cube = GadgetCube( center_pos - [0.5xyz_size, 0.5xyz_size, 0.5xyz_size], 
                       center_pos + [0.5xyz_size, 0.5xyz_size, 0.5xyz_size] )

    data = read_particles_in_geometry(filename, blocks, cube, use_keys=false)

    data["MASS"] = h.massarr[1] .* ones(length(data["HSML"])) .|> Float32

    ν_arr, P = get_synch_power(data, h)

    write_synch_spectrum(P_filename, ν_arr, P)

end


folders = [ "KR13d_thetaB", 
            "KR13t", 
            "KR13t_thetaB", 
            "Ryu19t",
            "Ryu19t_thetaB",
            "Ryu19t_thetaB_pinj", 
            "Ryu19t_thetaB_pinj_q"
            ]

"""
    Northern Relic
"""

# snap = 23


# for folder ∈ folders
#     sim_path = "$suite_path/$folder/"
#     P_filename = sim_path * "derived_data/NR_spectrum.dat"
#     filename = sim_path * "snap_$(@sprintf("%03i", snap))"
#     h = read_header(filename)
#     center_pos = h.boxsize .* [0.3, 0.475, 0.5]
#     xyz_size = 0.3h.boxsize
#     compute_synch_power(sim_path, center_pos, xyz_size, snap, P_filename)
# end

# error("done!")

"""
    Southern Relic
"""

# snap = 20

# for folder ∈ folders
#     sim_path = "$suite_path/$folder/"
#     P_filename = sim_path * "derived_data/SR_spectrum.dat"
#     filename = sim_path * "snap_$(@sprintf("%03i", snap))"
#     h = read_header(filename)
#     center_pos = h.boxsize .* [0.5, 0.475, 0.5]
#     xyz_size = 0.2h.boxsize
#     compute_synch_power(sim_path, center_pos, xyz_size, snap, P_filename)
# end

# error("done!")


"""
    Plot
"""

using PyPlot, PyPlotUtility

function read_spectrum(filename)

    f = open(filename, "r")
    Nnu = read(f, Int64)

    ν_arr    = read!(f, Vector{Float64}(undef, Nnu)) 
    P_sim    = read!(f, Vector{Float64}(undef, Nnu))
    P_10sim  = read!(f, Vector{Float64}(undef, Nnu))
    P_100sim = read!(f, Vector{Float64}(undef, Nnu)) 
    P_1muG   = read!(f, Vector{Float64}(undef, Nnu))

    close(f)

    return ν_arr, P_sim, P_10sim, P_100sim, P_1muG

end

"""
    Obs. Data: Stroe A. et al., 2016, MNRAS, 455, 2402, Table 3
"""

const ν_obs = [0.15,  0.325,  0.61,   1.2,   1.4,  1.7, 2.25,  2.3, 2.75, 3.25, 3.75, 16.0, 30.0] .* 1.e9
const Jν    = [780.4, 315.7, 222.3, 125.7, 117.3, 91.2, 61.0, 54.3, 50.0, 36.1, 27.2,  1.2,  0.6]
const σJν   = [80.0,   32.4,  22.4,   12.6, 11.8,  9.2,  3.6,  5.6,  2.7,  1.9,  2.6,  0.5,  0.3]

const z = 0.1921



"""
    plot_synch_spectrum(filenames, names, colors, plot_name)

Main plot routine.
"""
function plot_synch_spectrum(filenames, names, colors, plot_name)

    Nfiles = length(filenames)

    c = cosmology(h=0.674,
                 OmegaM=0.315,
                 Neff=2.99)

    fig = get_figure()
    plot_styling!()

        ax = gca()
        axis_ticks_styling!(ax)
        ax.set_xlim(1.e8, 1.e11)
        ax.set_ylim(1.e-1, 4.e3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        xlabel("Obs. Frequency  " * L"\nu" * " [Hz]")
        ylabel("Synchrotron Surface Brightness  " * L"S_{\mathrm{syn}, \nu_0}" * " [mJy]")

        for i = 1:Nfiles 
            ν_arr, P_sim, P_10sim, P_100sim, P_1muG = read_spectrum(filenames[i] * "/synchrotron/NR_spectrum.dat" )

            P_in_mJy = P_10sim ./ mJy_to_W(c, 1.0, z)
            P_norm = 1.e3 / P_in_mJy[1]
            
            plot(ν_arr, P_in_mJy, color=colors[i], label=names[i])
            plot(ν_arr, P_norm .* P_in_mJy, color=colors[i], linestyle="--")
        end

        # observations
        errorbar(ν_obs, Jν, yerr=σJν, color = "k", fmt="D", label="Stroe+16", markersize=3.0,
                markerfacecolor="w", markeredgecolor="k")

        # dummy lines for legend
        plot([0.0], [0.0], color="k", label=L"10 B_\mathrm{sim}")
        plot([0.0], [0.0], color="k", linestyle="--", label="Normalized")

        # legends
        handles, labels = ax.get_legend_handles_labels()
        ax.add_artist(legend(handles[1:Nfiles], labels[1:Nfiles], frameon=false, loc="upper right"))
        ax.add_artist(legend(handles[Nfiles+1:end], labels[Nfiles+1:end], frameon=false, loc="lower left"))


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

filenames = ["data/cluster_mergers/$folder" for folder ∈ folders]

names = ["KR13" * L"d \theta_B",
         "KR13" * L"t", 
         "KR13"  * L"t \theta_B",
         "Ryu19" * L"t",
         "Ryu19" * L"t \theta_B",
         "Ryu19" * L"t \theta_B p_\mathrm{inj}",
         "Ryu19" * L"t \theta_B p_\mathrm{inj} q_\alpha"
         ]
colors = ["#fc17f1", "purple", "darkblue", "dodgerblue", "teal", "#82f70c", "#dbc202"]

plot_name = "Plots/Fig12.pdf"
plot_synch_spectrum(filenames, names, colors, plot_name)

