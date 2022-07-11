# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using SPHtoGrid


"""
CRep + CRpP + Xcr + Synch
"""
function plot_cr_total_compare(map_paths, relic)

    if relic == "NR"
        vmin_arr = [1.e-18, 1.e-18, 1.e-5, 1.e15, 1.e15, -2.0]
        vmax_arr = [1.e-12, 1.e-13, 1.e0, 1.e20, 1.e20, -0.5]
        contour_levels = [4.6]
        time = 2.25 # this is hardcoded so you can plot the maps without the snapshots
        plot_name = "Plots/Fig11.pdf"
        Nsnap = 23

    else
        vmin_arr = [1.e-18, 1.e-18, 1.e-5, 1.e13, 1.e13, -2.0]
        vmax_arr = [1.e-13, 1.e-13, 1.e0, 1.e18, 1.e18, -0.5]
        contour_levels = [3.0]
        time = 1.96
        plot_name = "Plots/Fig14.pdf"
        Nsnap = 20
    end

    filenames = ["CRpP", "CReP", "Xcr", "synchP_144MHz", "synchP_1.4GHz", "alpha_144MHz_1.4GHz"]
    Ncols = length(filenames)
    Nrows = length(map_paths)

    files = [ map_path * "$(relic)_$(@sprintf("%03i", Nsnap)).$filename.fits" 
                for map_path ∈ map_paths, filename ∈ filenames ]

    contour_files = [ map_path * "$(relic)_$(@sprintf("%03i", Nsnap)).Mach.fits" 
                for map_path ∈ map_paths, filename ∈ filenames ]


    contours = trues(Nrows*Ncols)

    pixel_cut = 50:950

    map_ref, par, snap, units = read_fits_image(files[1])

    map_arr = Vector{typeof(map_ref)}(undef, Nrows*Ncols)

    i_col = 0
    for i = 1:Nrows*Ncols
        map_arr[i], par, snap, units = read_fits_image(files[i])
        map_arr[i] = map_arr[i][:,pixel_cut]

        if (i-1)%Nrows == 0
            i_col += 1
        end

        if i_col > 3
            map_arr[i][map_arr[i] .< vmin_arr[i_col]] .= NaN
        else
            map_arr[i][map_arr[i] .< vmin_arr[i_col]] .= vmin_arr[i_col]
        end
    end

    par_arr = [par for _ = 1:Nrows*Ncols]

    contour_arr = Vector{typeof(map_ref)}(undef, Nrows*Ncols)
    contour_par_arr = [par for _ = 1:Nrows*Ncols]

    for i = 1:Nrows*Ncols
        contour_arr[i], par, snap, units = read_fits_image(contour_files[i])
        contour_arr[i] = contour_arr[i][:,pixel_cut]
    end

    im_cmap = [ "Reds_r", "Purples_r", "jet", "cubehelix_r", "cubehelix_r", "jet" ]

    cb_labels = [ L"P_{CR,p}" * " [erg cm" * L"^{-3}" *"]",
                L"P_{CR,e}" * " [erg cm" * L"^{-3}" *"]",
                L"X_{CR,p} = \frac{P_{CR,p}}{P_{th}}",
                L"P_{ν,144 \mathrm{ MHz}}" * " [W Hz" * L"^{-1}" * "]",
                L"P_{ν,1.4 \mathrm{ GHz}}" * " [W Hz" * L"^{-1}" * "]",
                L"\alpha_\mathrm{144 MHz}^\mathrm{1.4 GHz}"]
                
    annotate_time = falses(Nrows*Ncols)
    annotate_time[1] = true

    time_labels = [L"" for _ = 1:Nrows*Ncols ]
    time_labels[1] = L"t \: = " * "$(@sprintf("%0.2f", time)) Gyrs"

    annotate_time = trues(Nrows*Ncols)
    time_labels = ["" for _ = 1:Nrows*Ncols ]
    time_labels[1] = "KR13" * L"d \theta_B"
    time_labels[2] = "KR13" * L"t" 
    time_labels[3] = "KR13"  * L"t \theta_B"
    time_labels[4] = "Ryu19" * L"t"
    time_labels[5] = "Ryu19" * L"t \theta_B"
    time_labels[6] = "Ryu19" * L"t \theta_B p_\mathrm{inj}"
    time_labels[7] = "Ryu19" * L"t \theta_B p_\mathrm{inj} q_\alpha"


    mask_bad=trues(Nrows*Ncols)
    bad_colors = ["w" for _ = 1:Ncols]
    log_map = trues(Ncols)
    log_map[end] = false

    coutoff_arr = [-1.0, -1.0, -1.0, vmin_arr[4], vmin_arr[5], vmin_arr[6]]

    cutoffs = [coutoff_arr[i] for map_path ∈ map_paths, i ∈ 1:Ncols]

    propaganda_plot_columns(Nrows, Ncols, [""], im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name, upscale=0.9,
                            contour_color="k",
                            scale_label="1 Mpc";
                            time_labels, annotate_time,
                            map_arr, par_arr,
                            contour_files,
                            contour_arr, contour_par_arr, contour_levels,
                            contours,
                            log_map,
                            mask_bad,
                            bad_colors,
                            cutoffs
                            )

end


folders = [ "KR13d_thetaB", 
        "KR13t", 
        "KR13t_thetaB",
        "Ryu19t",
        "Ryu19t_thetaB",
        "Ryu19t_thetaB_pinj", 
        "Ryu19t_thetaB_pinj_q"
        ]

map_paths = ["data/cluster_mergers/$folder/maps/" for folder ∈ folders ]


plot_cr_total_compare(map_paths, "SR")

