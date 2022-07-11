# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using SPHtoGrid
using PyPlot, PyPlotUtility
using Printf

# path to the fits files
map_path_1 = "data/cluster_mergers/KR13d_thetaB/"
map_path_2 = "data/cluster_mergers/KR13t_thetaB/"
sim_path = map_path_1

# unit conversion struct
GU = GadgetPhysical()

# image definition
Nrows = 2 
Ncols = 5
Nsnap = 23

# fits files to read
files = [map_path_1 * "maps/full_box_rho.000.a.z.pix",
        map_path_1 * "maps/full_box_rho.$(@sprintf("%03i", Nsnap)).a.z.pix",
        map_path_1 * "maps/full_box_T.000.a.z.pix",
        map_path_1 * "maps/full_box_T.$(@sprintf("%03i", Nsnap)).a.z.pix",
        map_path_1 * "maps/full_box_Xray.000.a.z.pix",
        map_path_1 * "maps/full_box_Xray.$(@sprintf("%03i", Nsnap)).a.z.pix",
        map_path_1 * "maps/full_box_B.000.a.z.pix",
        map_path_1 * "maps/full_box_B.$(@sprintf("%03i", Nsnap)).a.z.pix",
        map_path_2 * "maps/full_box_B.000.a.z.pix",
        map_path_2 * "maps/full_box_B.$(@sprintf("%03i", Nsnap)).a.z.pix"]


vmin_arr = [1.e-4, 1.e7, 1.e-8, 1.e-2, 1.e-2]
vmax_arr = [1.e-1, 5.e8, 1.e-4, 1.e1, 1.e1]

im_cmap = ["viridis", "hot", "plasma", "magma", "magma"]

cb_labels = [ "Surface Density  " * L"\Sigma_g" * " [ g cm" * L"^{-2}" * "]",
              "Temperature  T [K]",
              "Xray Surf. Brightn.  " * L"S_{X,bol}" * "[ 10" * L"^{44}" * " erg s" * L"^{-1}" * "]",
              "Magnetic Field  B [" * L"μ" * "G]",
              "Magnetic Field  B [" * L"μ" * "G]" ]

snap_base = sim_path * "snap_$(@sprintf("%03i", Nsnap))"


time_labels = [L"t \: = " * "$(@sprintf("%0.2f", 0.0)) Gyrs",
                L"t \: = " * "$(@sprintf("%0.2f", 2.25)) Gyrs",
                "", "", "", "", "", "", "", ""]
                
annotate_time = falses(Ncols * Nrows)
annotate_time[1:2] .= true

annotate_text = falses(Ncols * Nrows)
annotate_text[7] = true
annotate_text[9] = true

text_labels = ["", "", "", "", "", "", "B Dipole", "", "B Turb", ""]

plot_name = "Plots/Fig09.pdf"

propaganda_plot_columns(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name;
            time_labels = time_labels,
            upscale = 2.5,
            scale_label = L"5" * " Mpc",
            scale_kpc = 5000.0,
            read_mode=2, 
            annotate_time, 
            annotate_text,
            text_labels
        )
