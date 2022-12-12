# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO
using PyPlot, PyPlotUtility
using Printf
using DelimitedFiles

"""
    read_sinus_bounce_data_mean(path_in, N_output, Nbins)

Reads the test output of the adiabatic cycles test.

## Arguments
- `path_in`: Path to data file.
- `N_output`: Number of cycle step.
- `Nbins`: Number of CR bins used in the test. 
"""
function read_sinus_bounce_data_mean(path_in, N_output, Nbins)

    init_state = Float64.(readdlm(path_in * "CRp_adiabatic_Nbin-$(@sprintf("%03i", Nbins))_$(@sprintf("%03i", 0)).dat")[1:end-2, :])
    S_init = init_state[:, 2]
    E_init = init_state[:, 4]
    N_init = init_state[:, 5]

    state = Float64.(readdlm(path_in * "CRp_adiabatic_Nbin-$(@sprintf("%03i", Nbins))_$(@sprintf("%03i", N_output)).dat")[1:end-2, :])

    S = state[:, 2]
    E = state[:, 4]
    N = state[:, 5]

    S_error = abs(sum(S) - sum(S_init)) / sum(S_init)
    E_error = abs(sum(E) - sum(E_init)) / sum(E_init)
    N_error = abs(sum(N) - sum(N_init)) / sum(N_init)

    return S_error, E_error, N_error
end


"""
    plot_mean_errors(path_in, colors, Nbins, snap_range, plot_name)

Main plot routine
"""
function plot_mean_errors(path_in, colors, Nbins, snap_range, plot_name)

    N_end = snap_range[end]
    N_outputs = size(snap_range, 1)

    # allocate storage arrays
    S = Array{Float64,1}(undef, N_outputs)
    E = Array{Float64,1}(undef, N_outputs)
    N = Array{Float64,1}(undef, N_outputs)

    lw = 2.5

    fig = get_figure(1.0)
        plot_styling!()

        ax1 = gca()

        ax1.set_xlim([-1, N_end + 1])
        ax1.set_ylim([1.e-6, 3.e-2])

        ax1.set_yscale("log")

        axis_ticks_styling!(ax1)

        ylabel("Relative Error")
        xlabel("Number of Cycles  " * L"N_{\mathrm{cycle}}")

        # read data
        for i = 1:N_outputs
            S[i], E[i], N[i] = read_sinus_bounce_data_mean(path_in, snap_range[i], Nbins)
        end

        plot(snap_range, S, color = colors[1], lw = lw)
        plot(snap_range, E, color = colors[2], lw = lw)
        plot(snap_range, N, color = colors[3], lw = lw)

        plot(0.0, 0.0, color = colors[1], lw=lw, label = "Slope  " * L"q")
        plot(0.0, 0.0, color = colors[2], lw=lw, label = "Energy  " * L"E")
        plot(0.0, 0.0, color = colors[3], lw=lw, label = "Number  " * L"N")

        girichidis_solution = readdlm(path_in * "Girichidis2019_bounce.csv", ',')
        plot(girichidis_solution[:, 1], girichidis_solution[:, 2], color = "k", linestyle = ":", label = "Girichidis+20", lw = lw)

        legend(frameon = false, loc = "lower right")

    savefig(plot_name, bbox_inches = "tight")

    close(fig)
end

path_in = "data/tests/adiabatic/cycles/"
snap_range = collect(0:100)
Nbins = 12
colors = ["#000c8f", "#a700c9", "#c90000"]

plot_name = "Plots/Fig02.pdf"
plot_mean_errors(path_in, colors, Nbins, snap_range, plot_name)

