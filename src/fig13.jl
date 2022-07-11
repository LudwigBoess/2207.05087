# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using SpectralCRsUtility
using Printf
using ProgressMeter
using DelimitedFiles
using Statistics
using Base.Threads

"""
    get_spectral_data(ref_id, path_in, snap_range)

Obtains the spectrum of single particle over the simulation time.
"""
function get_spectral_data(ref_id, path_in, snap_range)

    N = size(snap_range, 1)

    CRp = Array{CRMomentumDistribution,1}(undef, N)
    CRe = Array{CRMomentumDistribution,1}(undef, N)
    t = Array{Float64,1}(undef, N)

    GU = GadgetPhysical()
    t_Gyr = GU.t_Myr / 1e3

    @showprogress for i = 1:N
        println("thread $(threadid()): i = $i")
        fi = path_in * "snap_$(@sprintf("%03i", snap_range[i]))"
        h = head_to_obj(fi)
        t[i] = h.time * t_Gyr
        CRp[i], CRe[i] = getCRMomentumDistributionFromPartID(fi, ref_id, pmin = 1.0e-1, pmax = 1.e5)
    end

    return t, CRp, CRe

end

"""
    plot_spectral_evolution(t, CRp, CRe, plot_name)

Main plot routine.
"""
function plot_spectral_evolution(t, CRp, CRe, plot_name)

    Nstart = 1
    Nsnaps = length(t)

    xmin = 1.0e-1
    xmax = 2.e5
    ymin = 1.e-40
    ymax = 1.e-10
    q0 = 4.0

    x_pixels = 600
    fig = get_figure(2.0; x_pixels)
    plot_styling!(x_pixels)

    gs = plt.GridSpec(1, 3, figure = fig, width_ratios = [1, 1, 0.05], wspace = 0.1)

    sm = plt.cm.ScalarMappable(cmap = PyPlot.cm.jet, norm = plt.Normalize(vmin = t[Nstart], vmax = t[Nsnaps]))
    sm.set_array([])

    subplot(get_gs(gs, 0, 0))

        ax = gca()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.set_xscale("log")
        ax.set_yscale("log")
        axis_ticks_styling!(ax)

        locmin = plt.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 20)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 12)
        ax.xaxis.set_major_locator(locmaj)

        xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_p c)^{-1}" * "]")
        ylabel("Distribution Function  " * L"f(p)")

        @showprogress "Protons..." for i = Nstart:Nsnaps
            plot(CRp[i].bound[1:end-1], CRp[i].norm, c = sm.to_rgba(t[i]))
        end

        text(1.e5, 1.e-13, "Protons", horizontalalignment = "right", fontsize = 20)

        plot([1.e3, 1.e4], [1.e-24, 1.e-24 * (1.e-1)^q0], color = "k", linestyle = "--", linewidth = 2)
        text(4.e3, 1.e-25, "q = -4")

        get_cr_energy_axis!(ax, "p")

    subplot(get_gs(gs, 0, 1))

        ax = gca()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.set_xscale("log")
        ax.set_yscale("log")

        axis_ticks_styling!(ax)

        locmin = plt.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 20)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 12)
        ax.xaxis.set_major_locator(locmaj)

        xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_e c)^{-1}" * "]")
        # ylabel("Distribution Function  " * L"f(p)")
        ax.set_yticklabels([])

        @showprogress "Electrons..." for i = Nstart:Nsnaps
            plot(CRe[i].bound[1:end-1], CRe[i].norm, c = sm.to_rgba(t[i]))
        end

        text(1.e5, 1.e-13, "Electrons", horizontalalignment = "right", fontsize = 20)

        plot([1.e3, 1.e4], [1.e-24, 1.e-24 * (1.e-1)^q0], color = "k", linestyle = "--", linewidth = 2, label = "q = -4")
        #legend(frameon = false, loc = "upper right")
        get_cr_energy_axis!(ax, "e")

    #cb = colorbar(sm, fraction = 0.046, pad = 0.04)
    subplot(get_gs(gs, 0, 2))
        cax = gca()
        cb = colorbar(sm, cax = cax, fraction = 0.046)
        cb.set_label("Time  t [Gyr]")
        cb.ax.tick_params(
            direction = "in",
            which = "major",
            size = 6, width = 1
        )


    savefig(plot_name, bbox_inches = "tight")
    close(fig)
end

ref_id = 4077968

path_out = "data/cluster_mergers/Ryu19t_thetaB_pinj_q/cr_spectra/"
output_file = path_out * "NR_spectra_$(ref_id)_"
t, CRp = read_cr_from_binary(output_file * "CRp.dat")
t, CRe = read_cr_from_binary(output_file * "CRe.dat")

plot_name = "Plots/Fig13.pdf"
plot_spectral_evolution(t, CRp, CRe, plot_name)


