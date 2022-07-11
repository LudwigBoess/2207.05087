# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO
using PyPlot, PyPlotUtility
using SpectralCRsUtility
using Printf
using StatsBase
using DelimitedFiles
using ProgressMeter

"""
    get_min_max_id(filename)

Reads the IDs of the particles with minimum and maximum density at t=0.
"""
function get_min_max_id(filename)

    id = read_block(filename, "ID", parttype=0)
    ρ  = read_block(filename, "RHO", parttype=0)

    k = findmax(ρ)[2][1]
    max_id = id[k]

    k = findmin(ρ)[2][1]
    min_id = id[k]

    return max_id, min_id
end

"""
    read_density(filename)

Reads density and position and returns both sorted by position.
"""
function read_density(filename)

    info_pos = InfoLine("POS", Float32, 3, [1, 0, 0, 0, 0, 0])
    info_rho = InfoLine("RHO", Float32, 1, [1, 0, 0, 0, 0, 0])

    x = read_block(filename, "POS", info=info_pos, parttype=0)
    j = sortperm(x[1,:])
    x = x[1,j]

    ρ = read_block(filename, "RHO", info=info_rho, parttype=0)[j]

    return x, ρ
end

"""
    density_of_id(filename, find_id)

Reads the density of the requested ID to track the tracer particles.
"""
function density_of_id(filename, find_id)

    id = read_block(filename, "ID", parttype=0)

    k  = findfirst(id .== find_id)[1]

    ρ  = read_block(filename, "RHO", parttype=0)[k]
    x  = read_block(filename, "POS", parttype=0)[1,k]

    return x, ρ
end


"""
    plot_sinus_spectra(path_in, snap_range, max_id, min_id, plot_name )

Main plot routine
"""
function plot_sinus_spectra(path_in, snap_range, max_id, min_id, plot_name )

    """
        Plot
    """
    xmin = -0.1
    xmax = 1.1

    markersize = 75
    
    sm = plt.cm.ScalarMappable(cmap=PyPlot.cm.jet, norm=plt.Normalize(vmin=0.0, vmax=20.0))
    sm.set_array([])

    upscale = 0.7
    #fig = figure(figsize=(22*upscale, 10*upscale), constrained_layout=true)
    fig = get_figure(2.5)
    plot_styling!()
    gs = plt.GridSpec(2, 3, hspace=0.0, wspace=0.3, width_ratios=[1, 1, 0.05], figure = fig)
    #gs = fig.add_gridspec(2, 3, hspace=0.0, wspace=0.3, width_ratios=[1, 1, 0.05])

    subplot(get_gs(gs, 0:2, 0))

        ax = gca()
        axis_ticks_styling!(ax)

        ax.set_xlim([xmin,xmax])

        xlabel("Position  " * L"x" * " [arb. units]")
        ylabel("Density  " * L"ρ" * " [arb. units]")

        @showprogress "Sinus" for i ∈ snap_range
            filename = path_in * "snap_$(@sprintf("%03i", i))"
            x, ρ = read_density(filename)
            h    = GadgetIO.read_header(filename)

            # all particles
            scatter(x, ρ, c=sm.to_rgba(h.time), s=pixel_size(fig), 
                    rasterized=true)

            # upper tracer particle
            x_id, ρ_id = density_of_id(filename, max_id)
            scatter(x_id, ρ_id, c=sm.to_rgba(h.time), 
                    marker="o", s=markersize, 
                    linewidth=1.0, 
                    edgecolor="k")

            # lower tracer particle
            x_id, ρ_id = density_of_id(filename, min_id)
            scatter(x_id, ρ_id, c=sm.to_rgba(h.time), 
                    marker="X", s=markersize, 
                    linewidth=1.0, 
                    edgecolors="k")
        end


    # upper spectrum
    subplot(get_gs(gs, 0, 1))

        ax = gca()
        axis_ticks_styling!(ax)

        ax.set_xlim([0.8,2.e6])
        ax.set_ylim([-14.7,-13.7])
        ax.set_xscale("log")
        ax.minorticks_on()
        locmin = plt.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 20)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 12)
        ax.xaxis.set_major_locator(locmaj)
        ax.set_xticklabels([])


        ylabel("log" * L"_{10}" * "( " * L"p^{-q_0} \: f(p)" * " )")

        @showprogress "particle O " for i ∈ snap_range
            filename = path_in * "snap_$(@sprintf("%03i", i))"
            h    = GadgetIO.read_header(filename)
            cr = getCRMomentumDistributionFromPartID( filename, max_id, mode=3, electrons=false )
            plot(cr.bound[1:end-1], log10.(cr.bound[1:end-1].^4.5 .* cr.norm),
                 c=sm.to_rgba(h.time))

        end

        # attach marker
        scatter(0.0, 0.0, c="k", marker="o", s=2markersize, label=" ")
        legend(frameon=false, loc="upper right")

    subplot(get_gs(gs, 1, 1))

        ax = gca()
        axis_ticks_styling!(ax)

        ax.set_xlim([0.8,2.e6])
        ax.set_ylim([-14.7,-13.7])

        ax.set_xscale("log")

        ax.minorticks_on()
        locmin = plt.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 20)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 12)
        ax.xaxis.set_major_locator(locmaj)

        ylabel("log" * L"_{10}" * "( " * L"p^{-q_0} \: f(p)" * " )")
        xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_p c)^{-1}" * "]")

        @showprogress "particle X " for i ∈ snap_range
            filename = path_in * "snap_$(@sprintf("%03i", i))"
            h    = GadgetIO.read_header(filename)
            cr = getCRMomentumDistributionFromPartID( filename, min_id, mode=3, electrons=false )

            plot(cr.bound[1:end-1], log10.(cr.bound[1:end-1].^4.5 .* cr.norm),
                 c=sm.to_rgba(h.time))

        end

        # attach marker
        scatter(0.0, 0.0, c="k", marker="X", s=2markersize, label=" ")
        legend(frameon=false, loc="upper right")
    
    subplot(get_gs(gs, 0:2, 2))
        cax = gca()
        cb = colorbar(sm, cax=cax, fraction=0.046)

        cb.set_label("Time " * L"t" * " [arb. units]")
        cb.ax.tick_params(
                            direction="in",
                            which="major",
                            size=6, width=1
                        )

        cb.ax.tick_params(
                            direction="in",
                            which="minor",
                            size=3, width=1
                        )

        cax.minorticks_on()


    savefig(plot_name, bbox_inches="tight", dpi=400)

    close(fig)

end


path_in = "data/tests/adiabatic/sinus/"
snap_range = collect(0:8)
Nbins = 48

max_id, min_id = get_min_max_id(path_in * "snap_000")

plot_name = "Plots/FigC1.pdf"

plot_sinus_spectra(path_in, snap_range, max_id, min_id, plot_name )
