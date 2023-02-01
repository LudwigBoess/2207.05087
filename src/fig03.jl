# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility, PyCall
using ProgressMeter
using Printf
using DelimitedFiles
using SpectralCRsUtility

# needs to by imported by hand to make inset axis
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

# define constants
const me = 9.10953e-28    # electron mass
const σT = 6.65245e-25    # thompson cross section
const cL = 2.9979e10      # speed of light
const uIC = 4.00544e-13   # CMB energy density 
const synch_ic_const = 4σT / (3me^2 * cL^2)
const β = synch_ic_const * uIC  


"""
    get_f_0(p::T, f_i::T, p_i::T, q_0::T) where {T}

Returns initial value for f(p) at given p.
"""
function get_f_0(p::T, f_i::T, p_i::T, q_0::T) where {T}
    f_i * (p / p_i)^(-q_0)
end


"""
    maximum momentum at time t for CMB energy density uIC
"""
function get_p_cool(t::T) where {T}
    (3.0 / 4.0 * me^2 * cL^2) / (σT * uIC * t)
end



"""
    f_synch_ic(p::T, q_0::T, p_cool::T, t::T, f_i::T = 1.0, p_i::T = 1.0) where {T}

Analytic solution Ogridnik+21 f(p), Eq. 33
"""
function f_synch_ic(p::T, q_0::T, p_cool::T, t::T, f_i::T = 1.0, p_i::T = 1.0) where {T}
    if p < p_cool
        return get_f_0(p, f_i, p_i, q_0) * (1.0 - β * t * p)^(q_0 - 4.0)
    else
        return 0.0
    end
end


"""
    set_spectrum(f_0::T, q_0::T, Nbins::Int64) where {T}

Init powerlaw spectrum with given norm `f_0` and slope `q_0`.
"""
function set_spectrum(f_0::T, q_0::T, Nbins::Int64) where {T}

    bin_width = log10(1.e6 / 1.0) / Nbins

    bounds = Vector{T}(undef, Nbins + 1)
    for i = 1:Nbins+1
        bounds[i] = 10.0 .^ ((i - 1) * bin_width)
    end

    norm = get_f_0.(bounds[1:end-1], f_0, 1.0, q_0)

    return bounds, norm
end


"""
    read_analytic_solution(fi::String, f_0::T, q_0::T, Nbins::Int64) where {T}

Read the analytic solution for a given number of bins.
"""
function read_analytic_solution(filename::String, f_0::T, q_0::T, Nbins::Int64) where {T}

    h = read_header(filename)
    Δt = h.time

    bounds, norm = set_spectrum(f_0, q_0, Nbins)

    for i = 1:Nbins
        p_cool  = get_p_cool(Δt)
        norm[i] = f_synch_ic(bounds[i], q_0, p_cool, Δt, f_0)
    end

    return bounds[1:end-1], norm
end


"""
    multiplot_spectra(path_in, snap_range, id, q_0, plot_name)

Main plot routine.
"""
function multiplot_spectra(path_in, snap_range, id, q_0, plot_name)

    N = size(snap_range, 1)
    Nplots = length(path_in)

    GU = GadgetPhysical()


    fi = path_in[1] * "snap_$(@sprintf("%03i", snap_range[end]))"
    h = head_to_obj(fi)
    t_max = h.time * GU.t_Myr

    # analytic solution
    Nbins_fit = 10_000

    xmin = 0.7e3
    xmax = 2.e6
    ymin = [1.e-25, 1.e-37]
    ymax = [1.e-14, 1.e-18]

    xmin_ins = 7.e3
    xmax_ins = 1.5e4
    ymin_ins = [1.e-18, 1.e-28]
    ymax_ins = [1.e-15, 1.e-24]

    sm = plt.cm.ScalarMappable(cmap = PyPlot.cm.jet, norm = plt.Normalize(vmin = 0.0, vmax = t_max))
    sm.set_array([])


    x_pixels = 700

    fig = get_figure(2.0; x_pixels)

    plot_styling!(x_pixels)


    for plot_num = 1:Nplots

        subplot(1, 2, plot_num)
        ax = gca()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin[plot_num], ymax[plot_num]])
        ax.set_xscale("log")
        ax.set_yscale("log")
        axis_ticks_styling!(ax)

        xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_e c)^{-1}" * "]")
        ylabel("Distribution Function  " * L"f(\hat{p})")

        title("Initial Slope  " * L"q_0 \: = \: " * "- $(q_0[plot_num])")

        # define inset axis with zoomed window
        axins = inset_locator.inset_axes(ax, width = "40%", height = "40%", loc = "upper right")#[0.55, 0.55, 0.4, 0.4])
        axins.set_xlim([xmin_ins, xmax_ins])
        axins.set_ylim([ymin_ins[plot_num], ymax_ins[plot_num]])
        axins.set_xscale("log")
        axins.set_yscale("log")
        axis_ticks_styling!(axins)


        @showprogress "q = -$(q_0[plot_num]) " for i = 1:N

            fi = path_in[plot_num] * "snap_$(@sprintf("%03i", snap_range[i]))"
            h = head_to_obj(fi)

            CRp, CRe = getCRMomentumDistributionFromPartID(fi, id, pmin = 1.0, pmax = 1.e6)

            ax.plot(CRe.bound[1:end-1], CRe.norm, c = sm.to_rgba(h.time * GU.t_Myr))

            axins.plot(CRe.bound[1:end-1], CRe.norm, c = sm.to_rgba(h.time * GU.t_Myr))

        end

        # analytic solution
        fi = path_in[plot_num] * "snap_$(@sprintf("%03i", snap_range[end]))"
        CRp, CRe = getCRMomentumDistributionFromPartID(fi, id, pmin = 1.0, pmax = 1.e6)

        bounds, norm = read_analytic_solution(fi, CRe.norm[1], q_0[plot_num], Nbins_fit)
        ax.plot(bounds, norm, color = "w", lw = 4, alpha = 0.3)
        ax.plot(bounds, norm, color = "k", linestyle = "--", alpha = 0.8)

        axins.plot(bounds, norm, color = "w", lw = 4, alpha = 0.3)
        axins.plot(bounds, norm, color = "k", linestyle = "--", alpha = 0.8)

        axins.set_xticklabels([])
        axins.tick_params(which = "minor", label1On = false)
        axins.yaxis.set_ticklabels([])

        if plot_num == 1
            inset_locator.mark_inset(ax, axins, 2, 3, alpha = 0.8, zorder = 1000.0, linestyle = ":", lw = 1.5)
        else
            inset_locator.mark_inset(ax, axins, 2, 4, alpha = 0.8, zorder = 1000.0, linestyle = ":", lw = 1.5)
        end

    end

    cb = colorbar(sm, location = "right", anchor = (2.1, 0.5))
    cb.set_label("Time  " * L"t" * " [Myr]")
    cb.ax.tick_params(
        direction = "in",
        which = "major",
        size = 6, width = 1
    )

    savefig(plot_name, bbox_inches = "tight")
    close(fig)
end

path_0 = "data/tests/cooling/spectra/"
path_in = path_0 .* ["q_3/", "q_6/"]
plot_name = "Plots/Fig03.pdf"
snap_range = 0:100
q_0 = [3.0, 6.0]
multiplot_spectra(path_in, snap_range, 1, q_0, plot_name)