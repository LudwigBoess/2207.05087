# activate the environment
using Pkg; Pkg.activate(".")

using GadgetIO
using Printf
using PyPlot, PyPlotUtility
using Interpolations


"""
    get_q(q0, p)

Analytic solution for slope at momenum `p` for a spectrum with initial slope `q0` after the cooling time for critical momentum `pn` is reached.
"""
function get_q(q0, p, pn = 1.e4)
    return q0 + (q0 - 4.0)* ( (p/pn)/( 1.0 - p/pn))
end

"""
    read_Kardashev_sim(fi)

Reads the slope of the first SPH particle
"""
function read_Kardashev_sim(fi)
    read_block(fi, "CReS", parttype=0)[:,1]
end

"""
    get_p_sim(pmin, pmax, Nbins)

Construct momentum boundaries.
"""
function get_p_sim(pmin, pmax, Nbins)

    bin_width = log10(pmax/pmin)/Nbins
    p_sim = zeros(Nbins+1)

    for i = 1:Nbins+1
        p_sim[i] = pmin * 10.0^(bin_width*(i-1))
    end

    return p_sim[1:end-1]
end

"""
    function get_rel_error(q_id, q)

Computes the relative error between analytic slope and simulated slope.
"""
@inline function get_rel_error(q_id, q)
    abs(q_id - q)/q_id
end

"""
    plot_Kardashev_sim_error(snap, q0, plot_name, Nbins = 1000)

Main plotting routine.
"""
function plot_Kardashev_sim_error(snap, q0, plot_name, Nbins = 1000)

    # momentum limits used in the simulation
    pmin = 1.0
    pmax = 1.e6
    bin_width = log10(pmax / pmin) / Nbins

    # allocate arrays for momentum and slope of analytic solution
    q = zeros(Nbins)
    p = zeros(Nbins)

    # get analytic solution
    for i = 1:Nbins
        p[i] = pmin * 10.0^(bin_width * (i - 1))
        q[i] = get_q(q0, p[i])
    end

    q_interp = LinearInterpolation(p, q)

    bin_arr = [12,
        24,
        48,
        96,
        192,
        384,
        768
    ]

    sm = plt.cm.ScalarMappable(cmap = PyPlot.cm.jet_r, norm = plt.Normalize(vmin = 1, vmax = size(bin_arr, 1)))
    sm.set_array([])


    fig = get_figure(1.0, x_pixels=700)
    plot_styling!(700)

    gs = plt.GridSpec(2, 1, figure = fig, height_ratios = [1, 0.5])

    # slope plot
    subplot(get_gs(gs, 0, 0))

        ax = gca()
        ax.set_xlim([1.0e2, 2.e4])
        ax.set_ylim([-20.0, -4.0])
        ax.set_xscale("log")
        axis_ticks_styling!(ax)

        ylabel("Slope  " * L"q" * " of bin")

        ax.set_xticklabels([])

        for i = 1:size(bin_arr, 1)
            @info "Slope $(bin_arr[i]) bins"
            p_sim = get_p_sim(1.0, 1.e6, bin_arr[i])
            s_sim = read_Kardashev_sim(snap[i])

            plot(p_sim, -s_sim, c = sm.to_rgba(i), label = L"N_{bins} =" * " $(@sprintf("%3i", bin_arr[i]/6))", lw=2.0)
        end

        plot(p, -q, label = "Kardashev '62", color = "k", alpha = 0.8, linestyle = ":", lw = 2.0)

        legend(frameon = false)

    # rel. error
    subplot(get_gs(gs, 1, 0))

        ax = gca()
        ax.set_xlim([1.0e2, 2.e4])
        ax.set_ylim([3.e-5, 5.0])
        ax.set_xscale("log")
        ax.set_yscale("log")

        axis_ticks_styling!(ax)

        ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base = 10, numticks = 6))

        ylabel("Rel. Error")
        xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_e c)^{-1}" * "]")

        for i = 1:size(bin_arr, 1)
            @info "Error $(bin_arr[i]) bins"
            p_sim = get_p_sim(1.0, 1.e6, bin_arr[i])
            s_sim = read_Kardashev_sim(snap[i])

            q_id = q_interp.(p_sim)

            l1_sim = get_rel_error.(q_id, s_sim)

            sel = findall(1.e-6 .< l1_sim .< 10.0)

            plot(p_sim[sel], l1_sim[sel], c = sm.to_rgba(i), label = "Nbin = $(@sprintf("%3i", bin_arr[i]))", lw = 2.0)#, marker="x")
        end

    #legend(frameon=false, fontsize=legend_font_size)

    subplots_adjust(hspace = 0)
    savefig(plot_name,
        bbox_inches = "tight")

    close(fig)
end


snap_name_func(i) = "data/tests/cooling/convergence/$(i)bins/snap_000"

bin_arr = [ 12,
            24,
            48,
            96,
           192,
           384,
           768
          ]

snap = snap_name_func.(bin_arr)
plot_name = "Plots/Fig04.pdf"
q0 = 6.0

plot_Kardashev_sim_error(snap, q0, plot_name, 10000)