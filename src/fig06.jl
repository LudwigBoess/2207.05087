# activate the environment
using Pkg; Pkg.activate(".")

using GadgetIO
using AnalyticMHDTestSolutions
using PyPlot, PyPlotUtility
using Statistics, StatsBase
using Distributions
using PyCall

inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

function get_snap_mean_values(fi, pos)

    h = head_to_obj(fi)

    info = read_info(fi, verbose = false)

    x = read_block(fi, "POS",
        info = info[getfield.(info, :block_name).=="POS"][1],
        parttype = 0)

    j = findall(pos .<= x[1, :] .< pos + 1.5)

    ρ = read_block(fi, "RHO",
        info = info[getfield.(info, :block_name).=="RHO"][1],
        parttype = 0)[j]

    u = read_block(fi, "U",
        info = info[getfield.(info, :block_name).=="U"][1],
        parttype = 0)[j]

    CRpP = read_block(fi, "CRpP",
        info = info[getfield.(info, :block_name).=="CRpP"][1],
        parttype = 0)[j]

    CReP = read_block(fi, "CReP",
        info = info[getfield.(info, :block_name).=="CReP"][1],
        parttype = 0)[j]

    u_mean = mean(u)

    P_cr = mean(CRpP .+ CReP)

    ρ_mean = mean(ρ)

    γ_cr = 4.0 / 3.0

    U_cr = P_cr / ((γ_cr - 1.0) * ρ_mean)

    return U_cr, u_mean
end

function get_b_angle_eff_sim(theta::Int64)

    fi = "data/tests/injection/eta_B/theta$theta/snap_000"
    CRpE, u = get_snap_mean_values(fi, 84.5)
    U_E_ratio = CRpE / (u + CRpE)

    return U_E_ratio
end

function get_ideal_etaB(theta::Float64)

    theta /= (180.0 / π)

    THETA_CRIT = π / 4.0
    DELTA_THETA = π / 18.0
    eta_B4 = 0.5 * (tanh((THETA_CRIT - theta) / DELTA_THETA) + 1.0)

    return eta_B4
end


function get_snap_pressure(fi)

    h = head_to_obj(fi)

    info = read_info(fi, verbose = false)

    x = read_block(fi, "POS",
        info = info[getfield.(info, :block_name).=="POS"][1],
        parttype = 0)

    # sort along x direction
    j = sortperm(x[1, :])

    x = x[1, j]

    ρ = read_block(fi, "RHO",
        info = info[getfield.(info, :block_name).=="RHO"][1],
        parttype = 0)[j]

    u = read_block(fi, "U",
        info = info[getfield.(info, :block_name).=="U"][1],
        parttype = 0)[j]


    CRpP = read_block(fi, "CRpP",
        info = info[getfield.(info, :block_name).=="CRpP"][1],
        parttype = 0)[j]

    CReP = read_block(fi, "CReP",
        info = info[getfield.(info, :block_name).=="CReP"][1],
        parttype = 0)[j]

    P_th = @. (5.0 / 3.0 - 1) * ρ * u

    return x, CRpP + CReP, P_th + CRpP + CReP
end

function plot_b_angle_eff_and_postshock_region(path, plot_name)


    #data = readdlm(path, comments=true)
    θ = collect(0.0:0.01:90.0)
    η = get_ideal_etaB.(θ)

    angle_arr = [0, 15, 30, 45, 60, 75, 90]
    eta_B = get_b_angle_eff_sim.(angle_arr)
    eta_ref = eta_B[1] 
    eta_B ./= eta_ref

    eta_ideal = get_ideal_etaB.(Float64.(angle_arr))

    eta_error = abs.(eta_ideal .- eta_B) ./ eta_ideal

    #sim_color = "#7a059e"
    sim_color = "#b5160b"

    thetaB = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]

    xs_first_guess = [4.7, 4.5,         # 0
                        4.3, 4.3,    # 30 / 45
                        3.8,         # 60
                        3.8, 3.8]   # 75 / 90

    xmin = 50.0
    xmax = 100.0
    x_step = 0.01
    x_ideal = collect(xmin:x_step:xmax)
    P_left = 63.499
    P_right = 0.1


    color_map = "jet"
    sm = PyPlot.cm.ScalarMappable(cmap = color_map, norm = PyPlot.Normalize01)

    xmin = -5.0
    xmax = 95.0

    linestyles = ["-.", "-"]

    x_pixels = 600

    fig = get_figure(2.2; x_pixels)

    plot_styling!(x_pixels)

    subplot(1, 2, 1)
    ax = gca()
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([1.e-4, 2.0])
    ax.set_yscale("log")


    locmin = plt.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 20)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 20)
    ax.yaxis.set_major_locator(locmaj)

    locmaj_x = matplotlib.ticker.MultipleLocator(15.0)
    ax.xaxis.set_major_locator(locmaj_x)

    for tick in ax.xaxis.get_major_ticks()
        tick.label.set_fontsize(12)
    end
    for tick in ax.yaxis.get_major_ticks()
        tick.label.set_fontsize(12)
    end

    #xticks(x)
    axis_ticks_styling!(ax)
    ax.minorticks_on()
    #ax.xaxis.set_minor_locator(plt.MultipleLocator(1))


    xlabel("Angle  " * L"θ_B")
    ylabel("Angle dep. Efficiency  " * L"η(θ_B)")

    #plot(data[:,1], data[:,2], color="#4d62db", label="Ideal")
    plot(θ, η, #color = "#4d62db", 
        color="#000c8f",
        label = "Model", lw = 3.0, zorder = 0.0,
        alpha=0.8)

    scatter(angle_arr, eta_B, label = "Simulation", marker = "x", s=75.0, color = sim_color, lw = 2.0)

    legend(frameon = false, loc = "upper right")

    # inset plot
    large_inbox_size = 0.45
    left, bottom, width, height = [0.26, 0.21, 0.7 * large_inbox_size, 0.7 * large_inbox_size]
    ax_small2 = inset_locator.inset_axes(ax, width = "40%", height = "40%", loc = "lower left", borderpad = 3.5)# bbox_to_anchor = (left, bottom, width, height))# loc = "lower left")

    axis_ticks_styling!(ax_small2, size_minor_ticks = 2, tick_label_size = 10)
    ax_small2.set_xlim([-5.0, 95])
    ax_small2.set_ylim([1.e-5, 1.0])
    ax_small2.set_yscale("log")
    #ax_small2.tick_params(reset = true, direction = "in", which = "both")
    ax_small2.minorticks_on()

    # locmaj_x2 = matplotlib.ticker.IndexLocator(base=15.0, offset=5.0)
    ax_small2.xaxis.set_major_locator(locmaj_x)
    ax_small2.yaxis.set_major_locator(locmaj)
    ax_small2.yaxis.set_minor_locator(locmin)
    ax_small2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    #ax_small2.set_ylabel(L"L_1" * " Error")
    ax_small2.set_ylabel("Rel. Error")
    ax_small2.set_xlabel("Angle  " * L"θ_B")

    ax_small2.scatter(angle_arr, eta_error, marker = "x", color = sim_color, lw = 2.0)

    subplot(1, 2, 2)
    ax = gca()

    ax = gca()
    ax.set_xlim([80.0, 90.0])
    ax.set_ylim([1.e-5, 20.0])
    ax.set_yscale("log")

    axis_ticks_styling!(ax)

    ax.minorticks_on()

    xlabel("Position  " * L"x")
    ylabel("Downstr. CR Pressure  " * L"P_{cr,2}")

    axvline(84.5, color = "grey", linestyle = ":", alpha = 0.75)
    axvline(86.0, color = "grey", linestyle = ":", alpha = 0.75)

    for i = 1:size(thetaB, 1)

        println("Getting data for θ_B = $(thetaB[i])°")

        fi = path * "theta$(Int(thetaB[i]))/snap_000"
        x, P_cr, P_tot = get_snap_pressure(fi)

        par = RiemannParameters(Pl = P_left, Pr = P_right,
            thetaB = thetaB[i],
            xs_first_guess = xs_first_guess[i],
            t = 1.5, dsa_model = 4, Pe_ratio = 0.0)

        sol = AnalyticMHDTestSolutions.solve(x_ideal, par)

        plot(x, P_cr, color = sm.to_rgba(thetaB[i] / thetaB[end]))

        plot(sol.x, sol.P_cr_p, "--", color = sm.to_rgba(thetaB[i] / thetaB[end]),
            label = L"θ_B = " * "$(thetaB[i])°")


    end

    legend(frameon = false, bbox_to_anchor = (1, 0.75))#, loc="center")

    subplots_adjust(wspace = 0.3)
    savefig(plot_name, bbox_inches = "tight")

    close(fig)

end

path_in = "data/tests/injection/eta_B/"
plot_name = "Plots/Fig06.pdf"
plot_b_angle_eff_and_postshock_region(path_in, plot_name)
