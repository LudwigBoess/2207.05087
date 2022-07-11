# activate the environment
using Pkg; Pkg.activate(".")

using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using AnalyticMHDTestSolutions
using StatsBase
using Statistics
using LaTeXStrings
using Printf
using PyCall
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")


"""
    get_data_for_cluster_shocktube(fi)

Read the simulation data into a `Dict`.
"""
function get_data_for_cluster_shocktube(fi)

    # read positions
    x = read_block(fi, "POS", parttype=0)

    # sort along x direction
    j = sortperm(x[1,:])
    x = x[1,j]

    γ_th = 5.0/3.0
    γ_cr = 4.0/3.0

    hsml = read_block(fi, "HSML", parttype=0)[j]

    ρ = read_block(fi, "RHO", parttype=0)[j]

    u = read_block(fi, "U", parttype=0)[j]

    P_th = ( γ_th - 1.0 ) .* ρ .* u

    MACH = read_block(fi, "MACH", parttype=0)[j]

    CRpP = read_block(fi, "CRpP", parttype=0)[j]

    CReP = read_block(fi, "CReP", parttype=0)[j]

    M_A  = read_block(fi, "MALF", parttype=0)[j]  

    CRpE = CRpP ./ ( ( γ_cr - 1.0 ) .* ρ )
    CReE = CReP ./ ( ( γ_cr - 1.0 ) .* ρ )

    mach_selected = findall( 160.0 .<= x .< 180.0 )
    mach_max_pos = findmax(MACH[mach_selected])[2]
    mach_max = MACH[mach_selected[mach_max_pos]]
    mach_max_x = x[mach_selected[mach_max_pos]]
    mean_hsml = mean(hsml)

    GU = GadgetPhysical()

    return Dict("x" => x, 
                "rho" => ρ .* GU.rho_cgs * 1.e28,
                "P_th" => P_th  .* GU.P_th_cgs .* 1.e11,
                "P_crp" => CRpP .* GU.P_CR_cgs .* 1.e11,
                "P_cre" => CReP .* GU.P_CR_cgs .* 1.e11,
                "P_tot" => (P_th .* GU.P_th_cgs .+ 
                            CRpP .* GU.P_CR_cgs .+ 
                            CReP .* GU.P_CR_cgs )  .* 1.e11,
                "T" => u .* GU.T_K .* 1.e-8,
                "T_crp" => CRpE .* GU.T_K .* 1.e-8,
                "T_cre" => CReE .* GU.T_K .* 1.e-8,
                "M_s" => MACH, 
                "M_A" => M_A, 
                "mach_max" => mach_max, 
                "mach_max_x" => mach_max_x, 
                "mean_hsml" => mean_hsml
                )
end

"""
    get_ideal_solution_cluster_shocktube(t)

Calculate analytic solution at time `t`.
"""
function get_ideal_solution_cluster_shocktube(t)

    # read ICs
    fi = "data/tests/injection/NLDSA/snap_IC"
    rho = read_block(fi, "RHO", parttype = 0)
    u = read_block(fi, "U", parttype = 0)
    rhol = Float64(maximum(rho))
    rhor = Float64(minimum(rho))
    Ul = Float64(maximum(u))
    Ur = Float64(minimum(u))
    x_contact = 140.0

    # define model parameters
    dsa_model = 2
    xs_first_guess = 3.3
    Mach = 5.541278364096982
    par = RiemannParameters(; rhol, rhor, Ul, Ur, t, x_contact, dsa_model, xs_first_guess, Mach)

    x_ideal = collect(LinRange(0.0, 280.0, 10_000))
    sol = AnalyticMHDTestSolutions.solve(x_ideal, par)

    # unit conversion 
    GU = GadgetPhysical()

    sol.rho .*= GU.rho_cgs .* 1.e28
    sol.P_th .*= GU.P_th_cgs .* 1.e11
    sol.P_cr_p .*= GU.P_CR_cgs .* 1.e11
    sol.P_cr_e .*= GU.P_CR_cgs .* 1.e11
    sol.P_tot = @. sol.P_th + sol.P_cr_p + sol.P_cr_e

    sol.U_th .*= GU.T_K .* 1.e-8
    sol.E_cr_p .*= GU.T_K .* 1.e-8
    sol.E_cr_e .*= GU.T_K .* 1.e-8


    return sol, par
end


"""
    plot_cluster_cr_shocktube_2x2_plots(fi, plot_name)

Main plot routine.
"""
function plot_cluster_cr_shocktube_2x2_plots(fi, plot_name)


    GU = GadgetPhysical()

    h = head_to_obj(fi)
    time_end = h.time

    @info "reading simulation data..."
    data = get_data_for_cluster_shocktube(fi)
    @info "done!"

    @info "reading analytic solution..."
    sol, par = get_ideal_solution_cluster_shocktube(time_end)
    @info "done!"

    # limits and calculated quantities

    mach_ideal = sol.Mach
    # alfven velocity upstream is hard-coded for simplicity
    M_A_ideal = sol.vs / 285.43817

    mach_lim = max(mach_ideal, M_A_ideal)

    # config
    mach_lim += 0.1 * mach_lim

    P_ymax = 4.5 * GU.P_th_cgs .* 1.e11
    P_ymin = -0.2 * GU.P_th_cgs .* 1.e11

    # x limits for plot
    xmin = 90.0
    xmax = 200.0

    @info "Plotting..."
    fig = get_figure(x_pixels = 1000)
    plot_styling!()

    gs = plt.GridSpec(4, 2, figure = fig)

    label_pad = 0
    labelx = -0.1

    """
        Density
    """
    @info "  Density"
    subplot(get_gs(gs, 0:2, 0))

        ax = gca()
        ax.set_xlim([xmin, xmax])

        axis_ticks_styling!(ax)

        ylabel("Density  " * L"\rho" * " [ " * L"10^{-28}" * "g cm" * L"^{-3}" * " ]", labelpad = label_pad)

        ax.yaxis.set_label_coords(labelx, 0.5)

        # Turn off tick labels
        ax.set_xticklabels([])

        plot(data["x"], data["rho"], color = "darkblue")
        plot(sol.x, sol.rho, linestyle = "--", color = "red", alpha = 0.8)


    """
        Pressure
    """
    @info "  Pressure"
    subplot(get_gs(gs, 0:2, 1))
        ax = gca()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([P_ymin, P_ymax])

        axis_ticks_styling!(ax)
        ax.set_xticklabels([])

        ylabel("Pressure  " * L"P" * " [ " * L"10^{-11}" * " erg cm" * L"^{-3}" * " ]", labelpad = label_pad)
        ax.yaxis.set_label_coords(labelx, 0.5)

        # Turn off tick labels
        ax.set_xticklabels([])

        
        plot(data["x"], data["P_crp"], color = "darkblue", label = L"P_{CR,p}")
        plot(sol.x, sol.P_cr_p, "--", color = "darkblue")
        plot(data["x"], data["P_th"], color = "red", label = L"P_{th}")
        plot(sol.x, sol.P_th, "--", color = "red")
        plot(data["x"], data["P_tot"], color = "black", label = L"P_{tot}")
        plot(sol.x, sol.P_tot, "--", color = "black")

        legend(frameon = false, loc = "lower left")

        # # inset plot
        ax_small = inset_locator.inset_axes(ax, width = "50%", height = "50%")

        ax_small.set_xlim([160.0, 178.0])
        ax_small.set_ylim([1.e-4, 1.0])
        ax_small.set_yscale("log")

        axis_ticks_styling!(ax_small)

        ax_small.plot(data["x"], data["P_crp"], color = "darkblue", label = L"P_{CR,p}")
        ax_small.plot(sol.x, sol.P_cr_p, "--", color = "darkblue")
        ax_small.plot(data["x"], data["P_th"], color = "red", label = L"P_{th}")
        ax_small.plot(data["x"], data["P_tot"], color = "black", label = L"P_{tot}")
        ax_small.plot(sol.x, sol.P_th, "--", color = "red")
        ax_small.plot(sol.x, sol.P_tot, "--", color = "black")

        ax_small.set_xticklabels([])
        inset_locator.mark_inset(ax, ax_small, 3, 4, alpha = 0.8, zorder = 1000.0, linestyle = ":")

    """
        Temperature
    """
    @info "  Temperature"
    subplot(get_gs(gs, 2:4, 0))
        ax = gca()
        ax.set_xlim([xmin, xmax])
        axis_ticks_styling!(ax)

        xlabel("Position  " * L"x")
        ylabel("Temperature  " * L"T" * " [ " * L"10^8" * " K ]", labelpad = label_pad)
        ax.yaxis.set_label_coords(labelx, 0.5)

        plot(data["x"], data["T"], color="darkblue", label = L"T")
        plot(data["x"], data["T_crp"], color = "red", label = L"T_{CR,p}")
        plot(data["x"], data["T_cre"], color = "orange", label = L"T_{CR,e}")

        plot(sol.x, sol.U_th, linestyle = "--", color = "grey", alpha = 0.8)
        plot(sol.x, sol.E_cr_e, linestyle = "--", color = "grey", alpha = 0.8)
        plot(sol.x, sol.E_cr_p, linestyle = "--", color = "grey", alpha = 0.8)


        # inset plot
        ax_small = inset_locator.inset_axes(ax, width = "25%", height = "25%")#, loc=10 )
    
        γ_th = 5 / 3
        γ_cr = 4 / 3
        U4_ideal = sol.P4_th / ((γ_th - 1.0) * sol.rho4)
        Ecr4_ideal = sol.P4_cr / ((γ_cr - 1.0) * sol.rho4)
        η_ideal = Ecr4_ideal / U4_ideal
            
        ax_small.set_xlim([160, 176])
        ax_small.set_ylim([1.0e-3, 0.1])
        ax_small.set_yscale("log")

        axis_ticks_styling!(ax_small)
    
        U_E_ratio = (data["T_crp"] .+ data["T_cre"]) ./ data["T"] #( u .+ CRpE )
    
        U_E_ratio[U_E_ratio.<1.e-4] .= 1.e-4
    
        ax_small.plot(data["x"], U_E_ratio, color = "darkblue")
        ax_small.axhline(η_ideal, linestyle = "--", color = "grey", alpha = 0.8, label = "ideal")
        xlabel(L"x")
        ylabel(L"η", fontsize = 15, labelpad = label_pad)
        legend(frameon = false, loc = "upper left", fontsize = 10)
    
        # zoom on post-shock
        ax_small2 = inset_locator.inset_axes(ax, width = "25%", height = "25%", loc = "lower center")
    
        ax_small2.set_xlim([160, 176])
        ax_small2.set_xticks([])
        ax_small2.set_ylim([1.0e-3, 2.0])
        ax_small2.set_yscale("log")

        ax_small2.plot(sol.x, sol.U_th, linestyle = "--", color = "red", alpha = 0.8)
        ax_small2.plot(sol.x, sol.E_cr_e, linestyle = "--", color = "orange", alpha = 0.8)
        ax_small2.plot(sol.x, sol.E_cr_p, linestyle = "--", color = "darkblue", alpha = 0.8)
        plot(data["x"], data["T"], color="darkblue", label = L"T")
        plot(data["x"], data["T_crp"], color = "red", label = L"T_{CR,p}")
        plot(data["x"], data["T_cre"], color = "orange", label = L"T_{CR,e}")


        inset_locator.mark_inset(ax, ax_small2, 2, 4, alpha = 0.8, zorder = 1000.0, linestyle = ":")

    # Alfven Mach number
    @info "  Alvfen Mach number"
    subplot(get_gs(gs, 2, 1))

        ax = gca()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([0.0, 13.0])

        axis_ticks_styling!(ax)
        ax.set_xticklabels([])

        #xlabel(L"x", fontsize=axis_label_font_size)
        ylabel("Alfvén Machn.  " * L"\mathcal{M}_A", labelpad = label_pad)
        ax.yaxis.set_label_coords(labelx, 0.5)

        plot(data["x"], data["M_A"], color = "darkblue", label = L"M_A")
        axhline(M_A_ideal, linestyle = "--", color = "grey")

        # M_Alfven
        ax_small = inset_locator.inset_axes(ax, width = "25%", height = "50%", loc = 10)

        ax_small.set_xlim([172.0, 174.5])
        ax_small.set_ylim([11.0, 12.5])

        axis_ticks_styling!(ax_small)
        ax_small.axhline(M_A_ideal, linestyle = "--", color = "grey")

        ax_small.plot(data["x"], data["M_A"], color = "darkblue")

        ax_small.set_xticklabels([])
        ax_small.set_yticklabels([])
        inset_locator.mark_inset(ax, ax_small, 2, 4, alpha = 0.8, zorder = 1000.0, linestyle = ":")


    # Sonic Mach number
    @info "  Sonic Mach number"
    subplot(get_gs(gs, 3, 1))

        ax = gca()
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([0.0, 7.0])

        axis_ticks_styling!(ax)

        xlabel("Position  " * L"x")
        ylabel("Sonic Machn.  " * L"\mathcal{M}_s", labelpad = label_pad)
        ax.yaxis.set_label_coords(labelx, 0.5)

        plot(data["x"], data["M_s"], color = "darkblue", label = L"M_s")
        axhline(5.532955671044105, linestyle = "--", color = "grey")

        # inset plot
        ax_small = inset_locator.inset_axes(ax, width = "25%", height = "50%", loc = 10)

        ax_small.set_xlim([172.0, 174.5])
        ax_small.set_ylim([5.2, 5.7])

        axis_ticks_styling!(ax_small)

        ax_small.plot(data["x"], data["M_s"], color = "darkblue")
        ax_small.axhline(5.532955671044105, linestyle = "--", color = "grey")


        ax_small.set_xticklabels([])
        ax_small.set_yticklabels([])
        inset_locator.mark_inset(ax, ax_small, 2, 4, alpha = 0.8, zorder = 1000.0, linestyle = ":")

        axhline(M_A_ideal, linestyle = "--", color = "grey")

    subplots_adjust(hspace = 0, wspace = 0.3)

    @info "  saving..."
    savefig(plot_name,
        bbox_inches = "tight")
    close(fig)

    @info "done!"
end


fi = "data/tests/injection/NLDSA/nldsa_recalc/snap_000"
plot_name = "Plots/Fig07.png"
plot_cluster_cr_shocktube_2x2_plots(fi, plot_name)
