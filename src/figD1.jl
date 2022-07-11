# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using AnalyticMHDTestSolutions
using SpectralCRsUtility
using Printf
using PyCall
using ProgressMeter
using Base.Threads
patheffects = pyimport("matplotlib.patheffects")
inset_locator = pyimport("mpl_toolkits.axes_grid1.inset_locator")

"""
    get_synch_emission(fi)

Calculate the synchrotron emissivity of all particles.
"""
function get_synch_emission(fi)

    # cr setup 
    pmin = 1.e-1
    pmax = 1.e5
    Nbins = 48

    par = CRMomentumDistributionConfig(pmin, pmax, Nbins)

    GU = GadgetPhysical(hpar = 1.0)

    norm  = GU.CR_norm .* 10.0.^read_block(fi, "CReN", parttype=0)
    slope = read_block(fi, "CReS", parttype=0) .|> Float64
    cut   = read_block(fi, "CReC", parttype=0) .|> Float64

    Npart = length(cut)

    j_ν = Vector{Float64}(undef, Npart)

    @showprogress for i = 1:Npart 
        j_ν[i] = synchrotron_emission(norm[:,i], slope[:,i], cut[i], 5.e-6, par, ν0 = 1.44e9, 
                                      reduce_spectrum=true,
                                      integrate_pitch_angle = true)
    end

    return j_ν
end

"""
    get_ciza_data(fi)

Load the simulation data
"""
function get_ciza_data(fi)

    h = head_to_obj(fi)

    # read positions
    x = read_block(fi, "POS", parttype=0)

    # sort along x direction
    j = sortperm(x[1,:])
    x = x[1,j]

    γ_th = 5.0/3.0

    hsml = read_block(fi, "HSML", parttype=0)[j]

    ρ = read_block(fi, "RHO", parttype=0)[j]

    u = read_block(fi, "U", parttype=0)[j]

    P_th = ( γ_th - 1.0 ) .* ρ .* u

    MACH = read_block(fi, "MACH", parttype=0)[j]

    CRpP = read_block(fi, "CRpP", parttype=0)[j]

    CReP = read_block(fi, "CReP", parttype=0)[j]

    @info "calculating synch. emission"
    j_ν = get_synch_emission(fi)[j]
    @info "  done"

    mach_selected = findall(50.0 .<= x .< 100.0)
    mach_max_pos = findmax(MACH[mach_selected])[2]
    mach_max = MACH[mach_selected[mach_max_pos]]
    mach_max_x = x[mach_selected[mach_max_pos], 1]
    mean_hsml = hsml[mach_selected[mach_max_pos]]

    return x, ρ, u, MACH, P_th, CRpP, CReP, j_ν, mach_max, mach_max_x, mean_hsml
end

"""
    plot_ciza_shock_4x1(fi, plot_name, mach_ideal)

Main plot routine.
"""
function plot_ciza_shock_4x1(fi, plot_name, mach_ideal)

    # ciza shock properties
    GU = GadgetPhysical()
    M  = 4.6
    T2 = 2.588e8 # K 
    ρ2 = 9.4e-27 # g/cm^3
    γ  = 5/3

    xmin = 50.0
    xmax = 100.0

    ρr = ρ2 / ( (γ + 1 ) * M^2 / ( (γ - 1 ) * M^2 + 2) ) / GU.rho_cgs
    Ur = T2 / ( (2γ*M^2 - (γ - 1))*( (γ-1)*M^2 + 2) / ( (γ+1)^2 * M^2) ) / GU.T_K
    Pr = (γ-1) * ρr * Ur

    ρl = 8ρr
    Pl = AnalyticMHDTestSolutions.solvePlfromMach(ρl, ρr, Pr, M, γ)

    h = head_to_obj(fi)
    time_end = h.time

    par = RiemannParameters(rhol=ρl, Pl=Pl, rhor=ρr, Pr=Pr, t=time_end)

    # xrange for ideal solution
    x_ideal = collect(LinRange(xmin, xmax, 1000))

    x_pixels = 600
    fig = get_figure(4.5; x_pixels)
    plot_styling!(x_pixels)

    x, ρ, u, MACH, P_th, CRpP, CReP, j_ν, mach_max, mach_max_x, mean_hsml = get_ciza_data(fi)

    P_tot = P_th .+ CRpP .+ CReP

    sol = AnalyticMHDTestSolutions.solve(x_ideal, par)

    println(sol.Mach)

    # config
    mach_lim = 5.0
    text_x = 55.0
    text_ideal_y = mach_ideal - 0.2 * mach_ideal
    text_max_y = mach_ideal - 0.1 * mach_ideal


    subplot(1, 4, 1)
        ax12 = gca()
        ax12.set_xlim([xmin, xmax])
        ax12.set_ylim([0.0, 2.5])
        axis_ticks_styling!(ax12)

        xlabel("Position  " * L"x" * " [kpc]")
        ylabel("Density  " * L"ρ" * " [" * L"10^{-26}" * "g cm" * L"^{-3}" * "]")

        plot(x, 1.e26 .* ρ .* GU.rho_cgs)
        plot(sol.x, 1.e26 .* sol.rho .* GU.rho_cgs, linestyle = "--", color = "red", alpha = 0.8)

        plot([0.0], [0.0], color = "k", label = "Simulation")
        plot([0.0], [0.0], color = "k", linestyle = "--", label = "Analytic")

        legend(frameon = false, loc = "upper right")

    subplot(1, 4, 2)
        ax = gca()
        ax.set_xlim([xmin, xmax])

        axis_ticks_styling!(ax)

        xlabel("Position  " * L"x" * " [kpc]")
        ylabel("Pressure  " * L"P" * " [" * L"10^{-9}" * "erg cm" * L"^{-3}" * "]")

        fac = 1.e9

        plot(x, fac .* CReP .*  GU.P_CR_cgs, color = "orange", label = L"P_{CR,e}") #  color="#FEDB4A"
        plot(x, fac .* P_th .* GU.P_th_cgs, color = "red", label = L"P_{\mathrm{th}}")
        plot(sol.x, fac .* sol.P .* GU.P_th_cgs, "--", color = "red")
        plot(x, fac .* P_tot .* GU.P_th_cgs, color = "black", label = L"P_{\mathrm{tot}}")
        plot(sol.x, fac .* sol.P .* GU.P_th_cgs, "--", color = "black")

        legend(frameon = false, loc = "lower left")

        # inset plot
        ax_small = inset_locator.inset_axes(ax, width = "50%", height = "50%", loc = "upper right")
        ax_small.set_xlim([82.0, 92.0])
        ax_small.set_ylim(fac .* [1.e-14, 0.5e-9])
        ax_small.set_yscale("log")


        axis_ticks_styling!(ax_small)
        ax_small.set_xticklabels([])
        ax_small.plot(x, fac .* CReP .* GU.P_CR_cgs, color = "orange", label = L"P_{CR,e}") #  color="#FEDB4A"
        ax_small.plot(x, fac .* P_th .* GU.P_th_cgs, color = "red", label = L"P_{\mathrm{th}}")
        ax_small.plot(sol.x, fac .* sol.P .* GU.P_th_cgs, "--", color = "red")
        ax_small.plot(x, fac .* P_tot .* GU.P_th_cgs, color = "black", label = L"P_{\mathrm{tot}}")
        ax_small.plot(sol.x, fac .* sol.P .* GU.P_th_cgs, "--", color = "black")

        inset_locator.mark_inset(ax, ax_small, 3, 4, alpha = 0.8, zorder = 1000.0, linestyle = ":", lw = 1.5)

    subplot(1, 4, 3)

        ax4 = gca()
        ax4.set_xlim([xmin, xmax])
        ax4.set_ylim([1.e-41, 1.e-39])
        ax4.set_yscale("log")
        axis_ticks_styling!(ax4)

        xlabel("Position  " * L"x" * " [kpc]")
        ylabel("Synch. Emissivity  " * L"j_{ν,1.4\mathrm{ Ghz}}" * " [ erg cm" * L"^{-3}" * " s" * L"^{-1}" * " Hz" * L"^{-1}" * " ]")

        plot(x, j_ν, label="Simulation")
        #plot(x, j_nu_output, label="Gadget output")
        axhline(1.6301033678680908e-40, linestyle = "--", color = "grey", label="Donnert+17")
        legend(frameon = false, loc = "lower left")

    subplot(1, 4, 4)

        ax4 = gca()
        ax4.set_xlim([xmin, xmax])
        ax4.set_ylim([0.0, mach_lim])
        axis_ticks_styling!(ax4)

        xlabel("Position  " * L"x" * " [kpc]")
        ylabel("Sonic Mach Number  " * L"\mathcal{M}_s")

        plot(x, MACH)

        axhline(sol.Mach, linestyle = "--", color = "grey")

        errorbar(mach_max_x, mach_max, xerr = mean_hsml,
            ecolor = "r", capsize = 10,
            label = "Kernel Size")
        legend(frameon = false, loc = "lower left")

        text(text_x, text_ideal_y, L"\mathcal{M}_{\mathrm{ideal}} = " * "$(@sprintf("%0.3f", sol.Mach))")

        text(text_x, text_max_y, L"\mathcal{M}_{\mathrm{max}} = " * "$(@sprintf("%0.3f", mach_max))")

    subplots_adjust(wspace=0.25)

    savefig(plot_name,
        bbox_inches = "tight")
    close(fig)
end

path_in = "data/tests/injection/CIZA_shock/"
fi = path_in * "snap_001"

plot_name = "Plots/FigD1.pdf"
mach_ideal = 4.6

plot_single_efficiency_4x1_plots(fi, plot_name, mach_ideal)
