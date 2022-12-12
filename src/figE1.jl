using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Unitful, UnitfulAstro
using Cosmology

const kpc = 3.085678e21
const Myr = 1.e6u"yr" |> u"s" |> ustrip


const α_κ   = 0.5

# redshift of CIZA
z = 0.1921

# Planck 2018
c = cosmology(h=0.674,
              OmegaM=0.315,
              Neff=2.99)

θ = 5.0u"arcsecond" |> u"arcminute"

# diameter of the 
const L_beam = arcmin_to_kpc(c, θ, z) |> u"cm" |> ustrip


const me = 9.10953e-28
const σT = 6.65245e-25
const cL = 2.9979e10
const uIC = 4.00544e-13 * (1 + z)^4
const qe        = 1.602176487e-20 * cL
const C_crit_p   = 3qe / ( 4π * me * cL ) # Donnert+16, MNRAS 462, 2014–2032 (2016), Eq. 20 
                                                        #  -> converted to dimensionless momentum

# unit conversion struct                                                       
GU = GadgetPhysical()

"""
    smallest_synch_bright_p(ν::Real, B::Real)

Calculate the smallest synchrotron bright moment following Donnert+16, Eq. 22
"""
function smallest_synch_bright_p(ν::Real, B::Real)
    √( ν / (C_crit_p * B) )
end

κ(p, κ_10k) = κ_10k * ( p * 1.e-4 )^α_κ

τ_cool(p, uB) = (3.0 / 4.0 * me^2 * cL^2) / (σT * (uB + uIC) * p) * GU.t_s

τ_diff(p, L, κ_10k) = L^2 / κ(p, κ_10k)


function plot_timescales(plot_name)

    GU = GadgetPhysical()

    Nsteps = 1_000
        
    B = 5.e-6
    uB = B^2/8π

    p_range = 10.0.^LinRange(2, 7, Nsteps)

    t_cool = Vector{Float64}(undef, Nsteps)
    t_diff = Vector{Float64}(undef, Nsteps)


    fig = get_figure(1.0, x_pixels=700)
    plot_styling!(700)

    gs = plt.GridSpec(2, 1, figure = fig, height_ratios = [1, 0.5])

    subplot(get_gs(gs, 0, 0))
        
        ax = gca()
        axis_ticks_styling!(ax)

        ax.set_xscale("log")
        ax.set_yscale("log")

        ylabel("Time Scale  [ Myr ]")
        ax.set_xticklabels([])

        for i = 1:Nsteps
            t_cool[i] = τ_cool(p_range[i], uB) / Myr
            t_diff[i] = τ_diff(p_range[i], L_beam)  / Myr
        end

        plot(p_range, t_diff, color="darkblue", lw=3, label=L"\tau_\mathrm{diff}")
        plot(p_range, t_cool, color="teal",     lw=3, label=L"\tau_\mathrm{cool}")

        for i = 1:Nsteps
            t_cool[i] = τ_cool(p_range[i], 0.0) / Myr
            t_diff[i] = τ_diff(p_range[i], L_beam)   / Myr
        end

        plot(p_range, t_diff, color="darkblue", linestyle="--", lw=3)
        plot(p_range, t_cool, color="teal",     linestyle="--", lw=3)

        axvline(smallest_synch_bright_p(144.0e6, B), color="gray", linestyle="--")
        axvline(smallest_synch_bright_p(1.4e9, B),   color="gray", linestyle=":")

        # dummy plots for labels 
        plot([0.0], [0.0], color="k", lw=3, label=L"B = 5\mu\mathrm{G}")
        plot([0.0], [0.0], color="k", lw=3, label=L"B = 0", linestyle="--")

        legend(frameon=false, loc="lower left", fontsize=15)


    subplot(get_gs(gs, 1, 0))
        ax = gca()
        axis_ticks_styling!(ax)
        
        ax.set_xscale("log")
        ax.set_yscale("log")

        xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_e c)^{-1}" * "]")
        ylabel(L"\tau_\mathrm{diff}/\tau_\mathrm{cool}")

        for i = 1:Nsteps
            t_cool[i] = τ_cool(p_range[i], uB) / Myr
            t_diff[i] = τ_diff(p_range[i], L_beam)  / Myr
        end

        plot(p_range, t_diff ./ t_cool, color="k", lw=3)

        for i = 1:Nsteps
            t_cool[i] = τ_cool(p_range[i], 0.0) / Myr
            t_diff[i] = τ_diff(p_range[i], L_beam)   / Myr
        end

        plot(p_range, t_diff ./ t_cool, color="k", linestyle="--", lw=3)

        axvline(smallest_synch_bright_p(144.0e6, B), label=L"p_\mathrm{c}(144\mathrm{ MHz})", color="gray", linestyle="--")
        axvline(smallest_synch_bright_p(1.4e9, B),   label=L"p_\mathrm{c}(1.4\mathrm{ GHz})", color="gray", linestyle=":")
    
        legend(frameon=false, loc="lower right")
    subplots_adjust(hspace = 0)

    savefig(plot_name, bbox_inches="tight")
    close(fig)


end

plot_name = "/e/ocean2/users/lboess/Paper/CRESCENDO/Plots/FigE1.png"

plot_timescales(plot_name)


function plot_timescales_compare(plot_name)

    GU = GadgetPhysical()

    # as in Ogrodnik+20
    κ_10k_ogrodnik = 3.e26

    # as in Trotta+11
    κ_10k_trotta = 6.e28

    # as in Chen+19
    κ_10k_chan = 3.e29

    κ_10k    = [κ_10k_ogrodnik, κ_10k_trotta, κ_10k_chan]
    κ_labels = ["Ogrodnik+20", "Trotta+11", "Chan+19"]
    κ_labels = [L"\kappa = 3 \cdot 10^{26}", L"\kappa = 6 \cdot 10^{28}", L"\kappa = 3 \cdot 10^{29}"]

    Nsteps = 1_000
        
    B = 5.e-6
    uB = B^2/8π

    p_range = 10.0.^LinRange(2, 7, Nsteps)

    t_diff = Vector{Float64}(undef, Nsteps)

    t_cool_B  = [τ_cool(p, uB)  / Myr for p ∈ p_range]
    t_cool_IC = [τ_cool(p, 0.0) / Myr for p ∈ p_range]





    # define colors
    colors = [ "k", #"darkblue", 
            "dodgerblue", "teal",  "#82f70c",]


    fig = get_figure(1.0, x_pixels=700)
    plot_styling!(700, legend_font_size=12)

    gs = plt.GridSpec(2, 1, figure = fig, height_ratios = [1, 0.5])

    subplot(get_gs(gs, 0, 0))
        
        ax1 = gca()
        axis_ticks_styling!(ax1)

        ax1.set_xscale("log")
        ax1.set_yscale("log")

        ax1.set_ylabel("Time Scale  [ Myr ]")
        ax1.set_xticklabels([])

        # plot cooling times
        ax1.plot(p_range, t_cool_B,  color=colors[1], lw=3, label=L"\tau_\mathrm{cool}")
        ax1.plot(p_range, t_cool_IC, color=colors[1], lw=3, linestyle="--")

    subplot(get_gs(gs, 1, 0))
        ax2 = gca()
        axis_ticks_styling!(ax2)
        
        ax2.set_xscale("log")
        ax2.set_yscale("log")

        ax2.set_ylim([0.1, 1.e4])

        ax2.set_xlabel("Dimensionless Momentum  " * L" \hat{p}" * " [ " * L" (m_e c)^{-1}" * "]")
        ax2.set_ylabel(L"\tau_\mathrm{diff}/\tau_\mathrm{cool}")

        

        # tau_diff = tau_cool line
        ax2.axhline(1.0, color="lightgray")


    for i = 1:length(κ_10k)

        for j = 1:Nsteps
            t_diff[j] = τ_diff(p_range[j], L_beam, κ_10k[i])  / Myr
        end

        # upper plot
        ax1.plot(p_range, t_diff, color=colors[i+1], lw=3, label=κ_labels[i])
        
        # lower plot
        ax2.plot(p_range, t_diff ./ t_cool_B,  color=colors[i+1], lw=3)
        ax2.plot(p_range, t_diff ./ t_cool_IC, color=colors[i+1], lw=3, linestyle="--")
        
    end
        

        # dummy plots for labels 
        ax2.plot([0.0], [0.0], color="k", lw=3, label=L"B = 5\mu\mathrm{G}")
        ax2.plot([0.0], [0.0], color="k", lw=3, label=L"B = 0", linestyle="--")

        

        leg1 = ax1.legend(frameon=false, loc="lower left")

        ax2.legend(frameon=false, loc="lower right")

        # synchrotron momenta
        ax2.axvline(smallest_synch_bright_p(144.0e6, B), color="gray", linestyle="--")
        ax2.axvline(smallest_synch_bright_p(1.4e9, B),   color="gray", linestyle=":")

        # synchrotron momenta
        l1 = ax1.axvline(smallest_synch_bright_p(144.0e6, B), 
                    label=L"p_\mathrm{c}(144\mathrm{ MHz})", 
                    color="gray", linestyle="--")
        l2 = ax1.axvline(smallest_synch_bright_p(1.4e9, B),   
                    label=L"p_\mathrm{c}(1.4\mathrm{ GHz})", 
                    color="gray", linestyle=":")

        leg2 = ax1.legend([l1, l2], [L"p_\mathrm{c}(144\mathrm{ MHz})", L"p_\mathrm{c}(1.4\mathrm{ GHz})"], frameon=false, loc="upper right")
        ax1.add_artist(leg1)
        ax1.add_artist(leg2)

    subplots_adjust(hspace = 0)

    savefig(plot_name, bbox_inches="tight")
    close(fig)


end

plot_name = "/e/ocean2/users/lboess/Paper/CRESCENDO/Plots/FigE1.pdf"

plot_timescales_compare(plot_name)
