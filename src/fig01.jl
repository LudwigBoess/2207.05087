# activate the environment
using Pkg; Pkg.activate(".")

# load packages
using PyPlot, PyPlotUtility
using DSAModels

"""
    plot_dsa_models(plot_name)

Run plot of DSA models, Fig. 01
"""
function plot_dsa_models(plot_name)

    # range of Mach numbers
    M = LinRange(1.0, 100.0, 1000)

    # colors of models
    acc_colors = ["#000c8f", "#a700c9", "#c90000", "#f28500"]
    model_names = ["Kang+ 2007", 
                   "Kang&Ryu 2013", 
                   "Caprioli&Spitkovsky 2014", 
                   "Ryu+ 2019"]
    dsa_models = [0, 1, 3, 2]

    # plot styling
    lw = 3.0
    x_pixels = 700
    fig = get_figure(; x_pixels)
    plot_styling!(x_pixels)
    ax = gca()

    axis_ticks_styling!(ax)
    ax.set_xlim([1.0, 110.0])
    ax.set_ylim([0.001, 1.0])
    ax.set_xscale("log")
    ax.set_yscale("log")

    xlabel("Sonic Mach Number  " * L"\mathcal{M}_s")
    ylabel("Acceleration efficiency  " * L"η" * "  of model")

    # allocate storage arrays
    η = Vector{Float64}(undef, length(M))
    η_r = Vector{Float64}(undef, length(M))

    # loop over models
    for m = 1:length(dsa_models)

        # define efficiency structs
        if dsa_models[m] == 0
            η_model = Kang07()
        elseif dsa_models[m] == 1
            η_model = KR13()
        elseif dsa_models[m] == 3
            η_model = CS14()
        elseif dsa_models[m] == 2
            η_model = Ryu19()
        end

        # compute efficiency models
        @inbounds for i = 1:length(M)
            η[i]   = η_Ms(η_model, M[i], 0.0)
            η_r[i] = η_Ms(η_model, M[i], 1.0)
        end

        plot(M[η.>0.0], η[η.>0.0], lw = lw, color = acc_colors[m], label = model_names[m])
        plot(M[η_r.>0.0], η_r[η_r.>0.0], lw = lw, color = acc_colors[m], linestyle = "--")

    end

    plot([NaN], [NaN], color = "k", label = "Acc. Efficiency")
    plot([NaN], [NaN], color = "k", linestyle = "--", label = "Reacc. Efficiency")

    legend(frameon = false)

    savefig(plot_name,
        bbox_inches = "tight")

    close(fig)

end

plot_name = "Plots/Fig01.pdf"
plot_dsa_models(plot_name)