# activate the environment
using Pkg; Pkg.activate(".")

using PyPlot, PyPlotUtility
using Statistics
using Printf
using DelimitedFiles
using LsqFit

function get_performance(file)

    f = open(file)
    lines = readlines(f)
    close(f)

    stepcounter = Array{Int64,1}(undef,length(lines[32:end]))
    timing = Array{Float64,1}(undef,length(lines[32:end]))
    active = Array{Int64,1}(undef,length(lines[32:end]))

    for i = 1:(length(lines)-31)
        stepcounter[i] = parse(Int64,strip(lines[i+31][6:12]))
        timing[i] = parse(Float64,strip(lines[i+31][19:28]))
        active[i] = parse(Int64,strip(lines[i+31][34:44]))
    end

    return stepcounter, timing, active
end

function get_total_time(file)

    stepcounter, timing, active = get_performance(file)

    return sum(timing)#/60.0
end

function get_timing(file, N_cpus)

    stepcounter, timing, active = get_performance(file)
    t_end = sum(timing)

    t = zeros(length(timing))

    for i = 1:length(timing)
        t[i] = sum(timing[1:i])/t_end
    end

    t_per_particle = active ./ timing
    t_per_particle ./= N_cpus

    return t, t_per_particle
end

function fit_data(t, tp)

    N_bins = 50
    N_datapoints = length(t)
    i_size = floor(N_datapoints/N_bins)

    t_fit = zeros(N_bins)
    tp_fit = zeros(N_bins)
    tp_min = zeros(N_bins)

    for i = 0:N_bins-2
        t_fit[i+1] = (t[Int(1+i*i_size)] + t[Int((i+1)*i_size)])/2.0
        tp_fit[i+1] = mean(tp[Int(1+i*i_size):Int((i+1)*i_size)])
        tp_min[i+1] = minimum(tp[Int(1+i*i_size):Int((i+1)*i_size)])
    end

    i = N_bins-1
    t_fit[i+1] = (t[Int(1+i*i_size)] + t[end])/2.0
    tp_fit[i+1] = mean(tp[Int(1+i*i_size):end])
    tp_min[i+1] = minimum(tp[Int(1+i*i_size):end])

    return t_fit, tp_fit, tp_min
end

function get_curves(file, N_cpus)

    t, tp = get_timing(file, N_cpus)

    t_fit, tp_fit, tp_min = fit_data(t, tp)

    return t, tp, t_fit, tp_fit, tp_min
end



function plot_performance(files, colors_data, colors_fit, labels, title_text, plot_file, N_cpus)

    N_files = length(files)

    t = Dict()
    tp = Dict()
    t_fit = Dict()
    tp_fit = Dict()
    tp_min = Dict()
    t_tot = Dict()

    for i = 1:N_files
        t[i], tp[i], t_fit[i], tp_fit[i], tp_min[i] = get_curves(files[i], N_cpus)
        t_tot[i] = get_total_time(files[i])
    end


    axis_label_font_size = 20
    legend_font_size = 15
    tick_label_size = 15
    rc("font", family="serif")
    rc("mathtext", fontset="stix")
    #rc("lines", color="#003366")

    upscale = 0.7
    fig = figure(figsize=(10*upscale, 10*upscale))

    ax = gca()
    ax.set_yscale("log")
    ax.set_xlim([0.0,1.0])
    ax.set_ylim([1.0e3,1.0e5])
    ax.tick_params(reset=true, direction="in", axis="both",labelsize=tick_label_size,
                which="major", size=12, width=1)
    ax.tick_params(reset=true, direction="in", axis="both",labelsize=tick_label_size,
                which="minor", size=6, width=1)

    ax.minorticks_on()
    xlabel(L"t/t_{max}", fontsize=axis_label_font_size)
    ylabel(L"N_{Part}/s/CPU", fontsize=axis_label_font_size)


    # plot background
    for i = 1:N_files
        plot(t[i], tp[i], color=colors_data[i], alpha=0.4)
    end

    # plot fits
    for i = 1:N_files
        plot(t_fit[i], tp_fit[i], color=colors_fit[i], label= labels[i] * L"t_{tot}" * " = $(@sprintf("%0.2f", t_tot[i] )) min -> $(@sprintf("%0.2f", t_tot[i]/t_tot[1] ))x")
        #plot(t_fit[i], tp_min[i], color=colors_fit[i], linestyle="--")
    end

    legend(frameon=false, fontsize=12, loc="upper right")# bbox_to_anchor=(1, 0.75))
    title(title_text, fontsize=legend_font_size)
    savefig(plot_file, bbox_inches="tight")
    close(fig)

end


# Adiabatic changes
file0 = "brent_24.txt"
file1 = "newton_24.txt"
file2 = "brent_192.txt"
file3 = "newton_192.txt"
files = "data/tests/performance/" .* [file0, file1, file2, file3]

colors_data = ["#003366", "#0d9828",
            "#a900cc", "#0083cc"
            ]

colors_fit = ["#003366", "#096c1c",
            "#a900cc", "#0083cc"
            ]

labels = ["24 Bins - Brent\t\t\t    - ",
          "24 Bins - Bisect.+Newton    - ",
          "192 Bins - Brent\t\t   - ",
          "192 Bins - Bisect.+Newton  - "]

title_text =  "4 MPI ranks - 4 OpenMP threads"
plot_file = "Plots/FigB1.pdf"
N_cpus = 4
plot_performance(files, colors_data, colors_fit, labels, title_text, plot_file, N_cpus)


