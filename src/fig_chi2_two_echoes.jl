using Random, CairoMakie
import fcVARPRO as FCV

# =============== χ² for two echoes =======================

B0 = 3
ppm_fat = [-3.80, -3.40, -2.60, -1.95, -0.5, 0.60]
ampl_fat = [0.0875, 0.6998, 0.1206, 0.0062, 0.0389, 0.0471]

ϕs = [range(-π, π, 3601);]

rng = MersenneTwister()
nTE = 2
t0 = 1.02
ΔTE = 1.07
TEs = t0 .+ [ΔTE .* (0:nTE-1);]
c_t = 1.0
σ = 0.0
n_σ = 1
ϕ_t = 0
R2ss_t = [1/20]
fs_t = [0.02]

res_MC = FCV.MC_sim(;
    rng=rng,
    TEs=TEs,
    B0=B0,
    ppm_fat=ppm_fat,
    ampl_fat=ampl_fat,
    ϕ_t=ϕ_t,
    R2ss_t=R2ss_t,
    fs_t=fs_t,
    c_t=c_t,
    σ=σ,
    n_σ=n_σ,
    ϕs=ϕs,
    R2ss=[0.0],
)

inch = 96
pt = 4 / 3

width, height = 3.42inch, 3.42inch

fig = Figure(size=(width, height), fontsize=12pt)

ax = Axis(fig[1, 1],
    title=L"$\log_{10}\,\left[\,\chi^2(\varphi,\, R_2^\ast = 0)\,\right]$",
    titlesize=12pt,
    xlabel=L"$\varphi$ [rad]",
    xlabelsize=12pt,
    xticks=([-π, 0.0, π], ["-π", "0", "π"]),
    xticklabelsize=10pt,
    yticklabelsize=10pt,
)

rw = lines!(ax, ϕs, log10.(res_MC.χ2[:RW][:]), color=:red)
fc = lines!(ax, ϕs, log10.(res_MC.χ2[:FC][:]), color=:blue)
xlims!(ax, -π, π)
ylims!(ax, -6, 0)
axislegend(ax, [fc, rw], ["FC", "RW"], labelsize=10pt, position=:lb)

display(fig)
##

fig_name = "fig_4"
save(fig_name * ".svg", fig)
save(fig_name * ".eps", fig)
run(`epspdf $fig_name".eps"`)