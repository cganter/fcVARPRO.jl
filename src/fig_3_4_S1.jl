using Revise, HDF5, CairoMakie, Random, LinearAlgebra
import VP4Optim as VP
import B0Map as BM
import fcVARPRO as FCV

BLAS.set_num_threads(1)

# read the HDF5 file
fid = h5open("data/two_echoes/20241024_171954_702_ImDataParamsBMRR_subspace2comp_wfi.h5", "r");
obj_data = read(fid["ImDataParams"]);

signal = obj_data["signal"];
ss = size(signal);
data = zeros(ComplexF64, ss[3:5]..., ss[1])
for i in 1:2
    data[:, :, :, i] .= signal[i, 1, :, :, :]
end

# the supplied mask is too inclusive
# the following choice is better but far from perfect..
# the choice is insofar important as 
S = abs.(data[:, :, :, 2]) .> 0.5;

# echo times
TEs = 1000obj_data["TE_s"];  # the expected unit is [ms]

# field strength
B0 = Float64(attrs(fid["ImDataParams"])["fieldStrength_T"])

# fat model
ppm_fat = read(fid["AlgoParams"]["FatModel"]["freqs_ppm"])
ampl_fat = read(fid["AlgoParams"]["FatModel"]["relAmps"])

# close the HDF5 file
close(fid)

# scanner-dependent convention for the orientation of precession
precession = :counterclockwise

# set up GRE parameters
grePar_FC = VP.modpar(BM.GREMultiEchoWF;
    ts=TEs,
    B0=B0,
    ppm_fat=ppm_fat,
    ampl_fat=ampl_fat,
    precession=precession)

grePar_RW = VP.modpar(BM.GREMultiEchoWFRW;
    ts=TEs,
    B0=B0,
    ppm_fat=ppm_fat,
    ampl_fat=ampl_fat,
    precession=precession)

# we select a single slice ...
sl = 64
S_ = S[:, :, sl]
not_S = (!).(S_)
data_ = data[:, :, sl, :]

# generate fit parameters
fitpar_FC = BM.fitPar(grePar_FC, deepcopy(data_), deepcopy(S_))
fitpar_RW = BM.fitPar(grePar_RW, deepcopy(data_), deepcopy(S_))

# generate fit options
fitopt = BM.fitOpt()
fitopt.R2s_rng = [0.0, 0.0]
fitopt.optim = false # not needed, since R2s = 0 is fixed

fitopt_FC = deepcopy(fitopt)
fitopt_RW = deepcopy(fitopt)

# local fit
res_FC = BM.local_fit!(fitpar_FC, fitopt_FC)
res_RW = BM.local_fit!(fitpar_RW, fitopt_RW)

# extract maps
ϕ_RW = fitpar_RW.ϕ
ϕ_RW[not_S] .= NaN
f_RW = BM.fat_fraction_map(fitpar_RW, fitopt_RW)
f_RW[not_S] .= NaN
χ2_fit_RW = fitpar_RW.χ2
χ2_fit_RW[not_S] .= NaN
freq_RW = BM.freq_map(fitpar_RW);

ϕ_FC = fitpar_FC.ϕ
ϕ_FC[not_S] .= NaN
f_FC = BM.fat_fraction_map(fitpar_FC, fitopt_FC)
f_FC[not_S] .= NaN
χ2_fit_FC = fitpar_FC.χ2
χ2_fit_FC[not_S] .= NaN
freq_FC = BM.freq_map(fitpar_FC);

coil_FC = map(x -> x[1], fitpar_FC.c)
in_phase_FC = abs.(coil_FC)
coil_phase_FC = angle.(coil_FC);

## =============== generate phase maps =======================

r_op, c_op = 80, 45
r_in, c_in = 85, 33
kwargs_op = (; markersize=10, marker=:circle, strokewidth=2, color=(:red, 0.0), strokecolor=:red)
kwargs_in = (; markersize=10, marker=:diamond, strokewidth=2, color=(:red, 0.0), strokecolor=:red)

# size definitions
# these are relative to 1 CSS px
inch = 96
pt = 4 / 3

width, height = 6.9inch, 4.2inch
colmapO = :romaO

fig = Figure(size=(width, height), fontsize=12pt)

ax = Axis(fig[1, 1],
    title=L"$\varphi_{RW}$",
)

heatmap!(ax,
    ϕ_RW',
    colormap=colmapO,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r_op, c_op, ϕ_RW, transpose)...; kwargs_op...)
scatter!(ax, FCV.rc2xy(r_in, c_in, ϕ_RW, transpose)...; kwargs_in...)
hidedecorations!(ax)
Label(fig[1, 1, TopLeft()], "A",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[2, 1],
    title=L"$\varphi_{FC}$",
)

heatmap!(ax,
    ϕ_FC',
    colormap=colmapO,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r_op, c_op, ϕ_FC, transpose)...; kwargs_op...)
scatter!(ax, FCV.rc2xy(r_in, c_in, ϕ_FC, transpose)...; kwargs_in...)
hidedecorations!(ax)
Label(fig[2, 1, TopLeft()], "D",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig[1, 2],
    colorrange=(-π, π),
    colormap=colmapO,
    ticks=([-π, 0.0, π], ["-π", "0", "π"]),
    ticklabelsize=8pt,
)

Colorbar(fig[2, 2],
    colorrange=(-π, π),
    colormap=colmapO,
    ticks=([-π, 0.0, π], ["-π", "0", "π"]),
    ticklabelsize=8pt,
)

# =============== generate fat fraction maps =======================

colmap = :imola

ax = Axis(fig[1, 3],
    title=L"$f_{RW}$",
)

heatmap!(ax,
    f_RW',
    colormap=colmap,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r_op, c_op, f_RW, transpose)...; kwargs_op...)
scatter!(ax, FCV.rc2xy(r_in, c_in, f_RW, transpose)...; kwargs_in...)

hidedecorations!(ax)

Label(fig[1, 3, TopLeft()], "B",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[2, 3],
    title=L"$f_{FC}$",
)

heatmap!(ax,
    f_FC',
    colormap=colmap,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r_op, c_op, f_FC, transpose)...; kwargs_op...)
scatter!(ax, FCV.rc2xy(r_in, c_in, f_FC, transpose)...; kwargs_in...)

hidedecorations!(ax)

Label(fig[2, 3, TopLeft()], "E",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig[1, 4],
    colorrange=(0, 1),
    colormap=colmap,
    ticklabelsize=8pt,
)

Colorbar(fig[2, 4],
    colorrange=(0, 1),
    colormap=colmap,
    ticklabelsize=8pt,
)

# =============== generate χ² plots =======================

ϕs_HR = [range(-π, π, 3601);]
ϕs = ϕs_HR[1:4:end]
R2ss = [range(0, 0.4, 101);]

data_op = vec(data[r_op, c_op, sl, :])
data_in = vec(data[r_in, c_in, sl, :])

res_op_HR = FCV.fit_data_loc(;
    TEs=fitpar_FC.grePar.ts,
    B0=fitpar_FC.grePar.B0,
    data=data_op,
    ϕs=ϕs_HR,
    R2ss=R2ss,
)
res_in_HR = FCV.fit_data_loc(;
    TEs=fitpar_FC.grePar.ts,
    B0=fitpar_FC.grePar.B0,
    data=data_in,
    ϕs=ϕs_HR,
    R2ss=R2ss,
)
res_op = FCV.fit_data_loc(;
    TEs=fitpar_FC.grePar.ts,
    B0=fitpar_FC.grePar.B0,
    data=data_op,
    ϕs=ϕs,
    R2ss=R2ss,
)
res_in = FCV.fit_data_loc(;
    TEs=fitpar_FC.grePar.ts,
    B0=fitpar_FC.grePar.B0,
    data=data_in,
    ϕs=ϕs,
    R2ss=R2ss,
)

lχ_lim = (-5, 0.2)

ax = Axis(fig[1, 5],
    title=L"$\log_{10}\,\left[\,\chi^2_{RW}\,\right]$",
)

heatmap!(ax,
    log10.(abs.(χ2_fit_RW))',
    colormap=:batlow,
    colorrange=lχ_lim,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r_op, c_op, χ2_fit_RW, transpose)...; kwargs_op...)
scatter!(ax, FCV.rc2xy(r_in, c_in, χ2_fit_RW, transpose)...; kwargs_in...)

hidedecorations!(ax)

Label(fig[1, 5, TopLeft()], "C",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[2, 5],
    title=L"$\log_{10}\,\left[\,\chi^2_{FC}\,\right]$",
)

heatmap!(ax,
    log10.(abs.(χ2_fit_FC))',
    colormap=:batlow,
    colorrange=lχ_lim,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r_op, c_op, χ2_fit_FC, transpose)...; kwargs_op...)
scatter!(ax, FCV.rc2xy(r_in, c_in, χ2_fit_FC, transpose)...; kwargs_in...)

hidedecorations!(ax)

Label(fig[2, 5, TopLeft()], "F",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig[1, 6],
    colorrange=lχ_lim,
    colormap=:batlow,
    ticklabelsize=8pt,
)

Colorbar(fig[2, 6],
    colorrange=lχ_lim,
    colormap=:batlow,
    ticklabelsize=8pt,
)

display(fig)

## Supporting information ==============================================================

width, height = 6.9inch, 6.9inch
colmapO = :romaO

fig_4 = Figure(size=(width, height), fontsize=12pt)
fig_S1 = Figure(size=(width, height), fontsize=12pt)

colmap = :turbid
colmap_f = :berlin
lim = (-5, -2)

f_FC = deepcopy(res_in.f[:FC])
f_RW = deepcopy(res_in.f[:RW])
la_Δf = log.(abs.(f_RW .- f_FC) .+ eps())
S_no_f = f_FC .== 0 .|| f_FC .== 1
S_1 = deepcopy(S_no_f)
S_1[1:end-1, :] .= S_1[1:end-1, :] .&& (!).(S_1)[2:end, :]
S_2 = deepcopy(S_no_f)
S_2[2:end, :] .= S_2[2:end, :] .&& (!).(S_2)[1:end-1, :]
S_no_f = S_1 .|| S_2
ciS = CartesianIndices(S_no_f)
ciS = filter(ci -> S_no_f[ci], ciS)
xs = map(ci -> ϕs[ci[1]], ciS)
ys = map(ci -> R2ss[ci[2]], ciS)
χ2_FC = res_in.χ2[:FC]
χ2_RW = res_in.χ2[:RW]
laχ2_FC = log10.(abs.(χ2_FC) .+ eps())
laχ2_RW = log10.(abs.(χ2_RW) .+ eps())
f_RW[res_in.wf_opp] .= -f_RW[res_in.wf_opp]

# --------------------------------------------------------------

ax = Axis(fig_4[1, 1],
    title=L"$\log_{10}\,\left[\,\chi^2_{RW}\,\right]$",
    titlesize=12pt,
    ylabel=L"$R_2^\ast$ \;[1/ms$\,$]",
    yticklabelsize=8pt,
    xticks=([fitpar_RW.ϕ[r_in, c_in]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:red,
)

heatmap!(ax, ϕs, R2ss, laχ2_RW, colormap=colmap, colorrange=lim)
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, -0.55π, 0.1; kwargs_in...)

Label(fig_4[1, 1, TopLeft()], "A",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ax = Axis(fig_4[2, 1],
    title=L"$\log_{10}\,\left[\,\chi^2_{FC}\,\right]$",
    titlesize=12pt,
    ylabel=L"$R_2^\ast$ \;[1/ms$\,$]",
    yticklabelsize=8pt,
    xticks=([fitpar_FC.ϕ[r_in, c_in]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:blue,
)

heatmap!(ax, ϕs, R2ss, laχ2_FC, colormap=colmap, colorrange=lim)
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, -0.55π, 0.1; kwargs_in...)

Label(fig_4[2, 1, TopLeft()], "D",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ax = Axis(fig_S1[1, 1],
    title=L"$f_{RW}$",
    titlesize=12pt,
    ylabel=L"$R_2^\ast$ \;[1/ms$\,$]",
    xticks=([fitpar_RW.ϕ[r_in, c_in]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:red,
    yticklabelsize=8pt,
)

heatmap!(ax, ϕs, R2ss, f_RW, colormap=colmap_f, colorrange=(-1, 1))
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, 0.0, 0.2; kwargs_in...)

Label(fig_S1[1, 1, TopLeft()], "A",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ax = Axis(fig_S1[2, 1],
    title=L"$f_{FC}$",
    titlesize=12pt,
    ylabel=L"$R_2^\ast$ \;[1/ms$\,$]",
    xticks=([fitpar_FC.ϕ[r_in, c_in]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:blue,
    yticklabelsize=8pt,
)

heatmap!(ax, ϕs, R2ss, f_FC, colormap=colmap_f, colorrange=(-1, 1))
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, 0.0, 0.2; kwargs_in...)

Label(fig_S1[2, 1, TopLeft()], "D",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

f_FC = deepcopy(res_op.f[:FC])
f_RW = deepcopy(res_op.f[:RW])
la_Δf = log.(abs.(f_RW .- f_FC) .+ eps())
S_no_f = f_FC .== 0 .|| f_FC .== 1
S_1 = deepcopy(S_no_f)
S_1[1:end-1, :] .= S_1[1:end-1, :] .&& (!).(S_1)[2:end, :]
S_2 = deepcopy(S_no_f)
S_2[2:end, :] .= S_2[2:end, :] .&& (!).(S_2)[1:end-1, :]
S_no_f = S_1 .|| S_2
ciS = CartesianIndices(S_no_f)
ciS = filter(ci -> S_no_f[ci], ciS)
xs = map(ci -> ϕs[ci[1]], ciS)
ys = map(ci -> R2ss[ci[2]], ciS)
χ2_FC = res_op.χ2[:FC]
χ2_RW = res_op.χ2[:RW]
laχ2_FC = log10.(abs.(χ2_FC) .+ eps())
laχ2_RW = log10.(abs.(χ2_RW) .+ eps())
f_RW[res_op.wf_opp] .= -f_RW[res_op.wf_opp]

# --------------------------------------------------------------

ax = Axis(fig_4[1, 2],
    title=L"$\log_{10}\,\left[\,\chi^2_{RW}\,\right]$",
    titlesize=12pt,
    xticks=([fitpar_RW.ϕ[r_op, c_op]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:red,
)

heatmap!(ax, ϕs, R2ss, laχ2_RW, colormap=colmap, colorrange=lim)
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, -0.55π, 0.1; kwargs_op...)
arrows2d!(ax, [0.3π],[0.015],[-0.15π],[0],color=:white)
arrows2d!(ax, [0.75π],[0.015],[0.15π],[0],color=:white)
hideydecorations!(ax)

Label(fig_4[1, 2, TopLeft()], "B",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ax = Axis(fig_4[2, 2],
    title=L"$\log_{10}\,\left[\,\chi^2_{FC}\,\right]$",
    titlesize=12pt,
    xticks=([fitpar_FC.ϕ[r_op, c_op]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:blue,
)

heatmap!(ax, ϕs, R2ss, laχ2_FC, colormap=colmap, colorrange=lim)
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, -0.55π, 0.1; kwargs_op...)
arrows2d!(ax, [0.3π],[0.015],[-0.15π],[0],color=:white)
arrows2d!(ax, [0.75π],[0.015],[0.15π],[0],color=:white)
hideydecorations!(ax)

Label(fig_4[2, 2, TopLeft()], "E",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ax = Axis(fig_S1[1, 2],
    title=L"$f_{RW}$",
    titlesize=12pt,
    xticks=([fitpar_RW.ϕ[r_op, c_op]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:red,
)

heatmap!(ax, ϕs, R2ss, f_RW, colormap=colmap_f, colorrange=(-1, 1))
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, 0.0, 0.2; kwargs_op...)
hideydecorations!(ax)

Label(fig_S1[1, 2, TopLeft()], "B",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ax = Axis(fig_S1[2, 2],
    title=L"$f_{FC}$",
    titlesize=12pt,
    xticks=([fitpar_FC.ϕ[r_op, c_op]], [""]),
    xticksize=10,
    xtickwidth=2,
    xtickcolor=:blue,
)

heatmap!(ax, ϕs, R2ss, f_FC, colormap=colmap_f, colorrange=(-1, 1))
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, 0.0, 0.2; kwargs_op...)
hideydecorations!(ax)

Label(fig_S1[2, 2, TopLeft()], "E",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

# --------------------------------------------------------------

ϕ_t, R2s_t, f_t = fitpar_FC.ϕ[r_op, c_op], 0.05, 0.02

res_MC = FCV.MC_sim(;
    TEs=fitpar_FC.grePar.ts,
    B0=fitpar_FC.grePar.B0,
    ppm_fat=ppm_fat,
    ampl_fat=ampl_fat,
    ϕ_t=ϕ_t,
    R2ss_t=[R2s_t],
    fs_t=[f_t],
    c_t=abs.(fitpar_FC.c[r_op, c_op])[1],
    σ=0.0,
    n_σ=1,
    ϕs=ϕs,
    R2ss=R2ss,
)

f_FC = deepcopy(res_MC.f[:FC])
f_RW = deepcopy(res_MC.f[:RW])
la_Δf = log.(abs.(f_RW .- f_FC) .+ eps())
S_no_f = f_FC .== 0 .|| f_FC .== 1
S_1 = deepcopy(S_no_f)
S_1[1:end-1, :] .= S_1[1:end-1, :] .&& (!).(S_1)[2:end, :]
S_2 = deepcopy(S_no_f)
S_2[2:end, :] .= S_2[2:end, :] .&& (!).(S_2)[1:end-1, :]
S_no_f = S_1 .|| S_2
ciS = CartesianIndices(S_no_f)
ciS = filter(ci -> S_no_f[ci], ciS)
xs = map(ci -> ϕs[ci[1]], ciS)
ys = map(ci -> R2ss[ci[2]], ciS)
χ2_FC = res_MC.χ2[:FC]
χ2_RW = res_MC.χ2[:RW]
laχ2_FC = log10.(abs.(χ2_FC) .+ eps())
laχ2_RW = log10.(abs.(χ2_RW) .+ eps())
f_RW[res_MC.wf_opp] .= -f_RW[res_MC.wf_opp]

ax = Axis(fig_4[1, 3],
    title=L"$\log_{10}\,\left[\,\chi^2_{RW}\,\right]$",
    titlesize=12pt,
)

heatmap!(ax, ϕs, R2ss, laχ2_RW, colormap=colmap, colorrange=lim)
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, ϕ_t, R2s_t, markersize=8, color=:lime)
arrows2d!(ax, [0.3π],[0.015],[-0.15π],[0],color=:white)
arrows2d!(ax, [0.75π],[0.015],[0.15π],[0],color=:white)
hidedecorations!(ax)

Label(fig_4[1, 3, TopLeft()], "C",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig_4[1, 4],
    colorrange=lim,
    colormap=colmap,
    ticklabelsize=8pt,
)

# --------------------------------------------------------------

ax = Axis(fig_4[2, 3],
    title=L"$\log_{10}\,\left[\,\chi^2_{FC}\,\right]$",
    titlesize=12pt,
)

heatmap!(ax, ϕs, R2ss, laχ2_FC, colormap=colmap, colorrange=lim)
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, ϕ_t, R2s_t, markersize=8, color=:lime)
arrows2d!(ax, [0.3π],[0.015],[-0.15π],[0],color=:white)
arrows2d!(ax, [0.75π],[0.015],[0.15π],[0],color=:white)
hidedecorations!(ax)

Label(fig_4[2, 3, TopLeft()], "F",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig_4[2, 4],
    colorrange=lim,
    colormap=colmap,
    ticklabelsize=8pt,
)

# --------------------------------------------------------------

ax = Axis(fig_S1[1, 3],
    title=L"$f_{RW}$",
    titlesize=12pt,
)

heatmap!(ax, ϕs, R2ss, f_RW, colormap=colmap_f, colorrange=(-1, 1))
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, ϕ_t, R2s_t, markersize=8, color=:lime)
hidedecorations!(ax)

Label(fig_S1[1, 3, TopLeft()], "C",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig_S1[1, 4],
    colorrange=(-1, 1),
    colormap=colmap_f,
    ticklabelsize=8pt,
)

# --------------------------------------------------------------

ax = Axis(fig_S1[2, 3],
    title=L"$f_{FC}$",
    titlesize=12pt,
)

heatmap!(ax, ϕs, R2ss, f_FC, colormap=colmap_f, colorrange=(-1, 1))
scatter!(ax, xs, ys, markersize=2, color=:red)
scatter!(ax, ϕ_t, R2s_t, markersize=8, color=:lime)
hidedecorations!(ax)

Label(fig_S1[2, 3, TopLeft()], "F",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig_S1[2, 4],
    colorrange=(-1, 1),
    colormap=colmap_f,
    ticklabelsize=8pt,
)

for fig_S in [fig_4, fig_S1]

    # --------------------------------------------------------------

    local ax = Axis(fig_S[3, 1],
        title=L"$\log_{10}\,\left[\,\chi^2(\varphi,\, 0)\,\right]$",
        titlesize=12pt,
        xlabel=L"$\varphi$ [rad]",
        xlabelsize=12pt,
        xticks=([-π, 0.0, π], ["-π", "0", "π"]),
        yticks=([-5:0;]),
        xticklabelsize=8pt,
        yticklabelsize=8pt,
    )

    local rw_in = lines!(ax, ϕs_HR, log10.(res_in_HR.χ2[:RW][:, 1]), color=:red)
    local fc_in = lines!(ax, ϕs_HR, log10.(res_in_HR.χ2[:FC][:, 1]), color=:blue)
    scatter!(ax, -0.55π, -2.5; kwargs_in...)

    xlims!(ax, -π, π)
    ylims!(ax, -5, 0.1)
    axislegend(ax, [fc_in, rw_in], ["FC", "RW"], labelsize=10pt, position=:lb)

    Label(fig_S[3, 1, TopLeft()], "G",
        font=:bold,
        padding=(0, -20, 5, 0),
        halign=:right)

    # --------------------------------------------------------------

    ax = Axis(fig_S[3, 2],
        title=L"$\log_{10}\,\left[\,\chi^2(\varphi,\, 0)\,\right]$",
        titlesize=12pt,
        xlabel=L"$\varphi$ [rad]",
        xlabelsize=12pt,
        xticks=([-π, 0.0, π], ["-π", "0", "π"]),
        xticklabelsize=8pt,
    )

    local rw_op = lines!(ax, ϕs_HR, log10.(res_op_HR.χ2[:RW][:, 1]), color=:red)
    local fc_op = lines!(ax, ϕs_HR, log10.(res_op_HR.χ2[:FC][:, 1]), color=:blue)
    scatter!(ax, -0.55π, -2.5; kwargs_op...)

    xlims!(ax, -π, π)
    ylims!(ax, -5, 0.1)
    hideydecorations!(ax)
    axislegend(ax, [fc_op, rw_op], ["FC", "RW"], labelsize=10pt, position=:lb)

    Label(fig_S[3, 2, TopLeft()], "H",
        font=:bold,
        padding=(0, -20, 5, 0),
        halign=:right)

    # --------------------------------------------------------------

    ax = Axis(fig_S[3, 3],
        title=L"$\log_{10}\,\left[\,\chi^2(\varphi,\, 0)\,\right]$",
        titlesize=12pt,
        xlabel=L"$\varphi$ [rad]",
        xlabelsize=12pt,
        xticks=([-π, 0.0, π], ["-π", "0", "π"]),
        xticklabelsize=8pt,
    )

    rw_op = lines!(ax, ϕs, log10.(res_MC.χ2[:RW][:, 1]), color=:red)
    fc_op = lines!(ax, ϕs, log10.(res_MC.χ2[:FC][:, 1]), color=:blue)

    xlims!(ax, -π, π)
    ylims!(ax, -5, 0.1)
    hideydecorations!(ax)
    axislegend(ax, [fc_op, rw_op], ["FC", "RW"], labelsize=10pt, position=:lb)

    Label(fig_S[3, 3, TopLeft()], "I",
        font=:bold,
        padding=(0, -20, 5, 0),
        halign=:right)
    
    display(fig_S)
end

##

wf_par_2e = map(x -> x[1] * x[2] < 0 ? 1.0 : 0.0, fitpar_RW.c)
wf_par_2e[not_S] .= NaN
t_2e = transpose
r_op_2e, c_op_2e = r_op, c_op
r_in_2e, c_in_2e = r_in, c_in
kwargs_op_2e = kwargs_op
kwargs_in_2e = kwargs_in

##

fig_name = "fig_3"
save(fig_name * ".svg", fig)
save(fig_name * ".eps", fig)
run(`epspdf $fig_name".eps"`)

fig_name = "fig_4"
save(fig_name * ".svg", fig_4)
save(fig_name * ".eps", fig_4)
run(`epspdf $fig_name".eps"`)

##

fig_name = "fig_S1"
save(fig_name * ".svg", fig_S1)
save(fig_name * ".eps", fig_S1)
run(`epspdf $fig_name".eps"`)