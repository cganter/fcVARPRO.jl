using Revise, MAT, CairoMakie, LinearAlgebra
import VP4Optim as VP
import B0Map as BM
import fcVARPRO as FCV

BLAS.set_num_threads(1)

# coronal three-echo dataset

# threshold for mask
sl = 37
R2s_max = 5.0
optim = true
#thresh, data_set = 25, "151725_0302"
thresh, data_set = 70, "171032_0302"
#thresh, data_set = 60, "183316_0202"

# local fit
fitopt_FW = BM.fitOpt()
fitopt_FW.R2s_rng = [0.0, R2s_max]
fitopt_FW.optim = optim
@time res_FW = FCV.three_echoes(BM.GREMultiEchoWFFW, fitopt_FW; data_set=data_set, thresh=thresh, slice=sl);

fitopt_RW = BM.fitOpt()
fitopt_RW.R2s_rng = [0.0, R2s_max]
fitopt_RW.optim = optim
@time res_RW = FCV.three_echoes(BM.GREMultiEchoWFRW, fitopt_RW; data_set=data_set, thresh=thresh, slice=sl);

fitopt_FC = BM.fitOpt()
fitopt_FC.R2s_rng = [0.0, R2s_max]
fitopt_FC.optim = optim
@time res_FC = FCV.three_echoes(BM.GREMultiEchoWF, fitopt_FC; data_set=data_set, thresh=thresh, slice=sl);

##

not_S = (!).(res_FW.fitpar.S)

# extract maps
ϕ_FW = res_FW.fitpar.ϕ#[:,:,sl]
ϕ_FW[not_S] .= NaN
f_FW = BM.fat_fraction_map(res_FW.fitpar, fitopt_FW)#[:,:,sl]
f_FW[not_S] .= NaN
T2s_FW = 1 ./ res_FW.fitpar.R2s#[:,:,sl]

ϕ_RW = res_RW.fitpar.ϕ#[:,:,sl]
ϕ_RW[not_S] .= NaN
f_RW = BM.fat_fraction_map(res_RW.fitpar, fitopt_RW)#[:,:,sl]
f_RW[not_S] .= NaN
T2s_RW = 1 ./ res_RW.fitpar.R2s#[:,:,sl]

ϕ_FC = res_FC.fitpar.ϕ#[:,:,sl]
ϕ_FC[not_S] .= NaN
f_FC = BM.fat_fraction_map(res_FC.fitpar, fitopt_FC)#[:,:,sl]
f_FC[not_S] .= NaN
T2s_FC = 1 ./ res_FC.fitpar.R2s#[:,:,sl]

T2s_FW[not_S] .= NaN
T2s_RW[not_S] .= NaN
T2s_FC[not_S] .= NaN

## =============== generate phase maps =======================

r, c = 55, 135
kwargs = (; markersize = 10, marker = :circle, strokewidth = 2, color = (:red, 0.0), strokecolor = :red)
T2s_min, T2s_max = 0.0, 100.0

# size definitions
# these are relative to 1 CSS px
inch = 96
pt = 4 / 3

width, height = 6.9inch, 6.9inch
colmapO = :romaO

fig = Figure(size=(width, height), fontsize=12pt)

ax = Axis(fig[1, 1],
    title=L"$\varphi_{FW}$",
)

heatmap!(ax,
    rotr90(ϕ_FW),
    colorrange=(-π, π),
    colormap=colmapO,
    nan_color=:black,
)
hidedecorations!(ax)
Label(fig[1, 1, TopLeft()], "A",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[1, 2],
    title=L"$\varphi_{RW}$",
)

heatmap!(ax,
    rotr90(ϕ_RW),
    colorrange=(-π, π),
    colormap=colmapO,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r, c, ϕ_RW, rotr90)...; kwargs...)
hidedecorations!(ax)
Label(fig[1, 2, TopLeft()], "B",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[1, 3],
    title=L"$\varphi_{FC}$",
)

heatmap!(ax,
    rotr90(ϕ_FC),
    colorrange=(-π, π),
    colormap=colmapO,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r, c, ϕ_FC, rotr90)...; kwargs...)
hidedecorations!(ax)
Label(fig[1, 3, TopLeft()], "C",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig[1, 4],
    colorrange=(-π, π),
    colormap=colmapO,
    ticks=([-π, 0.0, π], ["-π", "0", "π"]),
    ticklabelsize=8pt,
)

# =============== generate fat fraction maps =======================

colmap = :imola

ax = Axis(fig[2, 1],
    title=L"$f_{FW}$",
)

heatmap!(ax,
    rotr90(f_FW),
    colorrange=(0, 1),
    colormap=colmap,
    nan_color=:black,
)
hidedecorations!(ax)
Label(fig[2, 1, TopLeft()], "D",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[2, 2],
    title=L"$f_{RW}$",
)

heatmap!(ax,
    rotr90(f_RW),
    colorrange=(0, 1),
    colormap=colmap,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r, c, f_RW, rotr90)...; kwargs...)
hidedecorations!(ax)
Label(fig[2, 2, TopLeft()], "E",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[2, 3],
    title=L"$f_{FC}$",
)

heatmap!(ax,
    rotr90(f_FC),
    colorrange=(0, 1),
    colormap=colmap,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r, c, f_RW, rotr90)...; kwargs...)
hidedecorations!(ax)
Label(fig[2, 3, TopLeft()], "F",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig[2, 4],
    colorrange=(0, 1),
    colormap=colmap,
    ticklabelsize=8pt,
)

# =============== generate R2* maps =======================

colmap = :batlow

ax = Axis(fig[3, 1],
    title=L"$T^\ast_{2,FW}$",
)

heatmap!(ax,
    rotr90(T2s_FW),
    colorrange=(T2s_min, T2s_max),
    colormap=colmap,
    nan_color=:black,
)
hidedecorations!(ax)
Label(fig[3, 1, TopLeft()], "G",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[3, 2],
    title=L"$T^\ast_{2,RW}$",
)

heatmap!(ax,
    rotr90(T2s_RW),
    colorrange=(T2s_min, T2s_max),
    colormap=colmap,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r, c, T2s_RW, rotr90)...; kwargs...)
hidedecorations!(ax)
Label(fig[3, 2, TopLeft()], "H",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[3, 3],
    title=L"$T^\ast_{2,FC}$",
)

heatmap!(ax,
    rotr90(T2s_FC),
    colorrange=(T2s_min, T2s_max),
    colormap=colmap,
    nan_color=:black,
)
scatter!(ax, FCV.rc2xy(r, c, T2s_FC, rotr90)...; kwargs...)

hidedecorations!(ax)
Label(fig[3, 3, TopLeft()], "I",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

Colorbar(fig[3, 4],
    colorrange=(T2s_min, T2s_max),
    colormap=colmap,
    ticklabelsize=8pt,
)

display(fig)

##

wf_par_3e = map(x -> x[1] * x[2] < 0 ? 1.0 : 0.0, res_RW.fitpar.c)
wf_par_3e[not_S] .= NaN
t_3e = rotr90
r_3e, c_3e = r, c
kwargs_3e = kwargs

##

fig_name = "fig_2"
save(fig_name * ".svg", fig)
save(fig_name * ".eps", fig)
run(`epspdf $fig_name".eps"`)