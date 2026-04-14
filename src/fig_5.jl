using Revise, CairoMakie, LinearAlgebra
import fcVARPRO as FCV

# size definitions
# these are relative to 1 CSS px
inch = 96
pt = 4 / 3

width, height = 6.9inch, 2.5inch
colmap = :imola

fig = Figure(size=(width, height), fontsize=12pt)

ax = Axis(fig[1, 1],
    title=L"4 Echoes$$",
)

heatmap!(ax,
    t_4e(wf_par_4e),
    colormap=colmap,
    nan_color=:black,
)

hidedecorations!(ax)

Label(fig[1, 1, TopLeft()], "A",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[1, 2],
    title=L"3 Echoes$$",
)

heatmap!(ax,
    t_3e(wf_par_3e),
    colormap=colmap,
    nan_color=:black,
)

scatter!(ax, FCV.rc2xy(r_3e, c_3e, wf_par_3e, t_3e)...; kwargs_3e...)

hidedecorations!(ax)

Label(fig[1, 2, TopLeft()], "B",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

ax = Axis(fig[1, 3],
    title=L"2 Echoes$$",
)

heatmap!(ax,
    t_2e(wf_par_2e),
    colormap=colmap,
    nan_color=:black,
)


hidedecorations!(ax)

Label(fig[1, 3, TopLeft()], "C",
    font=:bold,
    padding=(0, -20, 5, 0),
    halign=:right)

scatter!(ax, FCV.rc2xy(r_op_2e, c_op_2e, wf_par_2e, t_2e)...; kwargs_op_2e...)
scatter!(ax, FCV.rc2xy(r_in_2e, c_in_2e, wf_par_2e, t_2e)...; kwargs_in_2e...)

display(fig)

##

fig_name = "fig_5"
save(fig_name * ".svg", fig)
save(fig_name * ".eps", fig)
run(`epspdf $fig_name".eps"`)