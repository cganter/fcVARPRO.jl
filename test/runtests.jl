using fcVARPRO
using Test

@testset "fig_1_2.jl" begin
    include("fig_1_2.jl")
end;

@testset "fig_ismrm_challenge_ds_14_sl_3.jl" begin
    include("fig_ismrm_challenge_ds_14_sl_3.jl")
end;

@testset "fig_ismrm_challenge_ds_17_sl_2.jl" begin
    include("fig_ismrm_challenge_ds_17_sl_2.jl")
end;

@testset "fig_cor_two_echoes.jl" begin
    include("fig_cor_two_echoes.jl")
end;
