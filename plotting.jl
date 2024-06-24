using Plots
using LaTeXStrings

include("bounds.jl")

const colors = palette(:default)
default(lw=2)

bound_with_attr_list = [
    (exact_value, "Exact", colors[1]),
    (hoeffding_bound, "Hoeffding", colors[2]),
    (bernstein_bound, "Bernstein", colors[3]),
    (bennett_bound, "Bennett", colors[4]),
    (sanov_bound, "Sanov", colors[5]),
]

N_with_style_list = [
    (100, :solid),
    (1000, :dash),
    (10000, :dot),
]
p = 0.1
nϵ = 100
ϵs = range(0, 0.2, length=nϵ)

plt = plot(; xlabel=L"\epsilon", ylabel="log bound")

for bound_with_attr in bound_with_attr_list
    bound, _, color = bound_with_attr
    for N_with_style in N_with_style_list
        N, ls = N_with_style
        bounds = [bound(N, p, ϵ) for ϵ in ϵs]
        plot!(plt, ϵs, bounds, label=false, c=color, ls=ls)
    end
end

for bound_with_attr in bound_with_attr_list
    bound, name, color = bound_with_attr
    plot!(plt, Shape(Int[], Int[]), label=name, c=color, lw=0)
end

for N_with_style in N_with_style_list
    N, ls = N_with_style
    plot!(plt, [], [], label=L"N=%$(N)", c=:black, ls=ls)
end

savefig(plt, "comparison_bounds.png")