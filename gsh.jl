# Generalized Swift-Hohenberg model 
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using FFTW
using Plots
using Random
using LaTeXStrings

figdir = joinpath(@__DIR__, "Plots")
mkpath(figdir)

function plot_kwargs()
    fontsizes = (
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        legendfontsize = 14,
    )
    (; linewidth = 3, markersize = 6, markerstrokewidth = 0,
       fontsizes..., size = (600, 500))
end


function generalized_swift_hohenberg_2d(;
    Nx = 256,
    Ny = 256,
    Lx = 100.0,
    Ly = 100.0,
    dt = 1.0e-3,
    T = 100.0,
    epsilon = 0.7,
    q0 = 1.0,
    sigma = 1.0,
    gm = 50.0,
    c2 = 2.0,
    psi0 = nothing,
    seed = 1234,
    snapshot_every = 100,
    store_snapshots = false,
)
    dx = Lx / Nx
    dy = Ly / Ny
    Nt = Int(round(T / dt))

    x = range(-0.5Lx + dx, 0.5Lx, length = Nx)
    y = range(-0.5Ly + dy, 0.5Ly, length = Ny)

    kx = fftfreq(Nx, Nx / Lx) .* 2π
    ky = fftfreq(Ny, Ny / Ly) .* 2π
    KX = kx .+ 0 .* ky' # # (Nx×1) + (1×Ny) → Nx×Ny ->   x-wavenumber for psi_hat[i,j] , KX[i,j] = kx[i], 
    KY = 0 .* kx .+ ky'
    K2 = KX.^2 .+ KY.^2


    # 1/k² with the k = 0 mode set to zero, used for ζ later 
    invK2 = [k2 == 0 ? 0.0 : 1.0 / k2 for k2 in K2]

    # Linear operator for ψ,
    Lpsi = epsilon .- (q0^2 .- K2).^2 # linear part in fourier space
    denom_psi = 1 .- dt .* Lpsi # nonlinear part

    denom_chi = 1 .+ dt .* sigma .* (K2 .+ c2) # denominator of update formula for χ

    # initial condition
    if psi0 === nothing
        Random.seed!(seed)
        psi = 0.05 .* (2 .* rand(Nx, Ny) .- 1)
    else
        psi = [psi0(xi, yi) for xi in x, yi in y]
    end

    #start with fft of ψ
    psi_hat = fft(psi)
    chi_hat = zeros(ComplexF64, Nx, Ny)  # χ = ∇²ζ = 0  →  v = 0

    snapshots = store_snapshots ? Matrix{Float64}[] : nothing # just there to store some things for animation
    times = Float64[]

    # main part
    for step in 1:Nt
        psi = real.(ifft(psi_hat)) # do a inverse fourier transform

        # calculate derivatives of ψ and ∇²ψ
        psi_x = real.(ifft(im .* KX .* psi_hat))
        psi_y = real.(ifft(im .* KY .* psi_hat))
        lap_psi_hat = -K2 .* psi_hat
        lap_psi_x = real.(ifft(im .* KX .* lap_psi_hat)) # ∂x ∇²ψ
        lap_psi_y = real.(ifft(im .* KY .* lap_psi_hat)) # ∂y ∇²ψ

        # Recover ζ from χ = ∇²ζ  because of ζ_hat = -χ_hat / k²
        zeta_hat = -chi_hat .* invK2

        vx = real.(ifft(im .* KY .* zeta_hat)) # vx = ∂y ζ
        vy = real.(ifft(-im .* KX .* zeta_hat))  # vy = -∂x ζ
        ###############################################################
        # calculate nonlinear part in real space and then do fft for ψ and χ equation
        adv = vx .* psi_x .+ vy .* psi_y # v·∇ψ
        nonlinear_psi_hat = fft(-psi.^3 .- adv) # fft of nonlinear part for ψ equation
        
        # Mean-flow forcing: gm * (∂x ∇²ψ ∂yψ - ∂y ∇²ψ ∂xψ)
        forcing_hat = fft(gm .* (lap_psi_x .* psi_y .- lap_psi_y .* psi_x)) # gm·F for χ equation
        
        ###############################################################
        # use semi-implicit update formulas
        psi_hat = (psi_hat .+ dt .* nonlinear_psi_hat) ./ denom_psi
        chi_hat = (chi_hat .+ dt .* forcing_hat) ./ denom_chi
        chi_hat[1, 1] = 0.0 + 0.0im

        ##################
        # just some stuff for animating
        if store_snapshots && (step % snapshot_every == 0 || step == Nt)
            push!(snapshots, real.(ifft(psi_hat)))
            push!(times, step * dt)
        end
        ##################
    end
    # now calculate ψ and the rest
    psi = real.(ifft(psi_hat))
    zeta_hat = -chi_hat .* invK2
    vx = real.(ifft(im .* KY .* zeta_hat))
    vy = real.(ifft(-im .* KX .* zeta_hat))

    return (; x, y, psi, psi_hat, vx, vy, zeta_hat, chi_hat,
            epsilon, q0, sigma, gm, c2, T, dt, snapshots, times)
end

function plot_gsh_final(data; filename = joinpath(figdir, "gsh_meanflow_final.png"))
    p1 = heatmap(data.x, data.y, data.psi',
        aspect_ratio = 1,
        xlims = extrema(data.x), ylims = extrema(data.y),
        xlabel = "x", ylabel = "y",
        title = "ψ at t = $(round(data.T, digits = 2))",
        colorbar_title = "ψ",
        c = :viridis,
        framestyle = :box,
        right_margin = 5Plots.mm,
        titlefontsize = 16,
        xtickfontsize = 12, ytickfontsize = 12,
        xguidefontsize = 14, yguidefontsize = 14,
        size = (650, 600),
    )
    savefig(p1, filename)
    return p1
end

function generalized_swift_hohenberg_2d_anim(;
     Nx = 2*256,
    Ny = 2*256,
    Lx = 2*256 * 2π / 16,
    Ly = 2*256 * 2π / 16,
    dt = 1e-1,
    T = 100.0,
    epsilon = 0.7,
    q0 = 1.0,
    sigma = 1.0,
    gm = 50.0,
    c2 = 2.0,
    nframes = 250,
    fps = 20,
    filename = joinpath(figdir, "gsh_meanflow.gif"),
)
    Nt = Int(round(T / dt))
    save_every = max(1, Nt ÷ nframes)

    data = generalized_swift_hohenberg_2d(
        Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, dt = dt, T = T,
        epsilon = epsilon, q0 = q0, sigma = sigma, gm = gm, c2 = c2,
        snapshot_every = save_every, store_snapshots = true,
    )

    anim = Animation()
    for (ψ, t) in zip(data.snapshots, data.times)
        plt = heatmap(data.x, data.y, ψ', aspect_ratio = 1,
            xlabel = "x", ylabel = "y",
            title = "Generalized Swift-Hohenberg: t = $(round(t, digits = 2))",
            size = (800, 700), c = :viridis)
        frame(anim, plt)
    end
    gif(anim, filename, fps = fps)
    return filename
end


data = generalized_swift_hohenberg_2d(
    Nx = 2*256,
    Ny = 2*256,
    Lx = 2*256 * 2π / 16,
    Ly = 2*256 * 2π / 16,
    dt = 1e-1,
    T = 100.0,
    epsilon = 0.7,
    q0 = 1.0,
    sigma = 1.0,
    gm = 50.0,
    c2 = 2.0,
)
plot_gsh_final(data)
 generalized_swift_hohenberg_2d_anim(T = 100.0, nframes = 150)