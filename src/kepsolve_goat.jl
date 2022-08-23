using StaticArrays


"""
    PlanetOrbits.Goat()

Kepler solver implementation from https://arxiv.org/abs/2103.15829 and https://github.com/oliverphilcox/Keplers-Goat-Herd

It is here for comparison purposes only. In general, Markley() is more performant and accurate.
"""
struct Goat <: AbstractSolver end

# Implementation from https://arxiv.org/abs/2103.15829 and https://github.com/oliverphilcox/Keplers-Goat-Herd
function kepler_solver(𝓁, e, ::Goat)

    # This function implements the 🐐 GOAT algorithm for 
    # solving Kepler's equation. It is approximately
    # 4x faster than the other methods implemented
    # here.

    if isapprox(e, 0)
        return 𝓁
    end

    if isapprox(rem(𝓁,π), 0)
        return 𝓁
    end

    N_it = 15
    N_points = N_it-2
    N_fft = (N_it)*2

    radius = e / 2

    # Generate e^{ikx} sampling points and precompute real and imaginary parts
    # Keep these on the stack inside an MVector -> no allocations
    exp2R = @MVector zeros(typeof(e), N_points)
    exp2I = @MVector zeros(typeof(e), N_points)
    exp4R = @MVector zeros(typeof(e), N_points)
    exp4I = @MVector zeros(typeof(e), N_points)
    coshI = @MVector zeros(typeof(e), N_points)
    sinhI = @MVector zeros(typeof(e), N_points)
    ecosR = @MVector zeros(typeof(e), N_points)
    esinR = @MVector zeros(typeof(e), N_points)
    @inbounds for j in 1:N_points
        freq = 2π*j/N_fft
        cf = cos(freq)
        sf = sin(freq)
        exp2R[j] = cf
        exp2I[j] = sf
        exp4R[j] = cf*cf-sf*sf
        exp4I[j] = 2.0*cf*sf
        coshI[j] = cosh(radius*exp2I[j])
        sinhI[j] = sinh(radius*exp2I[j])
        ecosR[j] = e*cos(radius*exp2R[j])
        esinR[j] = e*sin(radius*exp2R[j])
    end

    esinRadius = e*sin(radius)
    ecosRadius = e*cos(radius)


    # Define contour center for each ell and precompute sin(center), cos(center)
    if 𝓁 < π
        center = 𝓁 + e/2
    else
        center = 𝓁 - e/2
    end
    sinC = sin(center)
    cosC = cos(center)
    output = center

    # Accumulate Fourier coefficients
    # NB: we halve the range by symmetry, absorbing factor of 2 into ratio

    #######
    # Separate out j = 0 piece, which is simpler

    # Compute z in real and imaginary parts (zI = 0 here)
    zR = center + radius

    # Compute e*sin(zR) from precomputed quantities
    tmpsin = sinC*ecosRadius+cosC*esinRadius # sin(zR)

    # Compute f(z(x)) in real and imaginary parts (fxI = 0)
    fxR = zR - tmpsin - 𝓁

    # Add to array, with factor of 1/2 since an edge
    ft_gx2 = 0.5/fxR
    ft_gx1 = 0.5/fxR

    #######
    # Compute for j = 1 to N_points
    @inbounds @simd for j in 1:N_points

        # Compute z in real and imaginary parts
        zR = center + radius*exp2R[j]
        zI = radius*exp2I[j]

        # Compute f(z(x)) in real and imaginary parts
        # can use precomputed cosh / sinh / cos / sin for this!
        tmpcosh = coshI[j] # cosh(zI)
        tmpsinh = sinhI[j] # sinh(zI)
        tmpsin = sinC*ecosR[j]+cosC*esinR[j] # e sin(zR)
        tmpcos = cosC*ecosR[j]-sinC*esinR[j] # e cos(zR)

        fxR = zR - tmpsin*tmpcosh-𝓁
        fxI = zI - tmpcos*tmpsinh

        # Compute 1/f(z) and append to array
        ftmp = fxR*fxR+fxI*fxI
        fxR /= ftmp
        fxI /= ftmp

        ft_gx2 += (exp4R[j]*fxR+exp4I[j]*fxI)
        ft_gx1 += (exp2R[j]*fxR+exp2I[j]*fxI)
    end

    #######
    # Separate out j = N_it piece, which is simpler

    # Compute z in real and imaginary parts (zI = 0 here)
    zR = center - radius

    # Compute sin(zR) from precomputed quantities
    tmpsin = sinC*ecosRadius-cosC*esinRadius # sin(zR)

    # Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
    fxR = zR - tmpsin-𝓁

    # Add to sum, with 1/2 factor for edges
    ft_gx2 += 0.5/fxR
    ft_gx1 += -0.5/fxR

    #######
    # Compute E(ell)
    output += radius*ft_gx2/ft_gx1;
    
    return output
end

