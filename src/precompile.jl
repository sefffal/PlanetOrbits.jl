using SnoopPrecompile

@precompile_setup begin
    @precompile_all_calls begin
        o1 = orbit(
           a = 1,
           i = π/4,
           Ω = 0.001,
           ω = 0.001,
           e = 0.5,
           τ = 0.5,
           M = 1,
           tref=0,
       )
       o2 = orbit(
            a = 1,
            i = π/4,
            Ω = 0.001,
            ω = 0.001,
            e = 0.5,
            τ = 0.5,
            M = 1,
            tref=0,
            plx=100.
        )
        o3 = orbit(
            A = 500,
            B = 600,
            F = -500,
            G = 300,
            e = 0.5,
            τ = 0.5,
            M = 1,
            tref=0,
            plx=100.
        )
        o4 = CartesianOrbit(orbitsolve(o1, 0.))
        # CartesianOrbit(orbitsolve(o2, 0.))
        # CartesianOrbit(orbitsolve(o3, 0.))
        for o in (o1,o2,o3,o4)
            hostmass(o)
            period(o)
            meanmotion(o)
            eccentricity(o)
            periastron(o)
            semiamplitude(o)
            radvel(o,0.0)
            posangle(o,0.0)
        end
        for o in (o2,o3)
            distance(o)
            raoff(o, 0.0)
            raoff(o, 0.0)
            decoff(o, 0.0)
            pmra(o, 0.0)
            pmdec(o, 0.0)
        end
    end
end
