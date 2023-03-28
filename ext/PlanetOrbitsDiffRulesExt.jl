module PlanetOrbitsDiffRulesExt

using PlanetOrbits, DiffRules

DiffRules.@define_diffrule PlanetOrbits.kepler_solver(M, e) = :(
    @show M;
        EA = PlanetOrbits.kepler_solver($M,$e);
        temp = 1 - e*cos(EA);
        d_dM = 1 / temp
    ),
    :(
        EA = PlanetOrbits.kepler_solver($M,$e);
        temp = 1 - e*cos(EA);
        d_de = sin(EA) / temp
    )
end


