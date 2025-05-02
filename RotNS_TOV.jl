import OrdinaryDiffEq

export TOVsolution
export TOVSchwarzschild!, TOVSchwarzschild
export TOVin!, integrateTOVin
export TOVisotropic!, TOVisotropic, TOVshootout
export TOVModel

#########################################################################

# Common TOV solution struct

"""
struct TOVsolution
    TotalMass       : dimensionless total gravitational energy
    CentralDensity  : dimensionless total central energy density
    CentralPressure : dimensionless central pressure
    SurfaceRadius   : dimneionless coordinate radius (gauge dependent)
    Solution        : Full solution returned by integrator
end
"""
struct TOVsolution
    TotalMass::Float64
    CentralDensity::Float64
    CentralPressure::Float64
    SurfacePressure::Float64
    SurfaceRadius::Float64
    Solution
end



#########################################################################
# Methods for integration the TOV equations in Schwarzschild coordinates.

"""
    TOVSchwarzschild!(du,u,p,r)

Return in `du` the vector of derivatives for the TOV equations in terms of the Schwarzschild radius ``\\tilde{r}``.  
The Schwarzschild TOV variables are:
- `u[1]` : ``m(\\tilde{r})``, the metric mass function.
- `u[2]` : ``P(\\tilde{r})``, the dimensionless pressure.
- `u[3]` : ``\\Phi(\\tilde{r})``, the metric function.

In Geometrized units, the metric in Schwarzschild radial gauge, and the TOV equations are:
```math
\\begin{align}
	\\rm{d}s^2 &= -e^{2\\Phi(\\tilde{r})}{\\rm d}t^2 + \\frac1{1-\\frac{2m(\\tilde{r})}{\\tilde{r}}}{\\rm d}\\tilde{r}^2+\\tilde{r}^2{\\rm d}^2\\Omega\\qquad\\text{for}\\quad 0\\le \\tilde{r}\\le \\tilde{R}, \\\\
		&= -\\left(1-\\frac{2M}{\\tilde{r}}\\right){\\rm d}t^2 + \\frac1{\\left(1-\\frac{2M}{\\tilde{r}}\\right)}{\\rm d}\\tilde{r}^2+\\tilde{r}^2{\\rm d}^2\\Omega\\qquad\\text{for}\\quad \\tilde{R}\\le \\tilde{r},
\\end{align}
```
```math
\\begin{align}
	\\frac{{\\rm d} m(\\tilde{r})}{{\\rm d}\\tilde{r}} &= 4\\pi \\tilde{r}^2\\epsilon(\\tilde{r}) &&\\text{with}\\qquad m(\\tilde{r}=0)=0,\\\\
	\\frac{{\\rm d} P(\\tilde{r})}{{\\rm d}\\tilde{r}} &= -\\frac{[\\epsilon(\\tilde{r})+P(\\tilde{r})][m(\\tilde{r})+4\\pi \\tilde{r}^3P(\\tilde{r})]}{\\tilde{r}[\\tilde{r}-2m(\\tilde{r})]} &&\\text{with}\\qquad\\left\\{\\begin{array}{l}P(\\tilde{r}=0)=P_c,\\\\
	P(\\tilde{r}=\\tilde{R})=0,\\end{array}\\right.\\\\
	\\frac{{\\rm d} \\Phi(\\tilde{r})}{{\\rm d}\\tilde{r}} &= \\frac{m(\\tilde{r})+4\\pi \\tilde{r}^3P(\\tilde{r})}{\\tilde{r}[\\tilde{r}-2m(\\tilde{r})]} &&\\text{with}\\qquad \\Phi(\\tilde{r}=\\tilde{R})=\\textstyle\\frac12\\ln\\left(1-\\frac{2M}{\\tilde{R}}\\right).
\\end{align}
```

`u` contains the initial values, and both `u` and `du` store values in the same order. 
`p` is a vector containing parameters needed during the integration.  The single parameter passed in through `p` is:
- p[1] : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
"""

function TOVSchwarzschild!(du,u,p,r)
    eos = p[1]
    m = u[1]
    P = max(0,u[2]) # max needed to prevent negative P in EOS conversion
    Φ = u[3]
    ed = EnergyDensityfromPressure(eos,P)
    fourpi = 4*pi
    du[1] = fourpi*r^2*ed
    du[3] = r==0 ? 0 : (m + fourpi*r^3*P)/(r*(r - 2*m))
    du[2] = -(ed + P)*du[3]
end

"""
TOVSchwarzschild(eos, Ec, MaxR; 
                 SurfacePressure=0.0, AbsTol=1.0e-10, RelTol=1.0e-10, InitialStep=1.0e-4)

Integrate the TOV equations in Schwarzschild radial coordinates using `TOVSchwarzschild!`, starting at ``\\tilde{r}=0``,
and finding the radius of the star as defined by the condition that the dimensionless pressure ``\\bar{P}`` equal the value 
specified `SurfacePressure`.  Integration stops when ``\tilde{r}`` reaches the surface radius `MaxR`, or when the pressure 
falls below `SurfacePressure`, whichever occurs first.
- `EOS`   : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
- `Ec`    : The dimensionless total energy density at ``\tilde{r}=0``.
- `MaxR`  : A dimensionless radial integration limit that is larger than the Schwarzschild radius of the star.

The solution is returned via the `TOVsolution` object.

##### Optional keyword arguements:
- `SurfacePressure` : The desired dimensionless pressure at the surface.
- `AbsTol`          : The absolute tolerance used during the integration.
- `RelTol`          : The relative tolerance used during the integration.
- `InitialStep`     : The initial step size for the integration relative to the surface radius.
"""
function TOVSchwarzschild(eos::EOS,Ec::Number,MaxR::Number; SurfacePressure::Number=0.0,
						  AbsTol::Number=1.0e-10, RelTol::Number=1.0e-10, InitialStep::Number=1.0e-4)
    u0 = [0, PressurefromEnergyDensity(eos,Ec), 0]            # initial value m(0), P(0), Φ(0)
    rspan = (0.0,MaxR)                                        # Maximum range of integration.  MaxR must be larger than radius of NS
                                                              # Set up functions to allow solver to find surface where P(R)=SurfacePressure
    findzeroP(u,r,integrator) = u[2]-SurfacePressure          # Callback to find the surface with pressure SurfaceP
                                                              # Set up the TOV problem
    TOVoutprob = OrdinaryDiffEq.ODEProblem(TOVSchwarzschild!,u0,rspan,[eos]) # [eos] is the list of parameters passed into the integrator
    sol = OrdinaryDiffEq.solve(TOVoutprob,OrdinaryDiffEq.Vern9(), dt=InitialStep*MaxR,      # Solve the system (tried DP5, BS5, Vern9)
							   abstol=AbsTol,reltol=RelTol,
    						   callback=OrdinaryDiffEq.ContinuousCallback(findzeroP,OrdinaryDiffEq.terminate!))
    return TOVsolution(sol[1,end],Ec,u0[2],sol[2,end],sol.t[end],sol)
end



#########################################################################
# Methods for integration the TOV equations in isotropic coordinates.

# define the system
"""
    TOVin!(du,u,p,r)

Return in `du` the vector of derivatives for the TOV equations in terms of the isotropic radius `r`.  The isotropic TOV variables are:
- `u[1]` : ``\\tilde{r}(r)``, the Schwarzschild radial coordinate.
- `u[2]` : ``\\beta(r)``, the metric function.
- `u[3]` : ``P(r)``, the dimensionless pressure.
- `u[4]` : ``m(r)``, the metric mass function.

In Geometrized units, the metric in isotropic radial gauge, and the transformed TOV equations are:
```math
\\begin{align}
    \\rm{d}s^2 &= -e^{2\\beta(r)}{\\rm d}t^2 + e^{2\\alpha(r)}({\\rm d}r^2+r^2{\\rm d}^2\\Omega)
    \\qquad\\text{with}\\qquad e^{\\alpha(r)}\\equiv\\frac{\\tilde{r}(r)}{r}\\quad\\text{and}\\quad\\beta(r)\\equiv\\Phi(\\tilde{r}(r))
\\end{align}
```
```math
\\begin{align}
    R = \\frac{\\tilde{R}}4\\left(1 + \\sqrt{1 - \\frac{2M}{\\tilde{R}}}\\right)^2
\\end{align}
```
```math
\\begin{align}
    \\text{for}\\quad r\\ge R\\quad:\\quad
    \\tilde{r}(r) = r\\left(1 + \\frac{M}{2r}\\right)^2 \\quad;\\quad
    e^{2\\alpha(r)} = \\left(1 + \\frac{M}{2r}\\right)^4 \\quad;\\quad
    e^{2\\beta(r)} = \\left(\\frac{1-\\frac{M}{2r}}{1+\\frac{M}{2r}}\\right)^2 
\\end{align}
```
```math
\\begin{align}
    \\frac{{\\rm d}\\tilde{r}(r)}{{\\rm d}r} &= \\frac{\\tilde{r}}{r}\\sqrt{1 - \\frac{2m(r)}{\\tilde{r}}} &&\\text{with}\\qquad \\tilde{r}(R)=\\tilde{R}, \\\\
    \\frac{{\\rm d}\\beta(r)}{{\\rm d}r} &= \\frac{m(r)+4\\pi \\tilde{r}^3P(r)}{\\tilde{r}r\\sqrt{1 - \\frac{2m(r)}{\\tilde{r}}}} &&\\text{with}\\qquad \\beta(R)=\\ln\\left(\\frac{1-\\tfrac{M}{2R}}{1+\\tfrac{M}{2R}}\\right), \\\\
    \\frac{{\\rm d}P(r)}{{\\rm d}r} &= -[\\epsilon(r)+P(r)]\\frac{{\\rm d}\\beta(r)}{{\\rm d}r} &&\\text{with}\\qquad P(R)=0, \\\\
    \\frac{{\\rm d}m(r)}{{\\rm d}r} &= 4\\pi\\epsilon(r)\\frac{\\tilde{r}^3}{r}\\sqrt{1 - \\frac{2m(r)}{\\tilde{r}}} &&\\text{with}\\qquad m(R)=M.
\\end{align}
```

`u` contains the initial values, and both `u` and `du` store values in the same order. 
`p` is a vector containing parameters needed during the integration.  The single parameter passed in through `p` is:
- p[1] : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.

This version of the isotropic coordinate derivatives is intended for inward integrations starting from surface of the star.
Inward integration becomes very inaccurate as ``r=0`` is approached.
"""
function TOVin!(du,u,p,r)
    eos = p[1]
    r_tilde = u[1]  # Schwarzschild radius
    beta = u[2]  # metric function
    P = u[3] # pressure
    M = max(0, u[4])  # mass


    if r == 0.0 || P <= 0.0   
        du[1] = 0.0
        du[2] = 0.0
        du[3] = 0.0
        du[4] = 0.0
    else
        fourpi = 4*pi
        a = sqrt(max(0,1 - (2*M / r_tilde))) # avoid sqrt(-1) errors.
        ed = EnergyDensityfromPressure(eos,P) 
        du[1] = (r_tilde / r) * a
        du[2] = (M + fourpi * r_tilde^3 * P) / (r^2 * du[1])
        du[3] = -(ed + P) * du[2]
        du[4] = fourpi*ed * r_tilde^2*du[1]
    end

end

"""
    misbehave(u, r, integrator)

Callback function used by integrateTOVin to stop the integration when it starts to become inaccurate as it approaches ``r = 0``. 
Closeness to the origin is determined by some fraction of the starting radius of the integration.

There are several criteria used to determin when the integration is becoming inaccurate:
1) The mass becomes non-positive.
2) ``2m(r) > \\tilde{r}`` : This causes the solution to become complex.
3) ``d\\tilde{r}/dr`` becomes too small.  Currently too small is hard coded to be 0.5.
4) ``d\\sqrt{1 - 2m(r)/\\tilde{r}}/dr`` becomes too large when sufficiently close to ``r=0``.  
    Currently too large is hard coded to be 1.0, and sufficiently close is hard coded to be ``r<0.05R``
    (where ``R`` is the surface radius of the star.)
5) Sufficiently close to ``r=0``, leading order power-law behavior of the mass function deviates too far from ``m(r)\\propto r^3``.
    Given an estimated power-law behavior ``r^n``, the behavior is hard coded to be out of range if ``n<2.5`` or ``n>10``.  
    Sufficiently close is hard coded to be ``r<0.02R`` (where ``R`` is the surface radius of the star.)
"""
function misbehave(u, r, integrator)
    eos = integrator.p[1] # Get the EOS for this integration
    R = integrator.p[2] # Get the starting radius for this integration
    r_tilde = u[1]  # areal radius
    β = u[2]  # metric function
    P = u[3] # Pressure
    M = max(0, u[4])  # mass
    
    if (u[4] <= 0)
        @debug "Inward integration terminated because M <=0"
        return true
    end
    
    a = sqrt(max(0,1 - (2*M / r_tilde))) # Terminate on sqrt(-1) error
    drtdr = (r_tilde/r)*a
    if drtdr < 0.5
        if 2*M > r_tilde
            @debug "Inward TOV ntegration terminated because 2*M > r_tilde"
        else
            @debug "Inward TOV integration terminated because d(r_tilde)/dr < 0.5"
        end
        return true
    end

    ed = EnergyDensityfromPressure(eos,P)

    if r < 0.05*R
        dsqrtdr = M/(r_tilde*r)-4*pi*ed*r_tilde^2/r
        if dsqrtdr > 1 
            @debug "Inward TOV integration terminated because d(sqrt(1-2M/r_tilde))/dr > 1"
            return true
        end
    end
    
    if r < 0.02*R
        slope = 4*pi * ed * r_tilde^3 * a/M
        if slope < 2.5 || slope > 10
            @debug "Inward TOV integration terminated because M ~ r^3 is being violated"
            return true
        else
            return false # Don't terminate
        end
    else
        return false # Don't terminate
    end
end

"""
    integrateTOVin(EOS, solution; AbsTol=1.0e-10, RelTol=1.0e-10, InitialStep=1.0e-5)

Integrate the TOV equations in isotropic radial coordinates using `TOVin!`, starting from the surface of the star and integrating 
toward``r=0``.  This integration is used to determine an initial guess for an initial guess for ``d\\tilde{r}/dr``.  Because 
inward integration becomes inaccurate as ``r=0`` is approached, a callback function is used to decide when the integration is
not behaving well, at which point the integration is stopped.
Integration stops when ``r`` reaches the surface radius `Rsurf`, or when the pressure falls below `SurfacePressure`,
whichever occurs first.
- `EOS`      : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
- `solution` : The dimensionless total energy density at ``r=0``.

The return value is the tuple `(drtdr0, sol)`, where
- `drtdr0` is the value of ``d\\tilde{r}/dr`` at ``r=0``, where ``\\tilde{r}`` is the Schwarzschild radial coordinate.
- `sol`  is the full solution object returned by the `OrdinaryDiffEq` function `solve`

##### Optional keyword arguements:
- `AbsTol`          : The absolute tolerance used during the integration.
- `RelTol`          : The relative tolerance used during the integration.
- `InitialStep`     : The initial step size for the integration relative to the surface radius.
"""
function integrateTOVin(EOS::EOS, solution::TOVsolution; 
                        AbsTol::Number=1.0e-10, RelTol::Number=1.0e-10, InitialStep::Number=1.0e-5)
    cd = solution.CentralDensity
    R_tilde = solution.SurfaceRadius # Schwarzschild radius of the surface from outward TOV integration
    M = solution.TotalMass # mass from outward TOV integration
    P = solution.SurfacePressure # surface pressure from outward TOV integration
     
    # Calculate initial conditions
    R = (R_tilde/4)*(1+ sqrt(1 - (2*M / R_tilde)))^2 # isotropic surface radius
    β_R = log((1 - M / (2 * R)) / (1 + M / (2 * R))) # surface value for metric function β

    u0 = [R_tilde, β_R, P, M]  # Initial conditions for the inward integration
    r_span = (R, 5.e-3*R)  # Time span from the outer metric radius inward - longer than where instability exists
    prob = OrdinaryDiffEq.ODEProblem(TOVin!, u0, r_span, (EOS,R))
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Vern9(),
    						   dt=InitialStep*R, abstol=AbsTol,reltol=RelTol,
    						   callback=OrdinaryDiffEq.DiscreteCallback(misbehave,OrdinaryDiffEq.terminate!))
    drtdr0 = sol[1,end]/sol.t[end]*sqrt(1-2*sol[4,end]/sol[1,end])
    return (drtdr0,sol)
end

"""
    TOVisotropic!(du,u,p,r)

Return in `du` the vector of derivatives for the TOV equations in terms of the isotropic radius `r`.  The isotropic TOV variables are:
- `u[1]` : ``\\tilde{r}(r)``, the Schwarzschild radial coordinate.
- `u[2]` : ``\\beta(r)``, the metric function.
- `u[3]` : ``P(r)``, the dimensionless pressure.
- `u[4]` : ``m(r)``, the metric mass function.

In Geometrized units, the metric in isotropic radial gauge, and the transformed TOV equations are:
```math
\\begin{align}
    \\rm{d}s^2 &= -e^{2\\beta(r)}{\\rm d}t^2 + e^{2\\alpha(r)}({\\rm d}r^2+r^2{\\rm d}^2\\Omega)
    \\qquad\\text{with}\\qquad e^{\\alpha(r)}\\equiv\\frac{\\tilde{r}(r)}{r}\\quad\\text{and}\\quad\\beta(r)\\equiv\\Phi(\\tilde{r}(r))
\\end{align}
```
```math
\\begin{align}
    R = \\frac{\\tilde{R}}4\\left(1 + \\sqrt{1 - \\frac{2M}{\\tilde{R}}}\\right)^2
\\end{align}
```
```math
\\begin{align}
    \\text{for}\\quad r\\ge R\\quad:\\quad
    \\tilde{r}(r) = r\\left(1 + \\frac{M}{2r}\\right)^2 \\quad;\\quad
    e^{2\\alpha(r)} = \\left(1 + \\frac{M}{2r}\\right)^4 \\quad;\\quad
    e^{2\\beta(r)} = \\left(\\frac{1-\\frac{M}{2r}}{1+\\frac{M}{2r}}\\right)^2 
\\end{align}
```
```math
\\begin{align}
    \\frac{{\\rm d}\\tilde{r}(r)}{{\\rm d}r} &= \\frac{\\tilde{r}}{r}\\sqrt{1 - \\frac{2m(r)}{\\tilde{r}}} &&\\text{with}\\qquad \\tilde{r}(R)=\\tilde{R}, \\\\
    \\frac{{\\rm d}\\beta(r)}{{\\rm d}r} &= \\frac{m(r)+4\\pi \\tilde{r}^3P(r)}{\\tilde{r}r\\sqrt{1 - \\frac{2m(r)}{\\tilde{r}}}} &&\\text{with}\\qquad \\beta(R)=\\ln\\left(\\frac{1-\\tfrac{M}{2R}}{1+\\tfrac{M}{2R}}\\right), \\\\
    \\frac{{\\rm d}P(r)}{{\\rm d}r} &= -[\\epsilon(r)+P(r)]\\frac{{\\rm d}\\beta(r)}{{\\rm d}r} &&\\text{with}\\qquad P(R)=0, \\\\
    \\frac{{\\rm d}m(r)}{{\\rm d}r} &= 4\\pi\\epsilon(r)\\frac{\\tilde{r}^3}{r}\\sqrt{1 - \\frac{2m(r)}{\\tilde{r}}} &&\\text{with}\\qquad m(R)=M.
\\end{align}
```

`u` contains the initial values, and both `u` and `du` store values in the same order. 
`p` is a vector containing parameters needed during the integration.  The parameters passed in through `p` are:
- p[1] : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
- p[2] : The value of ``d\\tilde{r}/dr`` at ``r=0``, where ``\\tilde{r}`` is the Schwarzschild radial coordinate.

This version of the isotropic coordinate derivatives is intended for integrations starting from ``r=0``, and requires a
value for ``d\\tilde{r}/dr`` at ``r=0``.
"""
function TOVisotropic!(du,u,p,r)
    eos = p[1]
    
    r_tilde = u[1]  # areal radius; r is the isotropic radius
    β = u[2]  # metric function
    P = u[3] # pressure
    M = max(0, u[4])  # mass

    if r == 0.0 || P <= 0.0  
        du[1] = r == 0.0 ? p[2] : 0.0 
        du[2] = 0.0
        du[3] = 0.0
        du[4] = 0.0
    else
        fourpi = 4*pi
        a = sqrt(max(0,1 - (2*M / r_tilde))) # Avoid sqrt(-1) error
        ed = EnergyDensityfromPressure(eos,P) 
        du[1] = (r_tilde / r) * a
        du[2] = (M + fourpi * r_tilde^3 * P) / (r^2 * du[1])
        du[3] = -(ed + P) * du[2]
        du[4] = fourpi*ed * r_tilde^2*du[1]
    end
end

"""
    TOVisotropic(EOS, Ec, Rsurf, drtdr; 
                 SurfacePressure=0.0, AbsTol=1.0e-10, RelTol=1.0e-10, InitialStep=1.0e-4)

Integrate the TOV equations in isotropic radial coordinates using `TOVisotropic!`, starting from ``r=0`` out to `Rsurf`.
Integration stops when ``r`` reaches the surface radius `Rsurf`, or when the pressure falls below `SurfacePressure`,
whichever occurs first.
- `EOS`   : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
- `Ec`    : The dimensionless total energy density at ``r=0``.
- `Rsurf` : The dimensionless surface radius in isotropic coordinates.
- `drtdr` : The value of ``d\\tilde{r}/dr`` at ``r=0``, where ``\\tilde{r}`` is the Schwarzschild radial coordinate.

The solution is returned via the `TOVsolution` object.

##### Optional keyword arguements:
- `SurfacePressure` : The desired dimensionless pressure at the surface.
- `AbsTol`          : The absolute tolerance used during the integration.
- `RelTol`          : The relative tolerance used during the integration.
- `InitialStep`     : The initial step size for the integration relative to the surface radius.
"""
function TOVisotropic(eos::EOS,Ec::Number,Rsurf::Number, drdr::Number; SurfacePressure::Number=0.0,
                      AbsTol::Number=1.0e-10, RelTol::Number=1.0e-10, InitialStep::Number=1.0e-4)

    # ------------------------------------------------------------------------------------------------------
    # Keep this code in case the ContinuousCallback fails due to reaching Rsurf before finding the surface
    # pressure.
    #  -> Keyword argument needed my alternate surface finder:  SmallStep::Number=1.0e-12
    # # negativeP is used by isoutofdomain to prevent the pressure from becoming too small
    # # smallstep is a Callback function which treminates the integration if the step size becomes too small
    # negativeP(u,p,r) = u[3]<SurfacePressure ? true : false
    # smallstep(u,r,integrator) = (abs((integrator.t - integrator.tprev) /Rsurf) < SmallStep) ? true : false
    # ------------------------------------------------------------------------------------------------------
    
    findzeroP(u,r,integrator) = u[3]-SurfacePressure    # Callback to find the surface with pressure SurfacePressure
    
    u0 = [0, 0, PressurefromEnergyDensity(eos, Ec), 0]  # Initial conditions at center: R_tilde, β, P, M
    r_span = (0,Rsurf)  # Integration range

    prob = OrdinaryDiffEq.ODEProblem(TOVisotropic!, u0, r_span, (eos, drdr)) 

    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Vern9(), dt=InitialStep*Rsurf,
							   abstol=AbsTol, reltol=RelTol,
							   callback=OrdinaryDiffEq.ContinuousCallback(findzeroP,OrdinaryDiffEq.terminate!))
    # ------------------------------------------------------------------------------------------------------
    # Keep this code in case the ContinuousCallback fails due to reaching Rsurf before finding the surface
    # pressure.
    #           #isoutofdomain=negativeP,
    #           #callback=DiscreteCallback(smallstep,terminate!))
    # ------------------------------------------------------------------------------------------------------
    return TOVsolution(sol[4,end],Ec,u0[3],sol[3,end],sol.t[end],sol)
end

"""
    TOVshootout(EOS, Ec, Rsurf, drtdr; AbsTol=1.0e-10, RelTol=1.0e-10, InitialStep=1.0e-4)

Integrate the TOV equations in isotropic radial coordinates using `TOVisotropic!`, starting from ``r=0`` out to `Rsurf`.
- `EOS`   : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
- `Ec`    : The dimensionless total energy density at ``r=0``.
- `Rsurf` : The dimensionless surface radius in isotropic coordinates.
- `drtdr` : The value of ``d\\tilde{r}/dr`` at ``r=0``, where ``\\tilde{r}`` is the Schwarzschild radial coordinate.

The return value is the tuple `(pressure, sol)`, where
- `pressure` is the pressure at `Rsurf`.
- `sol`  is the full solution object returned by the `OrdinaryDiffEq` function `solve`

##### Optional keyword arguements:
- `AbsTol` : The absolute tolerance used during the integration.
- `RelTol` : The relative tolerance used during the integration.
- `InitialStep` : The initial step size for the integration relative to the surface radius.
"""
function TOVshootout(EOS::EOS,Ec::Number,Rsurf::Number, drtdr::Number;
                        AbsTol::Number=1.0e-10, RelTol::Number=1.0e-10,
                        InitialStep::Number=1.0e-4)
    
    u0 = [0, 0, PressurefromEnergyDensity(EOS, Ec), 0]  # Initial conditions at center: R_tilde, B_0, P, M
    r_span = (0,Rsurf)  # Integration range

    prob = OrdinaryDiffEq.ODEProblem(TOVisotropic!, u0, r_span, (EOS, drtdr)) 

    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Vern9(), dt=InitialStep*Rsurf, abstol=AbsTol,reltol=RelTol)
    return (sol[3,end],sol)
end

"""
    TOVstart(drtdr,p)

Wrapper routine to call TOVshootout as a function of `drtdr` and returning the difference between the integrated surface
pressure and the desired surface pressure.  Intended for use with NonlinearSolve::IntervalNonlinearProblem to find the correct
value of `drtdr` for outward integration of the TOV equations in isotropic radial coordinates.

The parameters passed in through `p` are:
- p[1] : The Equation of State struct.  Either a `PolytropicEOS` or a `RealisticEOS`.
- p[2] : The dimensionless total energy density at ``r=0``.
- p[3] : The dimensionless surface radius in isotropic coordinates.
- p[4] : The desired dimensionless pressure at the surface.
- p[5] : The absolute tolerance used during the integration.
- p[6] : The relative tolerance used during the integration.
- p[7] : The initial step size for the integration relative to the surface radius.
"""
function TOVstart(drtdr,p)
    sol=TOVshootout(p[1],p[2],p[3],drtdr,AbsTol=p[5], RelTol=p[6], InitialStep=p[7])
    return(p[4]-sol[1])
end

"""
    TOVModel(EOS, Ec, MaxR; 
             {SurfacePressure,QSurface}={1.01e+10,1.0e-10},
             AbsTol=1.0e-10, RelTol=1.0e-10, InitialStep=1.0e-4)

Integrate the TOV equations for a star with a specified central energy density.
- `EOS`  : The Equation of State object.  Either a `PolytropicEOS` or a `RealisticEOS`.
- `Ec`   : The dimensionless total energy density at ``r=0``.
- `MaxR` : A dimensionless radial integration limit that is larger than the Schwarzschild radius of the star.

The solution, computed in isotropic radial coordinates, is returned via the `TOVsolution` object.  All quantities are
dimensionless in Geometrized units where ``G=c=1``, but the scaling depends on the type of EOS used.

TOVModel first integrates the TOV equations in Schwarzschild radial coordinates starting at ``\\tilde{r}=0``, and finding the
radius of the star as defined by the condition that the dimensionless pressure ``\\bar{P}`` equal the value specified by either
`SurfacePressure` or `QSurface`.  For models based on a `RealisticEOS`, the dimensionless surface pressure is computed as 
``\\bar{P} = {\\tt SurfacePressure}/(\\epsilon_0 c^2)`` where ``\\epsilon_0=10^{15}\\text{g}/\\text{cm}^3`` and ``c`` is
the speed of light.  For models based on a `PolytropicEOS`, the dimensionless surface pressure is computed as 
``\\bar{P} = {\\tt QSurface}^{n+1}`` where ``n`` is the polytropic index.

Given a TOV model in Schwarzschild radial coordinates, the solution is converted to isotropic radial coordinates.  The solution
is first integrated in from the surface of the star toward ``r=0``.  However, integration towards the coordinate singularity at 
``r=0`` is unstable.  Integration outward from ``r=0`` is stable, but requires a value for ``d\\tilde{r}/dr`` at ``r=0``.
The inward integration is used to obtain an initial guess for this derivative, and then a shooting method is used to determine
the correct value of ``d\\tilde{r}/dr`` at ``r=0`` which will place the surface of the star at the same location for the 
given central energy density.  A **Warning** is printed if the relative error between the total masses obtained from the 
Schwarzschild and final isotropic integrations is larger than ``10^{-5}``.  The resulting solution in isotropic radial 
coordinates is returned.

##### Optional keyword arguements:
- `SurfacePressure` : The desired pressure at the surface.  This value is given in units of ``\\text{dyne}/\\text{cm}^2``,
                        and is used when EOS is of type RealisticEOS.
- `QSurface`        : The desired value of ``q=P/(\\rho_0c^2)`` at the surface.  This value is used when the
                        EOS is of type PolytropicEOS.
- `AbsTol`          : The absolute tolerance used during the integration.
- `RelTol`          : The relative tolerance used during the integration.
- `InitialStep`     : The initial step size for the integration relative to the dimensionless surface radius.

### Dimensionful results
Dimensionful results can be obtained from dimensionless code values.  Returning to dimensionful values for ``G`` and ``c``
in some desired units, we define ``\\kappa^{1/2}`` to set a fundamental length scale where ``\\kappa\\equiv\\frac{c^2}{G\\epsilon_0}``.
* For models constricted from a `RealisticEOS`, ``\\epsilon_0=10^{15}\\text{g}/\\text{cm}^3`` in cgs units.
* For models constructed from a `PolytropicEOS`, ``\\epsilon_0=\\left(\\frac{c^2}{K}\\right)^n``, where
    ``K`` is the polytropic constant and ``n`` the polytropic index relating pressure ``P`` to the rest 
    mass density ``\\rho_0`` via ``P= K\\rho_0^{1+1/n}``.  Note that ``K^{n/2}`` has units of length in gravitational units, 
    so its units depend on ``n``.  For example, a cold, non-relativistic, degenerate, ideal Fermi gas of neutrons has
    ``n=3/2`` and ``K=5.3802\\times10^9\\text{cm}^4/(\\text{g}^{2/3}\\text{s}^2)`` in cgs units.
The dimensionless (barred) quantities returned via the `TOVsolution` object are then related to dimensionful quantities via
- Radius   : ``r = \\kappa^{1/2}\\bar{r}``
- Mass     : ``M = \\kappa^{1/2}\\frac{c^2}{G}\\bar{M}``
- Pressure : ``P = \\kappa^{-1}\\frac{c^4}{G}\\bar{P}``
- Total Energy Density : ``\\epsilon = \\kappa^{-1}\\frac{c^4}{G}\\bar\\epsilon``
- Rest Mass Density : ``\\rho_0 = \\kappa^{-1}\\frac{c^2}{G}\\bar\\rho_0``
"""
function TOVModel(eos::PolytropicEOS,Ec::Number,MaxR::Number; QSurface::Number=1.0e-10,
                  AbsTol::Number=1.0e-10, RelTol::Number=1.0e-10, InitialStep::Number=1.0e-4)
    Polyn = eos.index;
    spress=QSurface^(Polyn+1);
    sol1 = TOVSchwarzschild(eos,Ec,MaxR; SurfacePressure=QSurface^(Polyn+1),
                            AbsTol=AbsTol, RelTol=RelTol, InitialStep=InitialStep)
    R = (sol1.SurfaceRadius/4)*(1+ sqrt(1 - (2*sol1.TotalMass / sol1.SurfaceRadius)))^2 # isotropic surface radius
    sol2 = integrateTOVin(eos,sol1; 
                          AbsTol=AbsTol, RelTol=RelTol, InitialStep=InitialStep);
    newdrdr = OrdinaryDiffEq.solve(
			  OrdinaryDiffEq.IntervalNonlinearProblem(TOVstart,
													  (0.9*sol2[1],1.1*sol2[1]),
													  (eos,Ec,R,spress,AbsTol,RelTol,InitialStep)))
    sol3 = TOVisotropic(eos,Ec,R,newdrdr[1];SurfacePressure=spress,
                        AbsTol=AbsTol, RelTol=RelTol, InitialStep=InitialStep)
    if (sol1.TotalMass-sol3.TotalMass)/sol1.TotalMass > 1.0e-5
        @warn "TOVModel produced relative error in the mass greater than 1.0e-5!" CentralDensity=Ec TOVmass=sol1.TotalMass Mass=sol3.TotalMass TOVpressure=sol1.SurfacePressure Pressure=sol3.SurfacePressure
    end
    return(sol3)
end
function TOVModel(eos::RealisticEOS,Ec::Number,MaxR::Number; SurfacePressure::Number=1.01e+10,
                  AbsTol::Number=1.0e-10, RelTol::Number=1.0e-10, InitialStep::Number=1.0e-4)
    spress=SurfacePressure/(Epsilon0*Cspeed^2)
    sol1 = TOVSchwarzschild(eos,Ec,MaxR; SurfacePressure=spress,
                           AbsTol=AbsTol, RelTol=RelTol, InitialStep=InitialStep)
    R = (sol1.SurfaceRadius/4)*(1+ sqrt(1 - (2*sol1.TotalMass / sol1.SurfaceRadius)))^2 # isotropic surface radius
    sol2 = integrateTOVin(eos,sol1; 
                          AbsTol=AbsTol, RelTol=RelTol, InitialStep=InitialStep);
    newdrdr = OrdinaryDiffEq.solve(
    		  OrdinaryDiffEq.IntervalNonlinearProblem(TOVstart,
    												  (0.9*sol2[1],1.1*sol2[1]),
    												  (eos,Ec,R,spress,AbsTol,RelTol,InitialStep)))
    sol3 = TOVisotropic(eos,Ec,R,newdrdr[1];SurfacePressure=spress,
                        AbsTol=AbsTol, RelTol=RelTol, InitialStep=InitialStep)
    if (sol1.TotalMass-sol3.TotalMass)/sol1.TotalMass > 1.0e-5
        @warn "TOVModel produced relative error in the mass greater than 1.0e-5!" CentralDensity=Ec TOVmass=sol1.TotalMass Mass=sol3.TotalMass TOVpressure=sol1.SurfacePressure Pressure=sol3.SurfacePressure
    end
    return(sol3)
end