
const global Cspeed = 2.99792458e+10  # cm/s
const global GNewton = 6.6743e-8 # cm^3/(g s)
const global Neutronmass = 1.6749275e-24 # g
const global Epsilon0 = 1.e+15 # g/cm^3

global huntindex::Int64 = 0 # Save the latest index into the current EOS
                            # This may need to be reconsidered for multi-threading

import 	HDF5

# EOS types
export 	EOS, PolytropicEOS, RealisticEOS, Cspeed, Epsilon0

# EOS transformations based on EOS
export 	EnergyDensityfromPressure,
		BaryonMassDensityfromPressure, 
		EnthalpyfromPressure,
		PressurefromEnergyDensity,
		BaryonMassDensityfromEnergyDensity,
		EnthalpyfromEnergyDensity,
		PressurefromEnthalpy,
		EnergyDensityfromEnthalpy,
		BaryonMassDensityfromEnthalpy

export mihunt

abstract type EOS end

"""
    PolytropicEOS(index)

Create a Polytropic EOS object with polytropic index `index`.

Can be used with the methods:
- `EnergyDensityfromPressure`
- `BaryonMassDensityfromPressure`
- `EnthalpyfromPressure`
- `PressurefromEnergyDensity`
- `BaryonMassDensityfromEnergyDensity`
- `EnthalpyfromEnergyDensity`
- `PressurefromEnthalpy`
- `EnergyDensityfromEnthalpy`
- `BaryonMassDensityfromEnthalpy`
to convert between dimensionless values of energy density ``\\bar\\epsilon``, 
baryon-mass density ``\\bar\\rho_0``, pressure ``\\bar{P}``, and enthalpy ``h``.

##### Available fields:
- `index` : The polytropic index ``n``.

##### Background
A polytropic EOS relates pressure ``P`` to rest-mass density (or baryon-mass density) ``\\rho_0`` or 
internal energy density ``\\rho_i``.  In geometrized units where ``G=c=1``, we have:
```math
\\begin{align} 
	P &=K\\rho_0^{1+\\frac1n}, \\\\
	&= (\\Gamma-1)\\rho_i,
\\end{align}
```
where the polytropic index ``n`` and adiabatic index ``\\Gamma`` are related by
```math
\\begin{align}
	\\Gamma = 1 + \\frac1n.
\\end{align}
```
Here, ``K`` is the polytropic constant, and the relationship between pressure and internal energy density 
assumes that all changes are adiabatic.

The total energy density is ``\\epsilon=\\rho_0+\\rho_i``.  The enthalpy ``h(P)`` is defined by:
```math
\\begin{align}
 h(P) \\equiv \\int^P_0{\\frac{{\\rm d}\\acute{P}}{\\epsilon(\\acute{P})+\\acute{P}}}.
\\end{align}
```
##### Dimensionless variables:
For polytropic EOS, and assuming geometrized units, ``K^{n/2}`` has units of length and is used to set
the fundamental length scale of the system.  The following dimensionless (barred) quantities are used:
- ``\\bar\\rho_0 \\equiv K^n\\rho_0 = q^n``
- ``\\bar\\rho_i \\equiv K^n\\rho_i = nq^{n+1}``
- ``\\bar{P} \\equiv K^nP = q^{n+1}``
- ``h`` = ``\\ln\\left(1+(n+1)q\\right)``
where ``q\\equiv\\frac{P}{\\rho_0} = \\frac{\\bar{P}}{\\bar\\rho_0}``.  Enthalpy ``h`` is 
a dimensionless quantity.

##### Dimensionfull quantities:
To obtain dimensionfull quantities in some desired units requires a value of ``K`` in those units.  The
dimensions of ``K`` depend on ``n``.  For example, a cold, non-relativistic, degenerate, ideal Fermi gas 
of neutrons has ``n=3/2`` and ``K=5.3802\\times10^9\\text{cm}^4/(\\text{g}^{2/3}\\text{s}^2)`` in cgs units.
Given a value for ``K``, we have:
```math
\\begin{align}
	\\rho_0 &= c^{2n}K^{-n}\\bar\\rho_0 \\qquad\\ \\ :\\ \\text{with units of mass/volume} \\\\
	\\rho_i &= c^{2n}K^{-n}\\bar\\rho_i \\qquad\\ \\ :\\ \\text{with units of mass/volume} \\\\
	\\epsilon &= \\left\\{\\begin{array}{ll}
		c^{2n}K^{-n}\\bar\\epsilon &\\,:\\ \\text{with units of mass/volume} \\\\
		c^{2(n+1)}K^{-n}\\bar\\epsilon &\\,:\\ \\text{with units of energy/volume} 
		\\end{array}\\right. \\\\
	P &= c^{2(n+1)}K^{-n}\\bar{P} \\quad\\ :\\ \\text{with units of force/area}
\\end{align}
```
"""
struct PolytropicEOS <: EOS
	index::Float64
end

"""
    RealisticEOS("EOSfile","EOSname")

Create an EOS object based on the `"EOSname"` stored in the HDF5 file `"EOSfile"`.

Can be used with the methods:
- `EnergyDensityfromPressure`
- `BaryonMassDensityfromPressure`
- `EnthalpyfromPressure`
- `PressurefromEnergyDensity`
- `BaryonMassDensityfromEnergyDensity`
- `EnthalpyfromEnergyDensity`
- `PressurefromEnthalpy`
- `EnergyDensityfromEnthalpy`
- `BaryonMassDensityfromEnthalpy`
to convert between dimensionless values of energy density ``\\bar\\epsilon``, 
baryon-mass density ``\\bar\\rho_0``, pressure ``\\bar{P}``, and enthalpy ``h``.

##### Available fields:
- `name` : The name of the EOS - `"EOSname"`.
- `baryonmass` : The value of the baryon mass used for this EOS.
- `phasetransition` : Does the EOS contain any phase transitions.
- `nphasechange` : Number of phase transitions within the EOS.
- `phasePs` : Vector of values of the pressure at each phase change.
- `phaseHs` : Vector of values of the enthalpy at each phase change.
- `data` : The raw, dimensionfull data stored in the HDF5 file `"EOSfile"`.
- `bmd` : Vector of natural logs of dimensionless baryon mass density.
- `tmd` : Vector of natural logs of dimensionless total mass density.
- `P` : Vector of natural logs of dimensionless pressure.
- `H` : Vector of natural logs of enthalpy.

##### Background
A realistic EOS consists of tables of values for the pressure ``P``, total energy density ``\\epsilon``,
and baryon-mass density (or rest-mass density) ``\\rho_0``.  Currently used tables store the base-10 log
of these quantities in cgs units.  A tabulated EOS typically has about 500 entries which are roughly 
evenly spaced in the base-10 log of the total energy density.  A table starts at about 
``7.86\\;\\text{g}\\,\\text{cm}^{-3}`` (beginning of the FMT EOS) and ends wherever the given "high-density" EOS 
ends (typically around ``5\\sim10\\times10^{16}\\text{g}\\,\\text{cm}^{-3}``).

When a `RealisticEOS` object is created, the values in a tabulated EOS are converted to dimensionless values 
in geometrized units where ``G=c=1``, and stored as the natural log of the dimensionless values.  Linear
interpolation of these natural logs is used to treat each stored quantity as a function of any other.
The natural log of the enthalpy ``h(P)`` is also compute, with ``h(P)`` defined by:
```math
\\begin{align}
 h(P) \\equiv \\int^P_0{\\frac{{\\rm d}\\acute{P}}{\\epsilon(\\acute{P})+\\acute{P}}},
\\end{align}
```
For tabulated EOS which contain phase transitions, interpolation through a phase transition must
be handled carefully.  If any quantity is interpolated as a function of either pressure or enthalpy,
then that quantity may have multiple values at the phase transition.  Interpolation always returns the 
smallest value.


##### Dimensionless variables:
For realistic EOS, we choose the fundamental length scale to be defined by ``\\kappa^{1/2}``, where
```math
\\begin{align}
	k\\equiv\\frac{c^2}{G\\epsilon_0}
\\end{align}
```
where ``\\epsilon_0=10^{15}\\text{g}\\,\\text{cm}^{-3}``. The following dimensionless (barred) 
quantities are used:
- ``\\bar\\rho_0 \\equiv \\kappa\\frac{G}{c^2}\\rho_0``
- ``\\bar\\epsilon \\equiv \\kappa\\frac{G}{c^2}\\epsilon``
- ``\\bar{P} \\equiv \\kappa\\frac{G}{c^4}P``
Enthalpy ``h`` is a dimensionless quantity.  *Note that tabulated values of ``\\epsilon`` are stored as
total mass density, and not as total energy density*.

##### Dimensionfull quantities:
To obtain dimensionfull quantities in some desired units requires a value of ``\\kappa`` in those units.
```math
\\begin{align}
	\\rho_0 &= \\kappa^{-1}\\frac{c^2}{G}\\bar\\rho_0 \\quad:\\ \\text{with units of mass/volume} \\\\
	\\epsilon &= \\left\\{\\begin{array}{ll}
		\\kappa^{-1}\\frac{c^2}{G}\\bar\\epsilon &\\,:\\ \\text{with units of mass/volume} \\\\
		\\kappa^{-1}\\frac{c^4}{G}\\bar\\epsilon &\\,:\\ \\text{with units of energy/volume} 
		\\end{array}\\right. \\\\
	P &= \\kappa^{-1}\\frac{c^4}{G}\\bar{P} \\quad:\\ \\text{with units of force/area}
\\end{align}
```
"""
struct RealisticEOS <: EOS
	name::String
	baryonmass::Float64 # Value of the baryon mass used for this EOS
	phasetransition::Bool # Does the EOS contain any phase transitions
	nphasechange::Int64 # Number of phase changes in the EOS
	phasePs # Storage for values of pressure at each phase transition
	phaseHs # Storage for values of enthalpy at each phase transition
	phaseIndices # Storage for the indices for each phase change
	data # The full raw data from the HDF5 file 
	bmd::Vector{Float64} # Vector of logs of dimensionless Baryon mass density
	tmd::Vector{Float64} # Vector of logs of dimensionless total mass density
	P::Vector{Float64} # Vector of logs of dimensionless pressure
	H::Vector{Float64} # Vector of logs of dimensionless enthalpy
	RealisticEOS(name::String,barymass::Float64,phase::Bool,nph,phPs,phHs,phInd,data,bmd,tmd,P,H) = new(name,barymass,phase,nph,phPs,phHs,phInd,data,bmd,tmd,P,H)
	RealisticEOS(eosfile::String,name::String) = readEOSdata(eosfile,name)
end

"""
    readEOSdata("EOSfile","EOSname")

Construct a RealisticEOS object based on the `"EOSname"` stored in the HDF5 file `"EOSfile"`.
"""
function readEOSdata(eosfile::String,name::String)
	# Read the HDF file, find and import the EOS data, pass it to the constructor
	try
		HDF5.h5open(eosfile,"r") do h5f 
		eosgroup=0 # needed to put eosgroup in the proper scope
			try
				eosgroup = h5f[name]
			catch
				# For future, if group named name is not found in HDF file,
				# look for group name 'EOS' and try to read the EOS data from there
				println("Failed to read $name from $eosfile")
				exit()
			end
			try
				storedname = HDF5.attrs(eosgroup)["EOS Name"]
				if name != storedname
					println("EOS names stored name are different")
					exit()
				end
				barymass = HDF5.attrs(eosgroup)["BaryonMass"]
				phase = HDF5.attrs(eosgroup)["PhaseTransition"]
				data = permutedims(HDF5.read_dataset(eosgroup,"Table"))
				bmd = log(barymass/Epsilon0) .+ data[:,1]*log(10)
				tmd = -log(Epsilon0) .+ data[:,2]*log(10)
				P = -log(Epsilon0*Cspeed^2) .+ data[:,3]*log(10)
				H = logenthalpy(data)
				#
				# Set up to handle interpolation with phase transitions
				nphasechange::Int64 = 0 # Number of phase changes in the EOS
				phasePs = Float64[] # values of pressure at each phase transition
				phaseHs = Float64[] # values of enthalpy at each phase transition
				phaseIndices = Array{Int64,2} # indices for each phase change
				if phase 
					phaseIndices = HDF5.attrs(eosgroup)["TransitionIndices"]
					# phaseIndices = HDF5.read_attribute(eosgroup,"TransitionIndices")
					nphasechange = size(phaseIndices)[2]
					for i=1:nphasechange
						push!(phasePs,P[phaseIndices[1,i]])
						push!(phaseHs,H[phaseIndices[1,i]])
					end
				end
				#
				return(RealisticEOS(storedname,barymass,phase,nphasechange,phasePs,phaseHs,phaseIndices,data,bmd,tmd,P,H))
			catch
				println("Failed to read Table from $eosfile")
				exit()
			end
		end
	catch
		println("Failed to open $eosfile")
		exit()
	end
end

"""
    EnergyDensityfromPressure(eos,P)

Return the dimensionless energy density associated with the dimensionless pressure `P`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.

For `RealisticEOS` containing phase transitions, if `P` is at a phase transition, then 
the smallest value of the energy density is returned.
"""
function EnergyDensityfromPressure(eos::EOS,P::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    BaryonMassDensityfromPressure(eos,P)

Return the dimensionless baryon mass density associated with the dimensionless pressure `P`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.

For `RealisticEOS` containing phase transitions, if `P` is at a phase transition, then 
the smallest value of the baryon mass density is returned.
"""
function BaryonMassDensityfromPressure(eos::EOS,P::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    EnthalpyfromPressure(eos,P)

Return the dimensionless enthalpy associated with the dimensionless pressure `P`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.

For `RealisticEOS` containing phase transitions, if `P` is at a phase transition, then 
the smallest value of the enthalpy is returned.
"""
function EnthalpyfromPressure(eos::EOS,P::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    PressurefromEnergyDensity(eos,E)

Return the dimensionless pressure associated with the dimensionless energy density `E`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.
"""
function PressurefromEnergyDensity(eos::EOS,E::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    BaryonMassDensityfromEnergyDensity(eos,E)

Return the dimensionless baryon mass density associated with the dimensionless energy density `E`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.
"""
function BaryonMassDensityfromEnergyDensity(eos::EOS,E::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    EnthalpyfromEnergyDensity(eos,E)

Return the dimensionless enthalpy associated with the dimensionless energy density `E`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.
"""
function EnthalpyfromEnergyDensity(eos::EOS,E::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    PressurefromEnthalpy(eos,H)

Return the dimensionless pressure associated with the dimensionless enthalpy `H`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.

For `RealisticEOS` containing phase transitions, if `H` is at a phase transition, then 
the smallest value of the pressure is returned.
"""
function PressurefromEnthalpy(eos::EOS,H::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    EnergyDensityfromEnthalpy(eos,H)

Return the dimensionless energy density associated with the dimensionless enthalpy `H`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.

For `RealisticEOS` containing phase transitions, if `H` is at a phase transition, then 
the smallest value of the energy density is returned.
"""
function EnergyDensityfromEnthalpy(eos::EOS,H::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end
"""
    BaryonMassDensityfromEnthalpy(eos,H)

Return the dimensionless baryon mass density associated with the dimensionless enthalpy `H`
as defined by the EOS `eos`.  `eos` can be either a `PolytropicEOS` or a `RealisticEOS`.

For `RealisticEOS` containing phase transitions, if `H` is at a phase transition, then 
the smallest value of the baryon mass density is returned.
"""
function BaryonMassDensityfromEnthalpy(eos::EOS,H::Number)::Float64
	println("Need to construct EOS functions for new composite EOS type")
	exit()
end

"""
    Etoq(E,n)

For polytropic EOS with polytropic index `n`, return the value of `q` associated with 
the dimensionless total energy density `E`.  `q` is the ratio of the pressure to the 
rest-mass density.
"""
function Etoq(E,n)
	# Use Newton's method to solve E=q^n(1+n*q) for q
	q=E^(1/n)
	err=1
	while err > 1.e-10
		dq = E - (1+n*q)*(q^n)
		dq /= n*(1+(n+1)*q)*(q^(n-1))
		q += dq
		err = abs(dq/q)
	end
	return q
end
# EOS functions for composite PolytropicEOS
# set up conversions based on
# dimensionless rest-mass density (baryon-mass density) = q^n
# dimensionless internal energy density = n*q^(n+1)
# dimensionless total energy density = (1+n*q)q^n
# dimensionless pressure = q^(n+1)
# dimensionless enthalpy = ln(1+(n+1)q)
function EnergyDensityfromPressure(eos::PolytropicEOS,P::Number)::Float64
	n=eos.index
	P==0 ? 0 : n*P + P^(n/(n+1))
end
function BaryonMassDensityfromPressure(eos::PolytropicEOS,P::Number)::Float64
	n=eos.index
	P==0 ? 0 : P^(n/(n+1))
end
function EnthalpyfromPressure(eos::PolytropicEOS,P::Number)::Float64
	n=eos.index
	P==0 ? 0 : log1p((n+1)*(P^(1/(n+1))))
end
function PressurefromEnergyDensity(eos::PolytropicEOS,E::Number)::Float64
	n=eos.index
	E==0 ? 0 : (q=Etoq(E,n); q^(n+1))
end
function BaryonMassDensityfromEnergyDensity(eos::PolytropicEOS,E::Number)::Float64
	n=eos.index
	E==0 ? 0 : (q=Etoq(E,n); q^n)
end
function EnthalpyfromEnergyDensity(eos::PolytropicEOS,E::Number)::Float64
	n=eos.index
	E==0 ? 0 : (q=Etoq(E,n); log1p((n+1)*q))
end
function PressurefromEnthalpy(eos::PolytropicEOS,H::Number)::Float64
	n=eos.index
	H==0 ? 0 : (q=expm1(H)/(n+1); q^(n+1))
end
function EnergyDensityfromEnthalpy(eos::PolytropicEOS,H::Number)::Float64
	n=eos.index
	H==0 ? 0 : (q=expm1(H)/(n+1); (1+n*q)*(q^n))
end
function BaryonMassDensityfromEnthalpy(eos::PolytropicEOS,H::Number)::Float64
	n=eos.index
	H==0 ? 0 : (q=expm1(H)/(n+1); q^n)
end

"""
    logenthalpy(data)

Return the natural log of the enthalpy associated with the raw tabulated EOS `data`.
An enthalpy value is returned for each entry in the EOS `data`, assuming the first EOS entry
corresponds to the lowest density in the table.

The integrand is written as ``\\frac{P{\\rm d}\\ln{P}}{\\epsilon + P}``, where ``\\epsilon`` is 
the total energy density and ``P`` is the pressure.  The initial value is approximated by 
``\\frac12\\frac{P_1}{\\epsilon_1+P_1}``, and the trapezoidal summed integrand 
is expressed as ``\\frac12\\ln\\left(\\frac{P_i}{P_{i-1}}\\right)\\left(\\frac{P_i}{\\epsilon_i+P_i}+\\frac{P_{i-1}}{\\epsilon_{i-1}+P_{i-1}}\\right)``.
"""
function logenthalpy(data)
	dlnp = log(10.0)*(data[2:end,3]-data[1:end-1,3])
	enth = 0.5./(1 .+ (Cspeed^2)*exp10.(data[:,2]-data[:,3]))
	dlnp .*= enth[1:end-1] + enth[2:end]
	enth[2:end] = dlnp
	return log.(accumulate(+,enth))
end

"""
    mihunt(xdata,x,guess)

Given an array of monotonically increasing data `xdata`, return the index `i` into `xdata` such
that `xdata[i] <= x <xdata[i+1]`.  If `0 < guess <= length(xdata)`, then `guess` is used as a 
starting location for the search.
"""
function mihunt(EOSx::Vector{Float64},x::Number,jguess::Int)::Int
	neos::Int = length(EOSx)
	if (neos < 2) println("Need at least 2 data points for interpolation"); exit() end
	jl::Int = jguess; ju::Int = jguess; jm::Int = 0; inc::Int = 1
	if ((jl < 1) | (jl > neos)) jl = 1; ju = neos # no initial guess
	elseif x > EOSx[jl] # find bracket from below
		while true
			ju = jl + inc
			if (ju > neos) ju = neos; break
			elseif (x < EOSx[ju]) break
			else jl = ju; inc += inc
			end
		end
	else # find bracket from above
		ju = jl
		while true
			jl = jl - inc
			if (jl < 1) jl = 1; break
			elseif (x > EOSx[jl]) break
			else ju = jl; inc += inc
			end
		end
	end
	while ju - jl > 1 # Use bisection to find index below or at value
		jm = (ju+jl) >> 1 # same as floor((ju+jl)/2)
		if (x >= EOSx[jm]) jl = jm
		else ju = jm
		end
	end
	return clamp(jl,1,neos-1) # ensure value is 1<=jl<=neos-1
end
# EOS functions for composite RealisticEOS
# set up conversions based linear interpolation
function EnergyDensityfromPressure(eos::RealisticEOS,Pv::Number)::Float64
	logP = eos.P; logE = eos.tmd
	logPv = log(Pv)
	j = huntindex
	if (eos.nphasechange>0) & in(logPv,Set(eos.phasePs))
		# guarantee that the beginning index for a phase transition is used
		for i=1:eos.nphasechange
			if eos.phasePs[i]==logPv
				global huntindex = j = eos.phaseIndices[1,i]
				return exp(logE[j])
			end
		end
	else
		global huntindex = j = mihunt(logP,logPv,j)
	end
	return exp(logE[j] + (logE[j+1]-logE[j])*(logPv - logP[j])/(logP[j+1]-logP[j]))
end
function BaryonMassDensityfromPressure(eos::RealisticEOS,Pv::Number)::Float64
	logP = eos.P; logB = eos.bmd
	logPv = log(Pv)
	j = huntindex
	if (eos.nphasechange>0) & in(logPv,Set(eos.phasePs))
		# guarantee that the beginning index for a phase transition is used
		for i=1:eos.nphasechange
			if eos.phasePs[i]==logPv
				global huntindex = j = eos.phaseIndices[1,i]
				return exp(logB[j])
			end
		end
	else
		global huntindex = j = mihunt(logP,logPv,j)
	end
	return exp(logB[j] + (logB[j+1]-logB[j])*(logPv - logP[j])/(logP[j+1]-logP[j]))
end
function EnthalpyfromPressure(eos::RealisticEOS,Pv::Number)::Float64
	logP = eos.P; logH = eos.H
	logPv = log(Pv)
	j = huntindex
	if (eos.nphasechange>0) & in(logPv,Set(eos.phasePs))
		# guarantee that the beginning index for a phase transition is used
		for i=1:eos.nphasechange
			if eos.phasePs[i]==logPv
				global huntindex = j = eos.phaseIndices[1,i]
				return exp(logH[j])
			end
		end
	else
		global huntindex = j = mihunt(logP,logPv,j)
	end
	return exp(logH[j] + (logH[j+1]-logH[j])*(logPv - logP[j])/(logP[j+1]-logP[j]))
end

function PressurefromEnergyDensity(eos::RealisticEOS,Ev::Number)::Float64
	logE = eos.tmd; logP = eos.P
	logEv = log(Ev)
	j = huntindex
	global huntindex = j = mihunt(logE,logEv,j)
	return exp(logP[j] + (logP[j+1]-logP[j])*(logEv - logE[j])/(logE[j+1]-logE[j]))
end
function BaryonMassDensityfromEnergyDensity(eos::RealisticEOS,Ev::Number)::Float64
	logE = eos.tmd; logB = eos.bmd
	logEv = log(Ev)
	j = huntindex
	global huntindex = j = mihunt(logE,logEv,j)
	return exp(logB[j] + (logB[j+1]-logB[j])*(logEv - logE[j])/(logE[j+1]-logE[j]))
end
function EnthalpyfromEnergyDensity(eos::RealisticEOS,Ev::Number)::Float64
	logE = eos.tmd; logH = eos.H
	logEv = log(Ev)
	j = huntindex
	global huntindex = j = mihunt(logE,logEv,j)
	return exp(logH[j] + (logH[j+1]-logH[j])*(logEv - logE[j])/(logE[j+1]-logE[j]))
end

function PressurefromEnthalpy(eos::RealisticEOS,Hv::Number)::Float64
	logH = eos.H; logP = eos.P
	logHv = log(Hv)
	j = huntindex
	if (eos.nphasechange>0) & in(logHv,Set(eos.phaseHs))
		# guarantee that the beginning index for a phase transition is used
		for i=1:eos.nphasechange
			if eos.phaseHs[i]==logHv
				global huntindex = j = eos.phaseIndices[1,i]
				return exp(logP[j])
			end
		end
	else
		global huntindex = j = mihunt(logH,logHv,j)
	end
	return exp(logP[j] + (logP[j+1]-logP[j])*(logHv - logH[j])/(logH[j+1]-logH[j]))
end
function EnergyDensityfromEnthalpy(eos::RealisticEOS,Hv::Number)::Float64
	logH = eos.H; logE = eos.tmd
	logHv = log(Hv)
	j = huntindex
	if (eos.nphasechange>0) & in(logHv,Set(eos.phaseHs))
		# guarantee that the beginning index for a phase transition is used
		for i=1:eos.nphasechange
			if eos.phaseHs[i]==logHv
				global huntindex = j = eos.phaseIndices[1,i]
				return exp(logE[j])
			end
		end
	else
		global huntindex = j = mihunt(logH,logHv,j)
	end
	return exp(logE[j] + (logE[j+1]-logE[j])*(logHv - logH[j])/(logH[j+1]-logH[j]))
end
function BaryonMassDensityfromEnthalpy(eos::RealisticEOS,Hv::Number)::Float64
	logH = eos.H; logB = eos.bmd
	logHv = log(Hv)
	j = huntindex
	if (eos.nphasechange>0) & in(logHv,Set(eos.phaseHs))
		# guarantee that the beginning index for a phase transition is used
		for i=1:eos.nphasechange
			if eos.phaseHs[i]==logHv
				global huntindex = j = eos.phaseIndices[1,i]
				return exp(logB[j])
			end
		end
	else
		global huntindex = j = mihunt(logH,logHv,j)
	end
	return exp(logB[j] + (logB[j+1]-logB[j])*(logHv - logH[j])/(logH[j+1]-logH[j]))
end
