using DSAModels
using AnalyticMHDTestSolutions

# ciza shock properties
GU = GadgetPhysical()
M  = 4.6
T2 = 2.588e8 # K 
ρ2 = 9.4e-27 # g/cm^3
γ  = 5/3

xmin = 50.0
xmax = 100.0

ρr = ρ2 / ( (γ + 1 ) * M^2 / ( (γ - 1 ) * M^2 + 2) ) / GU.rho_cgs
Ur = T2 / ( (2γ*M^2 - (γ - 1))*( (γ-1)*M^2 + 2) / ( (γ+1)^2 * M^2) ) / GU.T_K
Pr = (γ-1) * ρr * Ur

ρl = 8ρr
Pl = AnalyticMHDTestSolutions.solvePlfromMach(ρl, ρr, Pr, M, γ)

h = head_to_obj(fi)
time_end = 1.5

# set up rieman problem
par = RiemannParameters(rhol=ρl, Pl=Pl, rhor=ρr, Pr=Pr, t=time_end, dsa_model = 0, xs_first_guess=3.8)

# find compression factor for shock
xs = find_xs_first_guess(par.Ul, 4.6, dsa_model = 0)

# solve riemann problem
sol = AnalyticMHDTestSolutions.solve([0.0], par)

# print modified Mach number due to CR injection
println(sol.Mach)