#-----------kitaev paper chiral soliton------------
using ITensors, ITensorMPS
using Plots
using LaTeXStrings
#using CSV
#using DataFrames
gr()
 
let
    N = 384 # lenght of chains
    sites = siteinds("S=1/2",N)
    os = OpSum()
    bondenergy = OpSum()
    J = 1 # J <-> coupling
    K = 1
    h_xy = 0.55 # magnetic field strength in xy direction(s)
    h_z = 0 # im pretty sure this is same as h_xy, but since theta_z = 0 the magnetic field term for z direction is zero anyways
 
    h_1c = 0.511 # critical point 1 
    h_2c = 0.726 # critical point 2
 
    halfway::Int = ceil(N/2)
    phi_xy = pi/4
    theta_z = 0

 
    a = cos(phi_xy)*cos(theta_z)
    b = sin(phi_xy)*cos(theta_z)
    c = sin(theta_z)
 
#kitaev spin chain hamiltonian (alternating interactions: odd i have Sx at i and i+1, even i have Sy at i and i+1)
for i in 1:N-1
    if isodd(i)
        os += J, "Sx", i, "Sx", i+1
    else
        os += K, "Sy", i, "Sy", i+1
    end
end
 
# kitaev magnetic field term, first half
for i in 1:N
        os -= a*h_xy, "Sx", i
        os -= b*h_xy, "Sy", i
        os -= c*h_xy, "Sz", i            
end
H = MPO(os,sites)
    #psi0 = random_mps(sites;linkdims=60)
    states  = [i <= halfway ? "Up" : "Dn" for i in 1:N]
    # initial MPS with NO randomness in bonds
    psi0    = random_mps(sites, states; linkdims = 200)
    #orthogonalize psi0 at halfway point
    #orthogonalize!(psi0, halfway)
 
    # DMRG parameters
    nsweeps = 1000 # increase number of sweeps until it converges. For lower h term, LOTS more sweeps are needed (500+)
    maxdim = [800, 800, 1000, 1000, 1200, 1200, 1200,2000,2000,3000] # maxdim for convergence. Kitaev paper uses "more than 1000" for the bond dim.
    cutoff = [1E-10] # adjust this to tune accuracy of convergence. 1E-12 is "almost exact." kitaev paper uses less than 1E-10.
    noise = [1E-6]
    print("\n------Ground State------\n")
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    print("\nGround State Energy =", energy,"\n")
 
e0 = energy/N # average energy per site
 
# computing bond energy e_i
ei = zeros(1, N-1)
# new opsum for local hamiltonian terms
for i in 1:N-1
    if isodd(i)
        bondenergy += J, "Sx", i, "Sx", i+1
    else
        bondenergy += K, "Sy", i, "Sy", i+1
    end
       bondenergy -= a*h_xy, "Sx", i
       bondenergy -= b*h_xy, "Sy", i
end
 
# magnetic field for endpoint of chain    
bondenergy -= a*h_xy, "Sx", N
bondenergy -= b*h_xy, "Sy", N
# initalize endpoints
#ei[1] = real(inner(psi' , MPO(0.5*bondenergy[1] + bondenergy[2] + bondenergy[3], sites) , psi))
#ei[N-1] = real(inner(psi' , MPO(0.5*bondenergy[length(bondenergy)-2] + bondenergy[length(bondenergy)-1] + bondenergy[length(bondenergy)], sites) , psi))
 
 
    
# construct e_i's (with expectation value)
for i in 4:3:length(bondenergy)-3
    k = Int(ceil(i/3)) 
    ei[k] = real(inner(psi' , MPO(0.5*bondenergy[i] + bondenergy[i+1] + bondenergy[i+2] + 0.5*bondenergy[i-3], sites) , psi))
end
 
# soliton mass
#soli_mass = 0 
#for i in 1:N-1
#    soli_mass += ei[i]
#end  
 
   
 
# local hamiltonian terms
#for i in 1:3:N
#    hi = MPO(bondenergy[i]+bondenergy[i+1]+bondenergy[i+2], sites)
#    print(hi)
#    ei[i] = real(inner(psi', hi, psi))
#end
 
 
function spin_expectation(psi, op::AbstractString, site::Int)
    return real(expect(psi, op; sites=site))
end
 
Sx_vals = zeros(N)
Sy_vals = zeros(N)
    for i in 1:N
        Sx_vals[i] = spin_expectation(psi, "Sx", i)
        Sy_vals[i] = spin_expectation(psi, "Sy", i)
    end
 
print("\n------Sx VALUES AND Sy VALUES-----\n")
print(Sx_vals)
print("\n")
print(Sy_vals)
print("\n")
odd_sites = 1:2:N
even_sites = 2:2:N
 
    
 
    p1 = plot(
        odd_sites, Sx_vals[odd_sites],
        label = L"\langle S^x_i \rangle", 
        linewidth = 2, 
        color = :green,
        xlabel = "Site i",
        ylabel = "Spin Expectation Value",
        title = "Spin Expectation vs. Site (Odd Sites)"
    )
    plot!(
        odd_sites, Sy_vals[odd_sites],
        label = L"\langle S^y_i \rangle",
        linewidth = 2,
        color = :red,
        linestyle = :dash
    )
display(p1)
    p2 = plot(
        even_sites, Sx_vals[even_sites],
        label = L"\langle S^x_i \rangle", 
        linewidth = 2, 
        color = :green,
        xlabel = "Site i",
        ylabel = "Spin Expectation Value",
        title = "even sites"
    )
    plot!(
        even_sites, Sy_vals[even_sites],
        label = L"\langle S^y_i \rangle",
        linewidth = 2,
        color = :red,
        linestyle = :dash
    )
    display(p2)
 
    
ei_modified = zeros(1,N-1)
# create modified ei function
for i in 1:N-1
    ei_modified[i] = ei[i]-e0
 
end
    p3 = plot(
        1:N-1, ei_modified[1:N-1],
        label = L"\langle e_i \rangle-e_0^{XY}", 
        linewidth = 2, 
        color = :green,
        xlabel = "Site i",
        ylabel = L"\langle e_i \rangle-e_0^{XY}",
        title = "soliton?"
    )
 
    display(p3)
 
end
