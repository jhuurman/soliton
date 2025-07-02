# testing hamiltonian construction
using ITensors, ITensorMPS
using Plots
using LaTeXStrings

gr()


function second_derivative_nonuniform(y::Vector{Float64}, x::Vector{Float64})
    n = length(y)
    d2y = zeros(n - 2) 

    for i in 1:(n-2)
        x_prev = x[i]
        x_curr = x[i+1]
        x_next = x[i+2]

        y_prev = y[i]
        y_curr = y[i+1]
        y_next = y[i+2]


        h_prev = x_curr - x_prev
        h_next = x_next - x_curr


        term_next = (y_next - y_curr) / h_next
        term_prev = (y_curr - y_prev) / h_prev
        
        d2y[i] = 2 * (term_next - term_prev) / (h_prev + h_next)
    end
    return d2y 
end


let
    N = 99 # length of chain
    sites = siteinds("S=1/2",N)
    J = 1.0 # J <-> coupling
    K = 1.0
    theta_z = 0.0 

    h_values = [0.6,0.65,0.7] 
    theta_values_denomenator_range = 3.5:0.025:4.75
    theta_values_denomenator = reverse(collect(theta_values_denomenator_range))


    #-----DATA STORAGE FOR LATER-----
    calculated_e0_data = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    susceptibility_results = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    
    # the big loop
    for h in h_values
        print("Calculating for h = ", h, "\n")
        local_d_values = Float64[]       # denominators d for this h 
        local_e0_values = Float64[]      # e0 (energy/N) values for this h

        for d_val in theta_values_denomenator 
            phi_xy = pi / d_val 
            print("  d_val = ", d_val, ", phi_xy = ", round(phi_xy, digits=4), " rad (", round(rad2deg(phi_xy), digits=2)," deg)\n")
            
            halfway::Int = ceil(N/2)
            
            os_total_H = OpSum() 
            
            a = cos(phi_xy) * cos(theta_z)
            b = sin(phi_xy) * cos(theta_z)
            c = sin(theta_z)
    
            # kitaev spin chain hamiltonian
            for i in 1:N-1
                if isodd(i)
                    os_total_H += J, "Sx", i, "Sx", i+1
                else
                    os_total_H += K, "Sy", i, "Sy", i+1
                end
            end
            
            # magnetic field term for total hamiltonian
            for i in 1:N
                os_total_H -= a*h, "Sx", i
                os_total_H -= b*h, "Sy", i
                os_total_H -= c*h, "Sz", i
            end
            
            H_mpo = MPO(os_total_H, sites)
            
            states = [i <= halfway ? "Up" : "Dn" for i in 1:N]
            psi0 = random_mps(sites, states; linkdims=30)
            orthogonalize!(psi0, halfway)
            
            # DMRG parameters
            nsweeps = 1000 
            maxdim = [600, 600, 800, 800, 1200, 1200, 1200, 2000, 2000, 3000]
            cutoff = [1E-10] 
            noise = [1E-6, 1E-7, 1E-8, 1E-9, 0.0] 
            
            energy, psi = dmrg(H_mpo, psi0; nsweeps, maxdim, cutoff, noise, outputlevel=0)
        
            e0 = energy / N 
        
            push!(local_d_values, d_val) 
            push!(local_e0_values, e0)
            
        end # end of d_val (phi_xy) loop
        print("Finished calculations for h = ", h, "\n")
            
        calculated_e0_data[h] = (local_d_values, local_e0_values)
            
    end # end of h loop
    
    # -------- END OF DMRG LOOPS --------
    
    print("\nDMRG FINISHED. CALCULATING ENERGETIC SUSCEPTIBILITY (from e0)...\n")
    
    h_plot_values_sorted = sort(collect(keys(calculated_e0_data))) 
    
    for h_val in h_plot_values_sorted
        denominators_decreasing, e0_for_h = calculated_e0_data[h_val] 
        
        thetas_phi_xy_increasing = pi ./ denominators_decreasing 


        d2_e0_dtheta2 = second_derivative_nonuniform(e0_for_h, thetas_phi_xy_increasing)

        susceptibility_chi_e = -d2_e0_dtheta2

        thetas_for_derivative = thetas_phi_xy_increasing[2:end-1] 
        susceptibility_results[h_val] = (thetas_for_derivative, susceptibility_chi_e)
        print("calculated energetic susceptibility for h = ", round(h_val, digits=3), "\n")
    end
        
    print("\nPLOTTING ENERGETIC SUSCEPTIBILITY (derived from e0)...\n")
    
    colors_chi_e = cgrad(:viridis, length(h_plot_values_sorted), categorical=true)

    final_plot_chi_e = plot(
        xlabel = L"\phi_{xy}", 
        ylabel = L"\chi^e_{\phi_{xy}}", 
        title = "Energetic Susceptibility vs. Field Angle", 
        legend = :topright 
    )
    
    plot_idx = 1
    for h_val in h_plot_values_sorted
        if haskey(susceptibility_results, h_val)
            thetas_rad, susceptibility = susceptibility_results[h_val]
            thetas_deg = rad2deg.(thetas_rad) 
            plot!(
                final_plot_chi_e,
                thetas_deg,
                susceptibility,
                label = "h = $(round(h_val, digits=3))", 
                linewidth = 2.0, 
                color = colors_chi_e[plot_idx]
            )
            plot_idx += 1
        end
    end
    display(final_plot_chi_e) 

    println("\nAll calculations and plotting finished.")
end
