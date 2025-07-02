# testing hamiltonian construction
using ITensors, ITensorMPS
using Plots
using LaTeXStrings

gr()

# Second derivative function for non-uniformly spaced points
# (Robust even if points happen to be uniform)
function second_derivative_nonuniform(y::Vector{Float64}, x::Vector{Float64})
    n = length(y)
    if n < 3
        error("Input vectors must have at least 3 entries to calculate second derivative.")
    end
    if n != length(x)
        error("Input vectors y and x must be of the same length.")
    end

    d2y = zeros(n - 2) # Second derivative is calculated for interior points

    for i in 1:(n-2)
        x_prev = x[i]
        x_curr = x[i+1]
        x_next = x[i+2]

        y_prev = y[i]
        y_curr = y[i+1]
        y_next = y[i+2]

        if !(x_prev < x_curr && x_curr < x_next) # Ensure strictly increasing
             error("x values must be strictly increasing and distinct for the stencil. Problem at i=$i with x_prev=$x_prev, x_curr=$x_curr, x_next=$x_next")
        end

        h_prev = x_curr - x_prev
        h_next = x_next - x_curr

        if h_prev == 0.0 || h_next == 0.0
            error("Zero spacing detected in x values at i=$i. h_prev=$h_prev, h_next=$h_next")
        end
        
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
    
    # Fixed phi_xy
    phi_xy_fixed_rad = pi / 4.0 # 45 degrees

    # Scan range for theta_z
    theta_z_min_deg = -5.0
    theta_z_max_deg = 5.0
    num_theta_z_points = 51 # e.g., for a step of 0.1 degrees ( (5 - (-5))/0.1 + 1 = 101 )
    theta_z_scan_deg = range(theta_z_min_deg, stop=theta_z_max_deg, length=num_theta_z_points)
    theta_z_scan_rad = deg2rad.(theta_z_scan_deg) # Convert to radians for calculations

    h_values = [0.6, 0.65, 0.7] 

    #-----DATA STORAGE FOR LATER-----
    # Stores (list of theta_z_rad values, list of e0 values) for each h
    calculated_e0_data = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    susceptibility_results_theta_z = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    # The big loop
    for h in h_values
        print("Calculating for h = ", h, "\n")
        local_theta_z_rad_values = Float64[] # theta_z values in radians for this h 
        local_e0_values = Float64[]          # e0 (energy/N) values for this h

        for theta_z_current_rad in theta_z_scan_rad 
            current_theta_z_deg = rad2deg(theta_z_current_rad)
            print("  theta_z = ", round(current_theta_z_deg, digits=2)," deg (", round(theta_z_current_rad, digits=4), " rad)\n")
            
            halfway::Int = ceil(N/2)
            
            os_total_H = OpSum() 
            
            # Field components 'a', 'b', 'c'
            # phi_xy is fixed, theta_z_current_rad is from the scan
            a = cos(phi_xy_fixed_rad) * cos(theta_z_current_rad)
            b = sin(phi_xy_fixed_rad) * cos(theta_z_current_rad)
            c = sin(theta_z_current_rad)
    
            # Kitaev spin chain hamiltonian
            for i in 1:N-1
                if isodd(i)
                    os_total_H += J, "Sx", i, "Sx", i+1
                else
                    os_total_H += K, "Sy", i, "Sy", i+1
                end
            end
            
            # Magnetic field term for total Hamiltonian
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
        
            push!(local_theta_z_rad_values, theta_z_current_rad) 
            push!(local_e0_values, e0)
            
        end # end of theta_z_current_rad loop
        print("Finished calculations for h = ", h, "\n")
            
        calculated_e0_data[h] = (local_theta_z_rad_values, local_e0_values)
            
    end # end of h loop
    
    # -------- END OF DMRG LOOPS --------
    
    print("\nDMRG FINISHED. CALCULATING ENERGETIC SUSCEPTIBILITY (w.r.t. theta_z from e0)...\n")
    
    h_plot_values_sorted = sort(collect(keys(calculated_e0_data))) 
    
    for h_val in h_plot_values_sorted
        # theta_z_rad_scanned is already sorted ascending
        theta_z_rad_scanned, e0_for_h = calculated_e0_data[h_val] 
        
        if length(e0_for_h) < 3
            println("Not enough data points for h = $h_val to calculate second derivative w.r.t. theta_z. Skipping...")
            continue
        end

        d2_e0_dtheta_z2 = second_derivative_nonuniform(e0_for_h, theta_z_rad_scanned)

        susceptibility_chi_e_theta_z = -d2_e0_dtheta_z2

        # theta_z values for which the derivative was calculated
        theta_z_rad_for_derivative = theta_z_rad_scanned[2:end-1] 
        susceptibility_results_theta_z[h_val] = (theta_z_rad_for_derivative, susceptibility_chi_e_theta_z)
        print("Calculated energetic susceptibility (vs theta_z) for h = ", round(h_val, digits=3), "\n")
    end
        
    print("\nPLOTTING ENERGETIC SUSCEPTIBILITY (vs theta_z, derived from e0)...\n")
    
    colors_chi_e_theta_z = cgrad(:viridis, length(h_plot_values_sorted), categorical=true)

    # Determine common y-axis limits for better comparison if needed, or let Plots auto-scale
    # all_sus_values = vcat([sus for (_, sus_vector) in values(susceptibility_results_theta_z) for sus in sus_vector]...)
    # y_min = isempty(all_sus_values) ? 0 : minimum(all_sus_values)
    # y_max = isempty(all_sus_values) ? 1 : maximum(all_sus_values)
    # ylims_val = (floor(y_min, digits=2), ceil(y_max, digits=2))


    final_plot_chi_e_theta_z = plot(
        xlabel = L"\theta_z ", 
        ylabel = L"\chi^e_{\theta_z}", 
        title = "Energetic Susceptibility vs. Field Angle in z Direction", 
        legend = :topright
        # ylims = ylims_val # Uncomment to enforce common y-limits
    )
    
    plot_idx = 1
    for h_val in h_plot_values_sorted
        if haskey(susceptibility_results_theta_z, h_val)
            theta_z_rad_deriv, susceptibility = susceptibility_results_theta_z[h_val]
            theta_z_deg_deriv = rad2deg.(theta_z_rad_deriv) 

            plot!(
                final_plot_chi_e_theta_z,
                theta_z_deg_deriv,
                susceptibility,
                label = "h = $(round(h_val, digits=3))", 
                linewidth = 2.0, 
                color = colors_chi_e_theta_z[plot_idx]
            )
            plot_idx += 1
        end
    end
    display(final_plot_chi_e_theta_z) 

    println("\nAll calculations and plotting finished.")
end
