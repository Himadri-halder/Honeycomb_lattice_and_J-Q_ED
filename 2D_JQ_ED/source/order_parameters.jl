#************************* Staggered Spin-Spin Correlation***********************
function Stagg_spin_corr(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,basis_list::Vector{Int64},
                            Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                            site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},basis_copy::BitVector,
                            x_tstate::BitVector,xy_tstate::BitVector)
    
    H_size = length(basis_list)
    data,rows,cols = Vector{ComplexF64}(),Vector{Int64}(),Vector{Int64}()
    
    N = Nx*Ny
    kxx = 2*pi*kx/Nx
    kyy = 2*pi*ky/Ny
    
    for dec_basis in basis_list
        ind1 = find_state_index(dec_basis,basis_list)   
        basis = state_bin(dec_basis,N) 
        
        diag_elem = 0.0+0.0im
        s0 = basis[1]  
        for r in keys(site_coords)
            sr = basis[r+1]
            x,y = site_coords[r][1], site_coords[r][2]   # Access the coordinates
            stagg_phase = (-1)^(x+y)

            if s0 == sr
                #*********** Diagonal terms *************
                diag_elem += stagg_phase/4
            else
                #*********** Diagonal terms *************
                diag_elem += -stagg_phase/4
                
                #*********** Off-diagonal terms *************
                ii1 = flip(s0,sr)  
                basis_copy .= copy(basis)     #In-place modification
                basis_copy[1] = ii1[1]
                basis_copy[r+1] = ii1[2]
                rep_basis,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
                ind2 = find_state_index(rep_basis,basis_list)
                if ind2 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind2]
                    mat_elem = stagg_phase*(1/2)*(z^g)*sqrt(N_b/N_a)*exp(-1im*((kxx*lx)+(kyy*ly)))
                    
                    push!(data,mat_elem)
                    push!(rows,ind1)
                    push!(cols,ind2)
                end
            end
        end    
        push!(data,diag_elem)
        push!(rows,ind1)
        push!(cols,ind1)
    end
        
    data_per_site = data./N     #Per site
    Mz_matrix = sparse(rows,cols,data_per_site,H_size,H_size)
        
    return Mz_matrix
end

#************************Staggered x-directed Dimer-Dimer Correlation***********************
function Stagg_x_dimer_corr(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,basis_list::Vector{Int64},
                            Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                            site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},basis_copy::BitVector,
                            x_tstate::BitVector,xy_tstate::BitVector)
    
    H_size = length(basis_list)
    data,rows,cols = Vector{ComplexF64}(),Vector{Int64}(),Vector{Int64}()
    
    N = Nx*Ny
    kxx = 2*pi*kx/Nx
    kyy = 2*pi*ky/Ny
    
    for dec_basis in basis_list
        ind1 = find_state_index(dec_basis,basis_list)   
        basis = state_bin(dec_basis,N) 
        
        diag_elem = 0.0+0.0im
        s0 = basis[1]
        s1 = basis[2]
        
        for r in keys(site_coords)
            x = site_coords[r][1]
            stagg_phase = (-1)^x
            
            a = div(r,Nx)
            r1_site = Nx*a+(r+1)%Nx
            sr = basis[r+1]        
            sr1 = basis[r1_site+1] 

            #***************** S-(d(0,x)_d(r,r+x)) term ****************
            diag_elem += stagg_phase*(1/4)*(1/4)*((-1)^(s0+s1+sr+sr1))
            
            #***************** S-(d(0,x)_od(r,r+x)) term ****************
            ii1 = flip(sr,sr1)  
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[r+1] = ii1[1]
            basis_copy[r1_site+1] = ii1[2]
            rep_basis1,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind2 = find_state_index(rep_basis1,basis_list)
            if ind2 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind2]
                mat_elem1 = stagg_phase*(1/8)*((-1)^(s0+s1))*(z^g)*sqrt(N_b/N_a)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem1)
                push!(rows,ind1)
                push!(cols,ind2)
            end
            
            #***************** S-(od(0,x)_d(r,r+x)) term ****************
            ii2 = flip(s0,s1)  
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[1] = ii2[1]
            basis_copy[2] = ii2[2]
            rep_basis2,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind3 = find_state_index(rep_basis2,basis_list)
            if ind3 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind3]
                mat_elem2 = 
                           stagg_phase*(1/8)*((-1)^(sr+sr1))*(z^g)*sqrt(N_b/N_a)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem2)
                push!(rows,ind1)
                push!(cols,ind3)
            end
            
            #***************** S-(od(0,x)_od(r,r+x)) term ****************
            ii3 = flip(sr,sr1)
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[r+1] = ii3[1]
            basis_copy[r1_site+1] = ii3[2]
            
            ii4 = flip(s0,s1)
            flipped_basis .= copy(basis_copy)    #In-place modification
            flipped_basis[1] = ii4[1]
            flipped_basis[2] = ii4[2]
            rep_basis3,lx,ly = find_rep(flipped_basis,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind4 = find_state_index(rep_basis3,basis_list)
            if ind4 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind4]
                mat_elem3 = stagg_phase*(1/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem3)
                push!(rows,ind1)
                push!(cols,ind4)
            end
        end
        push!(data,diag_elem)
        push!(rows,ind1)
        push!(cols,ind1)
    end
        
    data_per_site = data./N     #Per site
    matrix = sparse(rows,cols,data_per_site,H_size,H_size)
    
    return matrix
end

#************************Staggered y-directed Dimer-Dimer Correlation***********************
function Stagg_y_dimer_corr(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,basis_list::Vector{Int64},
                            Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                            site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},basis_copy::BitVector,
                            x_tstate::BitVector,xy_tstate::BitVector)
    
    H_size = length(basis_list)
    data,rows,cols = Vector{ComplexF64}(),Vector{Int64}(),Vector{Int64}()
    
    N = Nx*Ny
    kxx = 2*pi*kx/Nx
    kyy = 2*pi*ky/Ny
    
    for dec_basis in basis_list
        ind1 = find_state_index(dec_basis,basis_list)   
        basis = state_bin(dec_basis,N)
        
        diag_elem = 0.0+0.0im
        s0 = basis[1]
        next_site = (0+Nx)%N
        s1 = basis[next_site+1]
        
        for r in keys(site_coords)
            y = site_coords[r][2]
            stagg_phase = (-1)^y
            
            r1_site = (r+Nx)%N
            sr = basis[r+1]        
            sr1 = basis[r1_site+1] 

            #***************** S-(d(0,x)_d(r,r+x)) term ****************
            diag_elem += stagg_phase*(1/4)*(1/4)*((-1)^(s0+s1+sr+sr1))
            
            #***************** S-(d(0,x)_od(r,r+x)) term ****************
            ii1 = flip(sr,sr1)  
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[r+1] = ii1[1]
            basis_copy[r1_site+1] = ii1[2]
            rep_basis1,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind2 = find_state_index(rep_basis1,basis_list)
            if ind2 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind2]
                mat_elem1 = stagg_phase*(1/8)*((-1)^(s0+s1))*(z^g)*sqrt(N_b/N_a)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem1)
                push!(rows,ind1)
                push!(cols,ind2)
            end
            
            #***************** S-(od(0,x)_d(r,r+x)) term ****************
            ii2 = flip(s0,s1)  
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[1] = ii2[1]
            basis_copy[2] = ii2[2]
            rep_basis2,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind3 = find_state_index(rep_basis2,basis_list)
            if ind3 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind3]
                mat_elem2 = 
                           stagg_phase*(1/8)*((-1)^(sr+sr1))*(z^g)*sqrt(N_b/N_a)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem2)
                push!(rows,ind1)
                push!(cols,ind3)
            end
            
            #***************** S-(od(0,x)_od(r,r+x)) term ****************
            ii3 = flip(sr,sr1)
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[r+1] = ii3[1]
            basis_copy[r1_site+1] = ii3[2]
            
            ii4 = flip(s0,s1)
            flipped_basis .= copy(basis_copy)    #In-place modification
            flipped_basis[1] = ii4[1]
            flipped_basis[2] = ii4[2]
            rep_basis3,lx,ly = find_rep(flipped_basis,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind4 = find_state_index(rep_basis3,basis_list)
            if ind4 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind4]
                mat_elem3 = stagg_phase*(1/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem3)
                push!(rows,ind1)
                push!(cols,ind4)
            end
        end
        push!(data,diag_elem)
        push!(rows,ind1)
        push!(cols,ind1)
    end
        
    data_per_site = data./N     #Per site
    matrix = sparse(rows,cols,data_per_site,H_size,H_size)
    
    return matrix
end

#****************************** Neel Structure Factor ****************************
function Neel_structure(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,basis_list::Vector{Int64},
                        Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                        site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},basis_copy::BitVector,
                        x_tstate::BitVector,xy_tstate::BitVector,temp::Vector{ComplexF64},
                        eigenvector::Vector{ComplexF64},qx::Real,qy::Real)
    
    H_size = length(basis_list)
    
    N = Nx*Ny
    kxx = 2*pi*kx/Nx
    kyy = 2*pi*ky/Ny
    
    value = 0.0+0.0im
    for r in keys(site_coords)
        x,y = site_coords[r][1], site_coords[r][2] 
        
        data,rows,cols = Vector{ComplexF64}(),Vector{Int64}(),Vector{Int64}()
        for dec_basis in basis_list
            ind1 = find_state_index(dec_basis,basis_list)  #Function that returns index of state
            basis = state_bin(dec_basis,N)                 #Converts Integer to binary array
            
            s0 = basis[1]
            sr = basis[r+1]
            diag_elem = 0.0+0.0im
                     
            if s0 == sr
                #*********** Diagonal terms *************
                diag_elem = 1/4
            else
                #*********** Diagonal terms *************
                diag_elem = -1/4
                
                #*********** Off-diagonal terms *************
                ii1 = flip(s0,sr)  
                basis_copy .= copy(basis)       #In-place modification
                basis_copy[1] = ii1[1]
                basis_copy[r+1] = ii1[2]
                rep_basis1,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
                ind2 = find_state_index(rep_basis1,basis_list)
                if ind2 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind2]
                    mat_elem = (1/2)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                    
                    push!(data,mat_elem)
                    push!(rows,ind1)
                    push!(cols,ind2)
                end
            end
            push!(data,diag_elem)
            push!(rows,ind1)
            push!(cols,ind1)
        end 
        
        data_per_site = data./N
        corr_matrix = sparse(rows,cols,data_per_site,H_size,H_size)
        mul!(temp,corr_matrix,eigenvector)      #In-place modification
        corr_avg = dot(eigenvector,temp)
        
        value += exp(-1im*((qx*x)+(qy*y)))*corr_avg
    end
    return value
end

#****************************** VBS (pi,0) Structure Factor ****************************
function VBS_p0_structure(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,basis_list::Vector{Int64},
                           Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                           site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},basis_copy::BitVector,
                           flipped_basis::BitVector,x_tstate::BitVector,xy_tstate::BitVector,
                           temp::Vector{ComplexF64},eigenvector::Vector{ComplexF64},qx::Real,qy::Real)
    
    H_size = length(basis_list)
    
    N = Nx*Ny
    kxx = 2*pi*kx/Nx
    kyy = 2*pi*ky/Ny
    
    value = 0.0+0.0im
    for r in keys(site_coords)
        x,y = site_coords[r][1], site_coords[r][2] 
        
        data,rows,cols = Vector{ComplexF64}(),Vector{Int64}(),Vector{Int64}()
        for dec_basis in basis_list
            ind1 = find_state_index(dec_basis,basis_list)  
            basis = state_bin(dec_basis,N)                 
            
            s0 = basis[1]
            s1 = basis[2]
            
            a = div(r,Nx)
            r1_site = Nx*a+(r+1)%Nx
            sr = basis[r+1]        
            sr1 = basis[r1_site+1]
            diag_elem = 0.0+0.0im
            
            #***************** S-(d(0,1)_d(r,r+1)) term ****************
            diag_elem = (1/4)*(1/4)*((-1)^(s0+s1+sr+sr1))
            
            #***************** S-(d(0,1)_od(r,r+1)) term ****************
            ii1 = flip(sr,sr1)  
            basis_copy .= copy(basis)       #In-place modification
            basis_copy[r+1] = ii1[1]
            basis_copy[r1_site+1] = ii1[2]
            rep_basis1,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind2 = find_state_index(rep_basis1,basis_list)
            if ind2 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind2]
                mat_elem1 = (1/8)*((-1)^(s0+s1))*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem1)
                push!(rows,ind1)
                push!(cols,ind2)
            end
            
            #***************** S-(od(0,1)_d(r,r+1)) term ****************
            ii2 = flip(s0,s1)  
            basis_copy .= copy(basis)    #In-place modification
            basis_copy[1] = ii2[1]
            basis_copy[2] = ii2[2]
            rep_basis2,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind3 = find_state_index(rep_basis2,basis_list)
            if ind3 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind3]
                mat_elem2 = (1/8)*((-1)^(sr+sr1))*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem2)
                push!(rows,ind1)
                push!(cols,ind3)
            end
            
            #***************** S-(od(0,1)_od(r,r+1)) term ****************
            ii3 = flip(sr,sr1)
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[r+1] = ii3[1]
            basis_copy[r1_site+1] = ii3[2]
            
            ii4 = flip(s0,s1)
            flipped_basis .= copy(basis_copy)    #In-place modification
            flipped_basis[1] = ii4[1]
            flipped_basis[2] = ii4[2]
            rep_basis3,lx,ly,g = find_rep(flipped_basis,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind4 = find_state_index(rep_basis3,basis_list)
            if ind4 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind4]
                mat_elem3 = (1/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem3)
                push!(rows,ind1)
                push!(cols,ind4)
            end
            push!(data,diag_elem)
            push!(rows,ind1)
            push!(cols,ind1)
        end
        
        data_per_site = data./N
        corr_matrix = sparse(rows,cols,data_per_site,H_size,H_size)
        mul!(temp,corr_matrix,eigenvector)      #In-place modification
        corr_avg = dot(eigenvector,temp)
        
        value += exp(-1im*((qx*x)+(qy*y)))*corr_avg
    end
    return value
end

#****************************** VBS (0,pi) Structure Factor ****************************
function VBS_0p_structure(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,basis_list::Vector{Int64},
                           Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                           site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},basis_copy::BitVector,
                           flipped_basis::BitVector,x_tstate::BitVector,xy_tstate::BitVector,
                           temp::Vector{ComplexF64},eigenvector::Vector{ComplexF64},qx::Real,qy::Real)
    
    H_size = length(basis_list)
    
    N = Nx*Ny
    kxx = 2*pi*kx/Nx
    kyy = 2*pi*ky/Ny
    
    value = 0.0+0.0im
    for r in keys(site_coords)
        x,y = site_coords[r][1], site_coords[r][2] 
        
        data,rows,cols = Vector{ComplexF64}(),Vector{Int64}(),Vector{Int64}()
        for dec_basis in basis_list
            ind1 = find_state_index(dec_basis,basis_list)  
            basis = state_bin(dec_basis,N)
            
            s0 = basis[1]
            next_site = (0+Nx)%N
            s1 = basis[next_site+1]
            
            r1_site = (r+Nx)%N
            sr = basis[r+1]        
            sr1 = basis[r1_site+1] 
            diag_elem = 0.0+0.0im
            
            #***************** S-(d(0,1)_d(r,r+1)) term ****************
            diag_elem = (1/4)*(1/4)*((-1)^(s0+s1+sr+sr1))
            
            #***************** S-(d(0,1)_od(r,r+1)) term ****************
            ii1 = flip(sr,sr1)  
            basis_copy .= copy(basis)       #In-place modification
            basis_copy[r+1] = ii1[1]
            basis_copy[r1_site+1] = ii1[2]
            rep_basis1,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind2 = find_state_index(rep_basis1,basis_list)
            if ind2 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind2]
                mat_elem1 = (1/8)*((-1)^(s0+s1))*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem1)
                push!(rows,ind1)
                push!(cols,ind2)
            end
            
            #***************** S-(od(0,1)_d(r,r+1)) term ****************
            ii2 = flip(s0,s1)  
            basis_copy .= copy(basis)    #In-place modification
            basis_copy[1] = ii2[1]
            basis_copy[2] = ii2[2]
            rep_basis2,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind3 = find_state_index(rep_basis2,basis_list)
            if ind3 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind3]
                mat_elem2 = (1/8)*((-1)^(sr+sr1))*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem2)
                push!(rows,ind1)
                push!(cols,ind3)
            end
            
            #***************** S-(od(0,1)_od(r,r+1)) term ****************
            ii3 = flip(sr,sr1)
            basis_copy .= copy(basis)     #In-place modification
            basis_copy[r+1] = ii3[1]
            basis_copy[r1_site+1] = ii3[2]
            
            ii4 = flip(s0,s1)
            flipped_basis .= copy(basis_copy)    #In-place modification
            flipped_basis[1] = ii4[1]
            flipped_basis[2] = ii4[2]
            rep_basis3,lx,ly = find_rep(flipped_basis,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
            ind4 = find_state_index(rep_basis3,basis_list)
            if ind4 >= 0
                N_a = Na_list[ind1]
                N_b = Na_list[ind4]
                mat_elem3 = (1/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                
                push!(data,mat_elem3)
                push!(rows,ind1)
                push!(cols,ind4)
            end
            push!(data,diag_elem)
            push!(rows,ind1)
            push!(cols,ind1)
        end
        
        data_per_site = data./N
        corr_matrix = sparse(rows,cols,data_per_site,H_size,H_size)
        mul!(temp,corr_matrix,eigenvector)      #In-place modification
        corr_avg = dot(eigenvector,temp)
        
        value += exp(-1im*((qx*x)+(qy*y)))*corr_avg
    end
    return value
end