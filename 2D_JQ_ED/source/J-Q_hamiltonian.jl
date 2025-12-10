#################################### J part of Hamiltonian ##########################################

function JQ_Ham_J(J::Float64,Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,z::Int64,basis_list::Vector{Int64},                          Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                   Jbonds::Vector{Vector{Int64}},use_sparse::Bool,basis_copy::BitVector,
                   x_tstate::BitVector,xy_tstate::BitVector)
    
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{ComplexF64} = Vector{ComplexF64}()  # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{ComplexF64} = zeros(ComplexF64,H_size,H_size) # Complex matrix for non-sparse case
    end
    
    N::Int64 = Nx*Ny
    kxx::Float64 = 2*pi*kx/Nx
    kyy::Float64 = 2*pi*ky/Ny
    
    for dec_basis::Int64 in basis_list
        diag_energy::ComplexF64 = 0.0+0.0im
        ind1::Int64 = find_state_index(dec_basis,basis_list)     #Function that returns index of state
        basis::BitVector = state_bin(dec_basis,N)                #Converts Integer to binary array      
        
        for bond in Jbonds
            s1 = basis[bond[1]+1]     #Julia uses 1-based indexing
            s2 = basis[bond[2]+1]
            
            if s1 != s2
            # *********** J-diagonal terms *************
                diag_energy += -J/2
                
            # *********** J-off diagonal terms *********
                ii1 = flip(s1,s2)  
                basis_copy .= copy(basis)       #In-place modification
                basis_copy[bond[1]+1] = ii1[1]
                basis_copy[bond[2]+1] = ii1[2]
                rep_basis,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate)
                ind2 = find_state_index(rep_basis,basis_list)
                if ind2 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind2]
                    mat_elem = (J/2)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                    if use_sparse != true
                        H_dense[ind1,ind2] += mat_elem
                    else
                        push!(data,mat_elem)
                        push!(rows,ind1)
                        push!(cols,ind2)
                    end
                end
            end
        end
        
        if use_sparse != true
            H_dense[ind1,ind1] += diag_energy
        else
            push!(data,diag_energy)
            push!(rows,ind1)
            push!(cols,ind1)
        end
    end
                            
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size,basis_list,Na_list)
    else
        return (H_dense,H_size,basis_list,Na_list)
    end
end           

#################################### Q part of Hamiltonian ##########################################

function JQ_Ham_Q(Q::Float64,Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,z::Int64,basis_list::Vector{Int64},
                    Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                    Qbonds::Vector{Vector{Vector{Int64}}},use_sparse::Bool,basis_copy::BitVector,
                    flipped_basis::BitVector,x_tstate::BitVector,xy_tstate::BitVector)
    
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{ComplexF64} = Vector{ComplexF64}()  # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{ComplexF64} = zeros(ComplexF64,H_size,H_size) # Complex matrix for non-sparse case
    end
    
    N::Int64 = Nx*Ny
    kxx::Float64 = 2*pi*kx/Nx
    kyy::Float64 = 2*pi*ky/Ny
    
    for dec_basis::Int64 in basis_list
        diag_energy::ComplexF64 = 0.0+0.0im
        ind1::Int64 = find_state_index(dec_basis,basis_list)  #Function that returns index of state
        basis::BitVector = state_bin(dec_basis,N)                 #Converts Integer to binary array        
        
        for bij_bkl in Qbonds
            b_ij = bij_bkl[1]     # ij bond
            b_kl = bij_bkl[2]     # kl bond

            s1_ij = basis[b_ij[1]+1]
            s2_ij = basis[b_ij[2]+1]

            s1_kl = basis[b_kl[1]+1]
            s2_kl = basis[b_kl[2]+1]
            
            if (s1_ij != s2_ij) && (s1_kl != s2_kl)
            
        #*********** Q-(d(ij)_d(kl)) term *************
                diag_energy += -Q/4
                
        #*********** Q-(d(ij)_od(kl)) term *************
                ii1 = flip(s1_kl,s2_kl)
                basis_copy .= copy(basis)       #In-place modification
                basis_copy[b_kl[1]+1] = ii1[1]
                basis_copy[b_kl[2]+1] = ii1[2]
                rep_basis1,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
                ind2 = find_state_index(rep_basis1,basis_list)
                if ind2 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind2]
                    mat_elem = (Q/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                    if use_sparse != true
                        H_dense[ind1,ind2] += mat_elem
                    else
                        push!(data,mat_elem)
                        push!(rows,ind1)
                        push!(cols,ind2)
                    end
                end
                
        # *********** Q-(od(ij)_d(kl)) term *************
                ii2 = flip(s1_ij,s2_ij)
                basis_copy .= copy(basis)    #In-place modification
                basis_copy[b_ij[1]+1] = ii2[1]
                basis_copy[b_ij[2]+1] = ii2[2]
                rep_basis2,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate)
                ind3 = find_state_index(rep_basis2,basis_list)
                if ind3 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind3]
                    mat_elem1 = (Q/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                    if use_sparse != true
                        H_dense[ind1,ind3] += mat_elem1
                    else
                        push!(data,mat_elem1)
                        push!(rows,ind1)
                        push!(cols,ind3)
                    end
                end
                
        # *********** Q-(od(ij)_od(kl)) term *************
                ii3 = flip(s1_ij,s2_ij)
                basis_copy .= copy(basis)     #In-place modification
                basis_copy[b_ij[1]+1] = ii3[1]
                basis_copy[b_ij[2]+1] = ii3[2]
                
                ii4 = flip(s1_kl,s2_kl)
                flipped_basis .= copy(basis_copy)    #In-place modification
                flipped_basis[b_kl[1]+1] = ii4[1]
                flipped_basis[b_kl[2]+1] = ii4[2]
                rep_basis3,lx,ly,g = find_rep(flipped_basis,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
                ind4 = find_state_index(rep_basis3,basis_list)
                if ind4 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind4]
                    mat_elem2 = -(Q/4)*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))
                    if use_sparse != true
                        H_dense[ind1,ind4] += mat_elem2
                    else
                        push!(data,mat_elem2)
                        push!(rows,ind1)
                        push!(cols,ind4)
                    end
                end
            end
        end
        
        if use_sparse != true
            H_dense[ind1,ind1] += diag_energy
        else
            push!(data,diag_energy)
            push!(rows,ind1)
            push!(cols,ind1)
        end
    end
                            
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size,basis_list,Na_list)
    else
        return (H_dense,H_size,basis_list,Na_list)
    end
end           