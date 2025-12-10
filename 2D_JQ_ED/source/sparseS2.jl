# Builds S^2 matrix
function Sparse_S2(Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,z::Int64,basis_list::Vector{Int64},
                    Na_list::Vector{Float64},x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                    basis_copy::BitVector,x_tstate::BitVector,xy_tstate::BitVector)

    H_size::Int64 = length(basis_list)
    data::Vector{ComplexF64} = Vector{ComplexF64}()   # Stores the matrix element values
    rows::Vector{Int64} = Vector{Int64}()             # Stores the row indices for sparse matrix
    cols::Vector{Int64} = Vector{Int64}()             # Stores the column indices for sparse matrix

    N::Int64 = Nx*Ny
    kxx::Float64 = 2*pi*kx/Nx
    kyy::Float64 = 2*pi*ky/Ny
    
    for dec_basis::Int64 in basis_list        
        ind1 = find_state_index(dec_basis,basis_list)   
        basis = state_bin(dec_basis,N)                 

        diag_elem::ComplexF64 = ((3/4)*N)+0.0im

        # Iterate over combinations of pairs (i,j)
        indices::Vector{Vector{Int64}} = collect(combinations(1:N,2))
        for (i::Int64,j::Int64) in indices
            s1 = basis[i]
            s2 = basis[j]

            if s1 == s2
                #*********** Diagonal terms *************
                diag_elem += 1/2
            else
                #*********** Diagonal terms *************
                diag_elem += -1/2

                #*********** Off-diagonal terms *********
                ii1 = flip(s1,s2)
                basis_copy .= copy(basis)     #In-place modification to reduce memory usage
                basis_copy[i] = ii1[1]
                basis_copy[j] = ii1[2]
                rep_basis,lx,ly,g = find_rep(basis_copy,Nx,Ny,x_lattice,y_lattice,x_tstate,xy_tstate) 
                ind2 = find_state_index(rep_basis,basis_list)

                if ind2 >= 0
                    N_a = Na_list[ind1]
                    N_b = Na_list[ind2]
                    mat_elem = 1*sqrt(N_b/N_a)*(z^g)*exp(-1im*((kxx*lx)+(kyy*ly)))

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

    S2 = sparse(rows,cols,data,H_size,H_size)
    elements_count = nnz(S2)  # Number of non-zero elements

    return S2,elements_count
end

# Calculates Total Spin(S) by counting degeneracy
function Total_spin(energy::Float64,mz1_energy::Vector{Float64},mz2_energy::Vector{Float64})
    match1::Bool = in(energy,mz1_energy)
    match2 = in(energy,mz2_energy)
    count::Int64 = 1
    if match1 == true
        count += 2
        if match2 == true
            count += 2
        end
    else
        if match2 == true
            count = -1
        end
    end
    return ((count-1)/2)::Float64
end