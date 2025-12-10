# Checks compatibility of a state for given (kx,ky,z)
function checkstate(state::Int64,Nx::Int64,Ny::Int64,kx::Int64,ky::Int64,z::Int64,x_lattice::Vector{Int64},
                     y_lattice::Vector{Int64},bin_state::BitVector,x_tstate::BitVector,xy_tstate::BitVector)
    N::Int64 = Nx*Ny
    D::Int64 = -1
    phase::ComplexF64 = 0.0+0.0im
    state_set::Set{Int64} = Set{Int64}()    #Set removes dual entries and retains only unique states
    d::Int64 = state
    bin_state .= state_bin(state,N)   #Converts decimal state to binary array
    x_tstate .= copy(bin_state)       #In-place modification to save memory
    
    for x::Int64 in 0:(Nx-1)
        for y::Int64 in 0:(Ny-1)
            t::Int64 = state_dec(x_tstate)
            
            if t < d
                return (state,D,phase,state_set)    #Breaks the process if smaller decimal is found
            else
                push!(state_set,t)
                if t == d
                    phase += exp(-1im*2*pi*((kx*x)/Nx+(ky*y)/Ny))
                end
            end
            
            # spin-invert the state
            t_inv::Int64 = invert_state(t,N)
            
            if t_inv < d
                return (state,D,phase,state_set)  #Breaks the process if smaller decimal is found
            else
                push!(state_set,t_inv)
                if t_inv == d
                    phase += z*exp(-1im*2*pi*((kx*x)/Nx+(ky*y)/Ny))
                end
            end
            
            xy_tstate .= y_translation(Nx,Ny,y_lattice,x_tstate)
            x_tstate .= xy_tstate
        end
        
        x_tstate .= x_translation(Nx,Ny,x_lattice,bin_state)  #In-place modification to save memory
        bin_state .= x_tstate                        #In-place modification to save memory
    end
    
    if (abs(phase) > 1e-6)        #Condition on D for compatibility with kx and ky
        D = length(state_set)     #Number of different states produced by translations
    end
    return (state,D,phase,state_set)
end

# Generates basis list for a particular mz, k and z block
function gen_basis_mz_k_z(Nx::Int64,Ny::Int64,mz::Float64,kx::Int64,ky::Int64,z::Int64,
                          x_lattice::Vector{Int64},y_lattice::Vector{Int64})
    N::Int64 = Nx*Ny
    rep_state_list::Vector{Int64} = Vector{Int64}()
    D_list::Vector{Int64} = Vector{Int64}()
    Na_list::Vector{Float64} = Vector{Float64}()
    
    n_up::Int64 = Int64((N+(2*mz))/2)  #No. of up spins in lowest integer of an mz-block basis
    initial_num::Int64 = (2^n_up)-1    #Lowest integer of an mz-block basis
    
    #Allocate memeory once for in-place modification later
    bin_state::BitVector = x_tstate::BitVector = xy_tstate::BitVector = similar(state_bin(initial_num,N))
    
    M::Int64 = 0
    while true
        rep_state,D,phase,state_set = checkstate(initial_num,Nx,Ny,kx,ky,z,x_lattice,y_lattice,
                                                                    bin_state,x_tstate,xy_tstate)
        if D::Int64 >= 0
            M += 1
            push!(rep_state_list,rep_state)
            push!(D_list,D)
            push!(Na_list,D*abs(phase)^2)
        end
        
        #Generate next integer
        next_num::Int64 = generate_next_binary(initial_num,N)
        if next_num <= initial_num   # If no higher combination is found, break the loop
            break
        end
        initial_num = next_num
    end
    
    return (M,rep_state_list,D_list,Na_list) 
end