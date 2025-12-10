# Performs spin-flip
function flip(l0::Bool,l1::Bool)  #'0' and '1's have been represented using Boolean rather than Int64
    return [l1,l0]
end

# Finds index of a basis in a basis list
function find_state_index(state::Int64,basis_list::Vector{Int64})
    for (i::Int64,s::Int64) in enumerate(basis_list)
        if s == state
            return i
        end
    end
    return -1
end

# Reads bit at a position
function read_bit(i::Int64,bit::Int64)
    return ((i>>(bit-1))&1)::Int64
end

# Sets bit 1 to a position
function set_bit(i::Int64,bit::Int64)
    return (i|(1<<(bit-1)))::Int64
end

# Sets bit 0 to a position
function clear_bit(i::Int64,bit::Int64)
    return (i&(~(1<<(bit-1))))::Int64
end

# Invert the bits
function invert_state(state::Int64,N::Int64)
    mask = (1<<N)-1
    return ((~state)&mask)::Int64
end

# Converts a binary number(all digits in array format) to a decimal number
function state_dec(binary::BitVector)
    decimal::Int64 = 0
    for digit::Bool in binary
        decimal = (decimal*2)+convert(Int64,digit)
    end
    return decimal
end

# Converts a decimal number to binary number of size N
function state_bin(decimal::Int64,N::Int64)
    #BitVector is more memory efficient for storing binary entries (1 bit for each element) 
    bin_vector::BitVector = BitVector(undef,N)    
    for i::Int64 in 1:N
        bin_vector[N-i+1] = convert(Bool,(decimal&1))
        decimal >>= 1
    end
    return bin_vector
end

# Through bit shuffling, it generates binary of smallest decimals in ascending order
function generate_next_binary(num::Int64,N::Int64)
    i::Int64 = 0
    bit::Int64 = 1

    while bit::Int64 < N
        if read_bit(num,bit) == 1
            if read_bit(num,bit+1) == 1
                i += 1
            else
                for j::Int64 in 1:i
                    num = set_bit(num,j)
                end
                
                for j::Int64 in (i+1):(bit)
                    num = clear_bit(num,j)
                end
                
                num = set_bit(num,bit+1)
                break
            end
        end

        bit += 1
    end
    return num
end

# Generates all possible mz values
function find_mz(Nx::Int64, Ny::Int64)
    N::Int64 = Nx*Ny
    mz_list::Vector{Float64} = collect(-N/2:N/2)  # Range from -N/2 to N/2, inclusive
    return sort(mz_list)            # Convert set to a sorted list (array in Julia)
end

# Combines energies from both inversion sectors and extracts lowest 10 energies with index of corresponding z-sector
function lowest_k_energies(energies1::Vector{NamedTuple{(:energy,:kx,:ky,:z),Tuple{Float64,Int64,Int64,Int64}}},
                   energies2::Vector{NamedTuple{(:energy,:kx,:ky,:z),Tuple{Float64,Int64,Int64,Int64}}},
                   z1::Int64,z2::Int64,k::Int64)
    
    # Combine the energies from both sectors
    combined_energies::Vector{NamedTuple{(:energy,:kx,:ky,:z),Tuple{Float64,Int64,Int64,Int64}}} =                                                                                                  vcat(energies1,energies2)
    
    # Sort energies based on the combined energies
    lowest_energies::Vector{NamedTuple{(:energy,:kx,:ky,:z),Tuple{Float64,Int64,Int64,Int64}}} =                                                                   sort(combined_energies, by = x -> x.energy)[1:min(k,end)]
    
    # Initialize lists to store indices for z = 1 and z = -1
    indices_z1::Vector{Int64} = Vector{Int64}()      # List for indices with z = z1
    indices_z2::Vector{Int64} = Vector{Int64}()      # List for indices with z = z2
    
    # Iterate over the lowest 10 energies and store the indices
    for (i::Int64,energy_tuple) in enumerate(lowest_energies)
        if energy_tuple.z == z1
            push!(indices_z1,i) 
        elseif energy_tuple.z == z2
            push!(indices_z2,i) 
        end
    end
    
    return lowest_energies,indices_z1,indices_z2
end

# Extracts eigenvectors according to corresponding inversion sector
function extract_eigenvectors(evecs_z1::Matrix{ComplexF64},indices_z1::Vector{Int64}, 
                              evecs_z2::Matrix{ComplexF64},indices_z2::Vector{Int64})
    # Initialize an OrderedDict to store the extracted eigenvectors for both sectors
    eigenvecs_combined::OrderedDict{Int64,Vector{ComplexF64}} = OrderedDict{Int64,Vector{ComplexF64}}()

    # Extract eigenvectors for z1 and z2 efficiently using views
    # Store z1 eigenvectors in the combined dictionary
    for (i::Int64,idx::Int64) in enumerate(indices_z1)
        eigenvecs_combined[idx] = evecs_z1[:,i]  # Store the i-th column of evecs1
    end

    # Store z2 eigenvectors in the same combined dictionary
    for (i::Int64,idx::Int64) in enumerate(indices_z2)
        eigenvecs_combined[idx] = evecs_z2[:,i]  # Store the i-th column of evecs2
    end

    return eigenvecs_combined
end