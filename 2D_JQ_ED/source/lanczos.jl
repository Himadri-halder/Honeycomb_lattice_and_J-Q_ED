function Act_H(H::SparseMatrixCSC{ComplexF64,Int64},vect::Vector{ComplexF64},outvect::Vector{ComplexF64})
    mul!(outvect,H,vect)  # Perform matrix-vector multiplication in-place by just modifying outvect
    return outvect::Vector{ComplexF64}
end

function overlap(vect1::Vector{ComplexF64},vect2::Vector{ComplexF64})  # both vector are of size M
    v::ComplexF64 = dot(vect1,vect2)   # dot() preserves the complex part 
    return v
end

#converts eigenvector from lanczos basis into original basis
function lanczos_eigenvec_mapper(lancbasis::Vector{Vector{ComplexF64}},lancvec::Vector{ComplexF64},
                                  M::Int64,lmda::Int64)
    original_coeff::Vector{ComplexF64} = Vector{ComplexF64}(undef,M)  # Pre-allocate the array for efficiency
    
    for i::Int64 in 1:M
        temp::ComplexF64 = 0.0+0.0im
        for j::Int64 in 1:lmda
            temp += lancvec[j]*lancbasis[j][i]
        end
        original_coeff[i] = temp
    end
    return original_coeff
end

#Creating Lanczos basis (with full reorthogonalization)
function Lanczos_reortho(M::Int64,lmda::Int64,H::SparseMatrixCSC{ComplexF64,Int64})
    # M: size of Lanczos vector
    # lmda: size of Lanczos basis
    # H: original Hamiltonian (of size M)
    # Initialization
    phi::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}(undef,lmda)  # Array of Lanczos vectors
    Trid_Ham::Matrix{ComplexF64} = zeros(ComplexF64,lmda,lmda)  # Tridiagonal Hamiltonian matrix of size lmda
    Random.seed!(10)
    initial_vect::Vector{ComplexF64} = (rand(M).+0.00001).+1im.*(rand(M).+0.00001)
    # A small offset of 10^-5 is set to prevent norm being 0 and dividing by 0 situation
    norm_initial::Float64 = norm(initial_vect)
    phi[1] = initial_vect/norm_initial
    
    temp_phi = similar(phi[1])               #create a duplicate array
    temp_phi = Act_H(H,phi[1],temp_phi)      #perform in-place multiplication 
    Trid_Ham[1,1] = overlap(phi[1],temp_phi)
    @.temp_phi -= (Trid_Ham[1,1])*phi[1]     #in-place operation to save memory
    norm_temp = norm(temp_phi)
    Trid_Ham[1,2] = Trid_Ham[2,1] = norm_temp
    phi[2] = temp_phi/norm_temp
    
    # Now the iterative part
    last_vect = phi[2] 
    last1_vect = phi[1]
    tempvect = similar(phi[1])               #create a duplicate array
    for n::Int64 in 2:(lmda-1)  
        temp_phi = Act_H(H,last_vect,temp_phi)   #Takes temp-phi and modifies it in-place in output
        Trid_Ham[n,n] = overlap(last_vect,temp_phi)
        @.temp_phi -= Trid_Ham[n,n]*last_vect + Trid_Ham[n-1,n]*last1_vect  #in-place modification
        Trid_Ham[n+1,n] = Trid_Ham[n,n+1] = norm(temp_phi)
        phi[n+1] = temp_phi/norm(temp_phi)

        # Extra code for reorthogonalization of Lanczos vectors
        for n1::Int64 in 1:n  
            temp = overlap(phi[n1],phi[n+1])
            @.tempvect = phi[n+1]-(temp*phi[n1])   #in-place modification to reduce memory usage 
            norm_temp = norm(tempvect)
            @.phi[n+1] = tempvect/norm_temp        #in-place modification to reduce memory usage
        end
        last1_vect = last_vect
        last_vect = phi[n+1]
    end
    
    Trid_Ham[lmda,lmda] = overlap(phi[lmda],Act_H(H,phi[lmda],temp_phi))
    
    return Trid_Ham,phi
end

function Lanczos_ED(H::SparseMatrixCSC{Complex{Float64},Int64},M::Int64,lmda::Int64,k::Int64)
    # Hamiltonian
    # M: Size of Hamiltonian
    # lmda: No. of iterations
    # k: No. of energies wanted
    
    Trid_Ham,phi = Lanczos_reortho(M,lmda,H)
    V,W = eigen(Trid_Ham)
    eigenvec::Matrix{ComplexF64} = Matrix{ComplexF64}(undef,M,k)
    for i::Int64 in 1:k
        eigenvec[:,i] = lanczos_eigenvec_mapper(phi,W[:,i],M,lmda)
    end
    
    # Return the first k eigenvalues and the corresponding eigenvectors
    return V[1:k],eigenvec
end