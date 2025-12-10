# Generates all nearest-neighbour bonds
function bond_list(Nx::Int64,Ny::Int64)
    N::Int64 = Nx*Ny
    Jbondlist::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,0) # Initialzing as empty array of []s
    Qbondlist::Vector{Vector{Vector{Int64}}} = Vector{Vector{Vector{Int64}}}(undef,0)
    if Nx==2 && Ny==2
        Jbondlist = [[0,1],[0,2],[1,3],[2,3]]         #Special Case for L=2
        Qbondlist = [[[0,1],[2,3]],[[0,2],[1,3]]]
    else
        for i::Int64 in 0:(N-1)
            x = div(i,Nx)
            
        # J-bonds
            xbond = [i,Nx*x+mod(i+1,Nx)]
            ybond = [i,mod(i+Nx,N)]
            push!(Jbondlist,xbond)
            push!(Jbondlist,ybond)
            
        # Q-bonds
            i1 = mod(i+Nx,N)
            x1 = div(i1,Nx)
        #-----------x bonds------------
            xij_bond = xbond
            xkl_bond = [i1,Nx*x1+mod(i1+1,Nx)]
            push!(Qbondlist,[xij_bond,xkl_bond])
        #-----------y bonds------------
            yij_bond = ybond
            i3 = Nx*x+mod(i+1,Nx)
            ykl_bond = [i3,mod(i3+Nx,N)]
            push!(Qbondlist,[yij_bond,ykl_bond])
        end
    end
    return (Jbondlist,Qbondlist)
end

# Generates lattice sites in order for x-translation
function lattice_x(Nx::Int64,Ny::Int64)
    site_list::Vector{Int64} = collect(Nx*Ny-1:-1:0)
    return site_list
end

# Generates lattice sites in order for y-translation
function lattice_y(Nx::Int64,Ny::Int64)
    site_list::Vector{Int64} = Vector{Int64}()
    A::Vector{Int64} = collect(Ny-1:-1:0)
    for i::Int64 in Nx-1:-1:0
        for j::Int64 in A
            k = i+Nx*j
            push!(site_list,k)
        end
    end
    return site_list
end

# Generates lattice sites and their (x,y) co-ordinates
function lattice_coords(Nx::Int64,Ny::Int64)
    coord::OrderedDict{Int64,Tuple{Int64,Int64}} = OrderedDict{Int64,Tuple{Int64,Int64}}()
    for i::Int64 in 0:(Nx*Ny-1)
        x = mod(i,Nx)
        y = div(i,Nx)
        coord[i] = (x,y)
    end
    return coord
end

# Performs x-translation
function x_translation(Nx::Int64,Ny::Int64,x_lattice::Vector{Int64},state::BitVector)
    state_c::BitVector = copy(state)
    temp = nothing         # Initialize to `nothing` to handle cases where it is not set
    
    for site::Int64 in x_lattice
        x = div(site,Nx)
        new_site = (Nx*x)+(site+1)%Nx
        
        if (site+1)%Nx == 0  # Takes care of boundary translation
            temp = state[site+1]
        else
            state_c[new_site+1] = state[site+1]
        end
        if site%Nx == 0     # Takes care of boundary translation
            state_c[site+1] = temp
        end
    end
    return state_c
end

# Performs y-translation
function y_translation(Nx::Int64,Ny::Int64,y_lattice::Vector{Int64},state::BitVector)
    N::Int64 = Nx*Ny
    state_c::BitVector = copy(state)
    temp = nothing      # Initialize to `nothing` to handle cases where it is not set
    
    for site::Int64 in y_lattice
        new_site = (site+Nx)%N
        
        if div(((site+Nx)%(Nx*Ny)),Nx) == 0   # Takes care of boundary translation
            temp = state[site+1]
        else
            state_c[new_site+1] = state[site+1]
        end
        if div(site,Nx) == 0              # Takes care of boundary translation
            state_c[site+1] = temp
        end
    end
    return state_c
end


# Finds representative state(smallest decimal) of a given state
function find_rep(state::BitVector,Nx::Int64,Ny::Int64,x_lattice::Vector{Int64},y_lattice::Vector{Int64},
                                                       x_tstate::BitVector,xy_tstate::BitVector)
    N::Int64 = Nx*Ny
    d::Int64 = state_dec(state)
    x_tstate .= copy(state)       #In-place modification to save memory
    lx::Int64,ly::Int64 = 0,0
    g::Int64 = 0

    for x::Int64 in 0:Nx-1
        for y::Int64 in 0:Ny-1
            t::Int64 = state_dec(x_tstate)
            
            if t < d
                d = t
                lx = x
                ly = y
                g = 0
            end
            
            # spin-invert the state
            t_inv::Int64 = invert_state(t,N)
            
            if t_inv < d
                d = t_inv
                lx = x
                ly = y
                g = 1
            end
                        
            # translate state in y direction
            xy_tstate .= y_translation(Nx,Ny,y_lattice,x_tstate) #In-place modification to save memory
            x_tstate .= xy_tstate
        end
        
        # translate state in x direction
        x_tstate .= x_translation(Nx,Ny,x_lattice,state)  
        state .= x_tstate                       
    end

    return (d,lx,ly,g)
end