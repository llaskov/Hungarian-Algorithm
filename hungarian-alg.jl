# Transformation $T_{uv}$
# Calculate the reduction defined in Lemma 4.5 (Mathematical models for optimization)
# in_mtrx       input / output matrix
function reduce!(in_mtrx)

    # reduce each row by the minimum element in the row
    for i in eachindex(in_mtrx[:, 1])
        in_mtrx[i, :] .-= minimum(in_mtrx[i, :])
    end

    # in the result, reduce each column by the minimum element in the column
    for j in eachindex(in_mtrx[1, :])
        in_mtrx[:, j] .-= minimum(in_mtrx[:, j])
    end
end

# Transformation $O_{X_1 Y_1}$
# in_mtrx       input / output matrix
# x_set         the $X_1$ set is a list of non negative integers
# y_set         the $Y_1$ set is a list of non negative integers
function zeroxy!(in_mtrx, x_set, y_set)
    # complement of the $Y_1$ set
    y_cmp = setdiff(collect(1:size(in_mtrx, 1)), y_set)
    
    # operation 1
    # find minimal element with indexes determined by $X_1$ and the complement of the $Y_1$
    d = in_mtrx[x_set[1], y_cmp[1]]
    for j in y_cmp, i in x_set        
        if in_mtrx[i, j] < d
            d = in_mtrx[i, j]           
        end
    end

    # operation 2
    # for each index i stored in $X_1$, subtract d from the i-th row of the matrix
    for i in x_set
        in_mtrx[i, :] .-= d   
    end

    # operation 3
    # for each index j stored in $Y_1$, add d from the j-th column of the matrix
    for j in y_set
        in_mtrx[:, j] .+= d   
    end
end

# Calculate the adjacency list of vertices of vertices connected with zero weight edges
# in_mtrx       input adjacency matrix
# return        adjacency list
function zeroadj(in_mtrx)
    result = fill(Vector{Int64}([]), size(in_mtrx, 1))
    
    for j in eachindex(in_mtrx[1, :]), i in eachindex(in_mtrx[:, 1]) 
        if in_mtrx[i, j] == 0
            if size(result[i], 1) == 0
                # the list is empty, then create a list with single element
                result[i] = [j]
            else
                # the list is not empty, then push the element
                push!(result[i], j)
            end
        end
    end

    return result
end

# DFS for a M-altering chain
# adj_lst_zero  input zero weight edges adjacency list
# my            input adjacency vector matching from set Y
# x_sta         input vertex from set X to start the search from
# px            input/output parents of the traversed vertices from X
# py            input/output parents of the traversed vertices from Y
# tr_x          input/output traversed vertices from X
# tr_y          input/output traversed vertices from Y
# found         return true if found and false otherwise 
# y_end         return the last vertex from the chain in the Y set
function malter!(adj_lst_zero, my, x_sta, px, py, tr_x, tr_y)
    found = false
    y_end = 0
    j = 1
    while j <= size(adj_lst_zero[x_sta], 1) && found == false
        y_end = adj_lst_zero[x_sta][j]
        if py[y_end] == 0
            py[y_end] = x_sta
            push!(tr_y, y_end)
            x_inr = my[y_end]
            if x_inr != 0
                px[x_inr] = y_end
                push!(tr_x, x_inr)
                (found, y_end) = malter!(adj_lst_zero, my, x_inr, px, py, tr_x, tr_y)
            else
                found = true
            end
        end
        j += 1
    end
    return found, y_end
end

# Hungarian algorithm
# adj_mtrx      input adjacency matrix of the input bipartite weighted graph
# mx            return matching adjacency vector from the X set
# my            return matching adjacency vector from the Y set
function hungalg(adj_mtrx)

    # adjacency vectors matching from set X and Y
    mx = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1))
    my = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1))

    # reduce the input matrix by transformation $T_{uv}$
    reduce!(adj_mtrx)

    # find the first 0 element in the first row, and set it in matching adj. vectors
    indx_zero = 1
    while adj_mtrx[1, indx_zero] != 0 && indx_zero <= size(adj_mtrx, 1)
        indx_zero += 1
    end
    mx[1] = indx_zero
    my[indx_zero] = 1

    # main loop of the algorithm   
    for x_sta in 2:size(adj_mtrx, 2)

        # adjacency list of vertices connected with zero weight edges
        adj_lst_zero = zeroadj(in_mtrx)

        # parents of the traversed vertices
        px = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1)) # from X
        py = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1)) # from Y

        # traversed vertices
        tr_x = [x_sta]                  # from X, starting vertex 
        tr_y = Vector{Int64}([])        # from Y

        # M-altering chain
        (found, y) = malter!(adj_lst_zero, my, x_sta, px, py, tr_x, tr_y)

        # if an M-altering chain is not found
        while found == false

            # transformation $O_{X_1 Y_1}$
            zeroxy!(in_mtrx, tr_x, tr_y)

            # adjacency list of vertices connected with zero weight edges
            # of the transformed matrix
            adj_lst_zero = zeroadj(in_mtrx)

            # complement of the traversed vertices from Y
            not_tr_y = setdiff(collect(1:size(adj_mtrx, 2)), tr_y)

            # find zero edge between traversed in X and not traversed from Y
            # that leads to an M-altering chain
            for i in tr_x, j in not_tr_y
                if in_mtrx[i, j] == 0                    
                    (found, y) = malter!(adj_lst_zero, my, i, px, py, tr_x, tr_y)
                    if found == true
                        break
                    end
                end
            end
        end

        # build the adjacency vectors matching from parent vectors
        x = py[y]
        mx[x] = y
        my[y] = x
        while x != x_sta
            y = px[x]
            x = py[y]
            mx[x] = y
            my[y] = x
        end        
    end

    return mx, my
end

# Calculate the price of the matching
# in_mtrx       input original adjacency matrix
# mtch          input matching adjacency vector
function mtchprice(in_mtrx, mtch)
    result = 0
    for i in mtch
        result += in_mtrx[i, mtch[i]]
    end

    return result
end

# Test transformation functions

# input matrices from text files
using DelimitedFiles

# input matrix to test transformations and M-altering chain
in_mtrx = readdlm("mtrx4-1.txt", Int64)
println("Transformations and M-altering chain test input matrix:")
display(in_mtrx)
print("\n")

# reduce the input matrix
reduce!(in_mtrx)
println("Reduced matrix:")
display(in_mtrx)
print("\n")

# transformation $O_{X_1 Y_1}$
x_set = [1, 2, 3]   # $X_1$ set
y_set = [1, 2]      # $Y_1$ set
zeroxy!(in_mtrx, x_set, y_set)
println("0X1Y1 transformation:")
display(in_mtrx)
print("\n")

# adjacency list of vertices connected with zero weight edges
adj_lst_zero = zeroadj(in_mtrx)
println("Zero adjacency list:")
display(adj_lst_zero)
print("\n")

# M-altering chain
my = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1)) # adjacency vector matching from set Y
px = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1)) # parents of the traversed vertices from X
py = zeros(typeof(in_mtrx[1, 1]), size(in_mtrx, 1)) # parents of the traversed vertices from Y
my[1] = 1                       # matching in the Y set
x_sta = 2                       # starting vertex in X
tr_x = [x_sta]                  # traversed vertices from X, starting vertex 
tr_y = Vector{Int64}([])        # traversed vertices from Y
(found, y_end) = malter!(adj_lst_zero, my, x_sta, px, py, tr_x, tr_y)
found ? print(y_end) : print("M-altering chain is not found.")
print("\nParent adjacency X: ")
display(px)
print("\nParent adjacency Y: ")
display(py)
print("\nTraversed list from X: ")
display(tr_x)
print("\nTraversed list from Y: ")
display(tr_y)
print("\n")

# Hungarian algorithm
in_mtrx = readdlm("mtrx50-1.txt", Int64)        # use it for the algorithm, and will modify
or_mtrx = copy(in_mtrx)                         # original input matrix
println("Hungarian algorithm input matrix:")
display(in_mtrx)
print("\n")
@time (mx, my) = hungalg(in_mtrx)
print("Matching vector in set X: ")
display(mx)
print("\nMatching vector in set Y: ")
display(my)
print("\n")
opt_prc = mtchprice(or_mtrx, mx)
print("Optimal price: ")
println(opt_prc)
