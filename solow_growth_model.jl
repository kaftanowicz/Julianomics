
using Calculus
using Plots

# Functions declaration

function CobbDouglas(K, AL, alpha, beta)
    # Cobb-Douglas production function
    K^alpha * AL^beta
end

function CobbDouglasCRS(K, AL, alpha)
    # Cobb-Douglas production function with constant returns to scale
    CobbDouglas(K, AL, alpha, 1-alpha)
end

function Production(K, A, L)
    # Production function
    CobbDouglasCRS(K, A*L, 0.7)
end

function CapitalChange(K, I, delta)
    # Capital in (t+1) as function of Capital in (t), Investment in (t)
    # and depreciation rate delta
    (1-delta)*K + I
end

function Investment(Y, s)
    # Investment (equal to savings) as function of Production in (t) and savings rate s
    s * Y
end

function LabourChange(L, n)
    # Labour in (t+1) as function of Labour in (t) and growth rate of population n
    (1+n)*L
end

function TechnologyChange(A, g)
    # Technology in (t+1) as function of Technology in (t) and growth rate of technology g
    (1+g)*A
end

function ExtractParameters(parameters_dict)
    # Takes a Dictionary and creates variables with names corresponding to its key names
    # and values corresponding to its values
           expr = quote end
    for (k, v) in parameters_dict
               push!(expr.args, :($(Symbol(k)) = $v))
           end
           eval(expr)
           return
       end

function SolveFirmsProblem(ProductionFunction, K, A, L)
    production_function_gradient = Calculus.gradient(x -> ProductionFunction(x[1], x[2], x[3]), [K, A, L])
    W = production_function_gradient[3]  # W = dF/dL
    Rk = production_function_gradient[1] # Rk = dF/dK
    W, Rk
end

function Solow(K0, A0, L0, parameters, T)
    ExtractParameters(parameters)
    
    # In t = 0:
    Y0 = Production(K0, A0*L0)
    I0 = Investment(Y0, s)
    C0 = Y0 - I0 # Goods market clearing condition
    W0, Rk0 = SolveFirmsProblem(Production, K0, A0, L0)
    
    # Initializing sequences of allocations and prices for times t=1..T
    for x = [:C, :Y, :K, :I, :Rk, :W, :A, :L]
        @eval $x = zeros(T)
    end
    
    t = 1
    K[t] = CapitalChange(K0, I0, delta)
    A[t] = TechnologyChange(A0, g)
    L[t] = LabourChange(L0, n)
    
    Y[t] = Production(K[t], A[t]*L[t])
    I[t] = Investment(Y[t], s)
    C[t] = Y[t] - I[t] # Goods market clearing condition
    W[t], Rk[t] = SolveFirmsProblem(Production, K[t], A[t], L[t])
    
    for t = 2:T
        K[t] = CapitalChange(K[t-1], I[t-1], delta)
        A[t] = TechnologyChange(A[t-1], g)
        L[t] = LabourChange(L[t-1], n)
        
        Y[t] = Production(K[t], A[t]*L[t])
        I[t] = Investment(Y[t], s)
        C[t] = Y[t] - I[t] # Goods market clearing condition
        W[t], Rk[t] = SolveFirmsProblem(Production, K[t], A[t], L[t])
    end
    
    C, Y, K, I, Rk, W, A, L
end

K0 = 100
A0 = 1.1
L0 = 300
parameters = Dict("s" => 0.04, "n" => 0.02, "g" => 0.01, "delta" => 0.01)
T = 100

C, Y, K, I, Rk, W, A, L = Solow(K0, A0, L0, parameters, T)

plot(C, color="blue")
