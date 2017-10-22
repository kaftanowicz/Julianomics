
using Calculus
using Plots

# Functions declaration

# Cobb-Douglas production function
# Default beta = 1-alpha gives constant return to scale (CRS)
CobbDouglas(K, AL, alpha, beta = 1-alpha) = K^alpha * AL^beta

# Production function to be used in model
Production(K, A, L) = CobbDouglas(K, A*L, 0.7)
    
# Capital in (t+1) as function of Capital in (t), Investment in (t)
# and depreciation rate delta
CapitalChange(K,Inv, delta) = (1-delta)*K + Inv

# Investment (equal to savings) as function of Production in (t) and savings rate s
Investment(Y, s) = s * Y

# Labour in (t+1) as function of Labour in (t) and growth rate of population n
LabourChange(L, n) = (1+n)*L

# Technology in (t+1) as function of Technology in (t) and growth rate of technology g
TechnologyChange(A, g) = (1+g)*A

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

mutable struct model_parameters
s::Float64
n::Float64
g::Float64
delta::Float64
end

parameters = model_parameters(0.04, 0.02, 0.01, 0.01)

function Solow(K0, A0, L0, parameters, T)
    # In t = 0:
    Y0 = Production(K0, A0, L0)
    Inv0 = Investment(Y0, parameters.s)
    C0 = Y0 - Inv0 # Goods market clearing condition
    W0, Rk0 = SolveFirmsProblem(Production, K0, A0, L0)
    
    # Initializing sequences of allocations and prices for times t=1..T
    for x = [:C, :Y, :K, :Inv, :Rk, :W, :A, :L]
        @eval $x = zeros(T)
    end
    
    t = 1
    K[t] = CapitalChange(K0, Inv0, parameters.delta)
    A[t] = TechnologyChange(A0, parameters.g)
    L[t] = LabourChange(L0, parameters.n)
    
    Y[t] = Production(K[t], A[t], L[t])
    Inv[t] = Investment(Y[t], parameters.s)
    C[t] = Y[t] - Inv[t] # Goods market clearing condition
    W[t], Rk[t] = SolveFirmsProblem(Production, K[t], A[t], L[t])
    
    for t = 2:T
        K[t] = CapitalChange(K[t-1], Inv[t-1], parameters.delta)
        A[t] = TechnologyChange(A[t-1], parameters.g)
        L[t] = LabourChange(L[t-1], parameters.n)
        
        Y[t] = Production(K[t], A[t], L[t])
        Inv[t] = Investment(Y[t], parameters.s)
        C[t] = Y[t] - Inv[t] # Goods market clearing condition
        W[t], Rk[t] = SolveFirmsProblem(Production, K[t], A[t], L[t])
    end
    
    C, Y, K, Inv, Rk, W, A, L
end

K0 = 100
A0 = 1.1
L0 = 300
T = 100

C, Y, K, Inv, Rk, W, A, L = Solow(K0, A0, L0, parameters, T)

plot(C, color="blue")
