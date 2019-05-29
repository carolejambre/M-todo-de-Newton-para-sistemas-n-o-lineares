using LinearAlgebra
"""`NewtonNaolinear(F,J,x; atol = 1e-6, rtol = 1e-6, fdertol = 1e-12, Max_iter = 1e6)`

Calcula o vetor solução do sistema não-linear `F` utilizando o método de Newton.


Entradas da Função:

`F` - Sistema Não-Linear - `F(x) = 0`;

`J` - Matriz Jacobiana da função `F` - Matriz que tem como entradas os vetores gradientes das `fi(x)'s` de `F`;
(Se `Det(J) = 0`, então retorna um erro)

`x` - vetor aproximação inicial para o sistema.

Saídas da Função:

`xₖ`    - Aproximação para um zero de `F`;

`F(xₖ)` - F aplicada nessa aproximação;

`k`     - Número de interações.

"""
function NewtonNaoLinear(F,J,x; atol = 1e-6, rtol = 1e-6, fdertol = 1e-12, Max_iter = 1e6)
    k = 0
    fx = F(x)
    ϵ = atol + rtol * norm(fx)
    jx = J(x)

    if det(jx) == 0
        error("O determinante da matriz Jacobiana não pode ser 0")
    end
    while k ≤ Max_iter && norm(fx) > ϵ
        k = k + 1
        s = jx\-fx
        x = x + s
        fx ,jx = F(x), J(x)
    end
println("Aproximação para um zero de F: x = $x")
println("F(x) = $(F(x))")
println("Número de interações necessárias: $k")

end
