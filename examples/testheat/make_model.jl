
function make_model_H20(N = 25)
    fens, fes = H20block(A, A, A, N, N, N)
    ir = GaussRule(3, 3)
    return fens, fes, ir
end

function make_model_T10(N = 25)
    fens, fes = T10block(A, A, A, N, N, N)
    ir = TetRule(5)
    return fens, fes, ir
end

make_model = make_model_H20
make_model = make_model_T10
