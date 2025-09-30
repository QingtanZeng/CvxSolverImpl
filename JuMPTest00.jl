using JuMP, ECOS

function solve_constrainedLSPgm(A::Matrix, b::Vector)
    m, n =size(A);
    model = Model(ECOS.Optimizer);
    @variable(model, x[1:n]);
    @variable(model, residuals[1:m]);
    @constraint(model, residuals == A*x-b);
    @constraint(model, sum(x)==1);
    @objective(model, Min, sum(residuals.^2));
    optimize!(model);
    return value.(x);
end

A, b = rand(10, 3), rand(10);

x = solve_constrainedLSPgm(A, b)