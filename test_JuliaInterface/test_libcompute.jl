 const LIB_PATH = "./test_JuliaInterface/libcompute.so"

 a = 3.0
 b = 4.0

 result =ccall((:compute_hypot, LIB_PATH), Cdouble, (Cdouble, Cdouble), a, b)
 println("hypot($a, $b) = $result")

data = [1.0, 2.0, 3.0, 4.0, 5.0]
factor = 10.0
size_int = Int32(length(data))

println("缩放前: ", data)
ccall((:scale_array, LIB_PATH), Cvoid, (Ptr{Cdouble}, Cint, Cdouble), data, size_int, factor)
println("缩放后: ", data)