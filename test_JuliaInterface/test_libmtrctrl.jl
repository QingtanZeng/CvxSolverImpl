const LIB_PATH = "./test_JuliaInterface/libmtrctrl.so"

mutable struct MtrCtrlWrapper
    ptr::Ptr{Cvoid}

    function MtrCtrlWrapper()
        ptr = ccall((:MtrCtrl_new, LIB_PATH), Ptr{Cvoid}, ())

        obj = new(ptr)

        finalizer(destroy!, obj)

        return obj
    end
end

function destroy!(m::MtrCtrlWrapper)
    if m.ptr == C_NULL
        ccall((:MtrCtrl_delete, LIB_PATH), Cvoid, (Ptr{Cvoid},), m.ptr)
        m.ptr = C_NULL
    end
end

function set_speed!(m::MtrCtrlWrapper, speed::Float64)
    # 将 Julia 对象的 ptr 传给 C++ 的 obj 参数
    ccall((:MtrCtrl_set_speed, LIB_PATH), Cvoid, (Ptr{Cvoid}, Cdouble), m.ptr, speed)
end

function get_speed(m::MtrCtrlWrapper)::Float64
    return ccall((:MtrCtrl_get_speed, LIB_PATH), Cdouble, (Ptr{Cvoid},), m.ptr)
end

mtrctrl = MtrCtrlWrapper()

set_speed!(mtrctrl, 5000.3)
current_speed = get_speed(mtrctrl)
println("[Julia] 读取到的电机速度: ", current_speed)

# 强制触发 Julia 的垃圾回收 (GC)，验证 C++ 的 Destructor 是否被正确调用
println("\n[Julia] 准备销毁对象，触发 GC...")
mtrctrl = nothing 
GC.gc() # 强制回收内存

set_speed!(mtrctrl, 5000.3)