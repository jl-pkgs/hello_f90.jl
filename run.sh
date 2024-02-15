# rm build -rf
mkdir -p build

cd build && cmake .. && make
cd ..

julia jl_call_fortran.jl
