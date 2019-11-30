using GpSelection

my_tests = ["test_GpSelection.jl"]

println("Running tests:")

for my_test in my_tests
    println(" * $(my_test)")
    include(my_test)
end
