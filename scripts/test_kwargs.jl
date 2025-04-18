#!/usr/bin/env julia

# Test script for keyword argument handling

println("Testing Keyword Argument Handling")
println("=================================")

# Define a function with keyword arguments
function test_function(a, b, c; option1=false, option2="default")
    println("Regular arguments: a=$a, b=$b, c=$c")
    println("Keyword arguments: option1=$option1, option2=$option2")
    return a + b + c
end

# Define a function that calls another function with keyword arguments
function call_with_semicolon(a, b, c; option1=false, option2="default")
    # This passes keyword args with semicolon syntax
    result = test_function(a, b, c; option1=option1, option2=option2)
    println("Called with semicolon - Result: $result")
    return result
end

function call_without_semicolon(a, b, c; option1=false, option2="default")
    # This passes keyword args without semicolon syntax
    result = test_function(a, b, c, option1=option1, option2=option2)
    println("Called without semicolon - Result: $result")
    return result
end

# Run the different versions
println("\nTest 1: Direct call with semicolon syntax")
test_function(1, 2, 3; option1=true, option2="custom")

println("\nTest 2: Direct call without semicolon syntax")
test_function(1, 2, 3, option1=true, option2="custom")

println("\nTest 3: Indirect call with semicolon syntax")
call_with_semicolon(1, 2, 3; option1=true, option2="custom")

println("\nTest 4: Indirect call without semicolon syntax")
call_without_semicolon(1, 2, 3; option1=true, option2="custom")