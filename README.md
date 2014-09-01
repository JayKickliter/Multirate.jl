[![Build Status](https://travis-ci.org/JayKickliter/Multirate.jl.svg?branch=master)](https://travis-ci.org/JayKickliter/Multirate.jl)
[![CoverageStatus](https://coveralls.io/repos/JayKickliter/Multirate.jl/badge.png)](https://coveralls.io/r/JayKickliter/Multirate.jl)

**Multirate** is a package for the creation and execution of state-preserving FIR filters which perform sample rate changes.


# Installation

```julia
Pkg.add( "Multirate" )
```

**Multirate** depends on **[DSP](https://github.com/JuliaDSP/DSP.jl)** for [windowing](http://en.wikipedia.org/wiki/Window_function) functions.

# Usage

## Stateless ##

Use these methods if you only want to resample a vector in one chunk.

In these examples, `x` and `h` are previously defined signal & filter-taps vectors.

```julia
using Multirate

# decimate by 3
y = filt( h, x, 1//3 )

# interpolate by 3 
y = filt( h, x, 3//3 )

# resample with a ratio of 3/4
# equivalent to interpolating by 3, then decimating by 4
# ( but much more efficient )
y = filt( h, x, 3//4 )
```

## State-preserving ##

To use **Multirate**'s stateful filtering functionality, you must first create a filter object. Every time you call `filt` with that object, the filtering picks up where it left off. This is good for processing large vectors from a file, or filtering an stream of samples with indefinite length.

Each filter object is of type `FIRFilter{Tk<:FIRKernel}`. Thera are four subtypes of `FIRKernel`:

* `FIRStandard`: normal single-rate FIR filter
* `FIRInterpolator`: resampling with a ratio of `L//1`
* `FIRDecimator`: resampling with a ratio of `1//M`
* `FIRRational`: resampling with a ratio of `L/M`

You do not need to specify the kernel type. It is chosen for you when you based on the resampling ratio you specify when creating a new `FIRFilter` object.

In the following example, we will be resampling with a ratio of `3//17`. Please note that the filter taps `h` in this example is contrived to show you the input to output progressions. It performs no useful signal filtering.

```jlcon
julia> x = [ 1.0:100 ]
100-element Array{Float64,1}:
   1.0
   2.0
   3.0
   4.0
   5.0
   6.0
   ⋮
  95.0
  96.0
  97.0
  98.0
  99.0
 100.0


julia> h = [ ones(3), zeros(6) ]
9-element Array{Float64,1}:
 1.0
 1.0
 1.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0

julia> myfilt = FIRFilter( h, 3//17 )
FIRFilter{FIRRational}(FIRRational(3x3 Array{Float64,2}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 1.0  1.0  1.0,3//17,3,3,0,1,1),[0.0,0.0],2)

julia> y1 = filt( myfilt, x[1:5] )
1-element Array{Float64,1}:
 1.0

julia> y2 = filt( myfilt, x[6:23] )
4-element Array{Float64,1}:
  6.0
 12.0
 18.0
 23.0

julia> y3 = filt( myfilt, x[24:100] )
13-element Array{Float64,1}:
 29.0
 35.0
 40.0
 46.0
 52.0
 57.0
 63.0
 69.0
 74.0
 80.0
 86.0
 91.0
 97.0

julia> y = [ y1, y2, y3 ]
18-element Array{Float64,1}:
  1.0
  6.0
 12.0
 18.0
 23.0
 29.0
  ⋮
 69.0
 74.0
 80.0
 86.0
 91.0
 97.0
```

Let's check that `y`, created by filtering three separate chunks of `x`, matches the result we would obtain from stateless filtering.

```jlcon
julia> sum( y .- filt( h, x, 3//17 ) )
0.0
```