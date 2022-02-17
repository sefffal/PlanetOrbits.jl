
# Image Warping
If you have an image of a system, you can warp the image as if each pixel were a test particle following Kepler's laws. 
This is an easy way to see what a disk or a system of planets would look like at a time other than when it was captured.

To make this possible, DirectOrbits.jl can create `OrbitalTransformation` objects. These follow the conventions set out
in CoordinateTransformations.jl and are compatible with ImageTransformations.jl.

Example:
```julia
ot = OrbitalTransformation(
    i = 0.3,
    e = 0.1,
    M = 1.0,
    ω = 0.5,
    Ω = 0.5,
    plx = 30.0,
    
    platescale=10.0, # mas/px
    dt = 3*365.25 # days forward in time
)

img_centered = centered(img)
img_future = warp(img_centered, ot, axes(i))

# Display with DirectImages.jl
using DirectImages
imshow2([img; img_future], clims=(0,1), cmap=:seaborn_icefire_gradient)
```

**Before, and After Orbital Transformation**

![image](https://user-images.githubusercontent.com/7330605/129625363-c0295432-47f4-4400-a5a7-7140a7e7d997.png)

Note the arguments `platescale` and `dt` are required, but `a` and `τ` are not. The position of the pixel in X/Y space uniquely determines the semi-major axis and epoch of periastron passage when the rest of the orbital parameters are known. `platescale` in units of milliarseconds/pixel is necessary to get the overall scale of the transform correct. This is because an orbital transformation is **not** linear (and therefore, care must be taken when composing an OrbitalTransformation with other CoordinateTransformations). Scaling an image will change the amount of rotation that occurs at each separation. `dt` is the the amount of time in days to project the image forward. It can also be negative to project the image into the past. 
