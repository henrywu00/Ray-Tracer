
# Ray Tracer

![alt text](https://github.com/henrywu00/Ray-Tracer/blob/main/test-.png)

Image output is test-.png

Used OpenMP for multiple CPU for calculation

Made an animation by revolving a light source around the scene. This is done in the main function. Please check the comment in main(). output file is output.mp4

Input file is a text document called cornell.txt.

# Input format
```
Background color
Depth

number of lights

FOR EACH_LIGHT
LIGHT_POS
LIGHT_DIFFUSE
LIGHT_SPECULAR

number of spheres

FOR EACH_SPHERE
SPHERE_POS
SPHERE_RADIUS
SPHERE_DIFFUSE
SPHERE_SPECULAR
SPHERE_SHININESS
SPHERE_TRANSLUCENT SPHERE_INDEX

number of quads

FOR EACH_QUAD
POINT_A
POINT_B
POINT_C
QUAD_DIFFUSE
QUAD_SPECULAR
QUAD_SHININESS

resX resY
```
