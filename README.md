# OpenGL_RayTracing
Project of Selective Course 'Computer Graphics'. Implement ray tracing which is a useful technique to generate realistic image of 3D scene.  

***Main Implementation of Interface:***  
There are five main functions for the pipeline except for the template: 
- `1.TraceRay(in main.cpp):`  
Given a ray, calculate the color at the nearest intersection point.  

- `2.Shade(in main.cpp):`  
Given the ray-surface intersection information, calculate the corresponding color.  

- `3.Triangle::Hit(in hittable.cpp):`  
Perform intersection calculation between a ray and a triangle.  

- `4.Sphere::Hit(in hittable.cpp):`  
Perform intersection calculation between a ray and a sphere.  

- `5.Quadric::Hit(in hittable.cpp):`  
Perform intersection calculation between a ray and a quadric surface.  

  
***Input format:***  
<p align="center">
<img src="https://github.com/Yu-1218/OpenGL_RayTracing/blob/main/Input_format.png" width="250">
</p>

***Output Rendering result:***  
<p align="center">
<img src="https://github.com/Yu-1218/OpenGL_RayTracing/blob/main/Output_format.png" width="500">
</p>
