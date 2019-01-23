using Images, Makie, ShapeFromShading

#generate synthetic image
img = generate_surface(0.5, [0.2,0,0.9], radius = 5)

#calculate the heightmap
Z = retrieve_surface(Pentland(), img)
@benchmark retrieve_surface(Pentland(), img)

#normalize to maximum of 1 (not necessary but makes displaying easier)
Z = Z./maximum(Z)

#display using Makie
r = 0.0:0.1:2
surface(r, r, Z)
