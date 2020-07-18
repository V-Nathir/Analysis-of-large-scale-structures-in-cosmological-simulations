# Tools-for-data-from-a-cosmological-simulation
My final degree project in physics: Inertia tensor of structures at different scales in cosmological simulations.<br/>
A compilation of diferent tools and packages to calculate some physics parameters of interest and tools like: reading Big Endian files, creating Aitoff's projections and a lot of plots (scatter3d, scatter2d, plot3d,barplots with modifications...), calculating eigenvalues and eigenvector, reduced inertia tensor, inertia ellipsoids, normalized angular momentum and morphology of large structures, creating and saving files, creating and reading directories, techniques to calculate center of mass with several tolerances...
<br/>
In these 5 volumes you can find different types of mechanics to plot, to implement tolerance in loops and conditional structures ... in Python.<br/>
You can use these mechanics to inspire your code.
<br/>
<br/>
***Warning: The Data cannot be upload for copyright. So you can not running the code.***
<br/>
<br/>
## Final Degree Project
<br/>
-- _Abstrac_ --
<br/>
Data from an N-Body cosmological simulation written in Big Endian code will be analysed. Such data will be interpreted and packed in useful packages to study particle structures on different scales. Particularly, the values used will be of mass, position, velocity, and type of baryonic particle (stars and gas). It will be discerned if the relative positions and gravitational interactions between structures could create new persistent or own orbital plans, or important instabilities that do not allow these and aggressively integrating in the cosmic web (CW). 
<br/>
The inertia tensor (reduced) of 4 structures (group of particles-> stars+gas-> galaxies) will be calculated: web 2, identified as cosmic web, web 3 attached to web 2, web of persistent satellites (SP), with hypothesis of a future orbital persistent plan, and web of non-persistent satellites (SNP), with opposing hypothesis. This will allow determining the main semi-axes and eigenvectors associated to each eigenvalues in the principal directions.  And this will reduce a structure into an inertia ellipsoid with which changes and evolutions can be dated. 
<br/>
A possible correlation confirming the previous hypothesis has been discovered. For SP: the entry angle of the smaller semi-axis is tangential with the one of web 2, a stable orbital plan and with its own morphology will exist. For SNP: the entry is done in parallel, the creation of an stable orbital plan is not ensured, and it will be swiftly accreted to the cosmic web, with a different repercussion in the changes of the mass centres, angular momentum (modified), morphology, and the star-gas relation.
<br/>
<br/>

## Code Info:
<br/>

  .-Vol I and VolI_launch : VolI_launch import the code of Vol I. Convert Big Endian code to numpy files and save them. 
<br/>

  .-Vol II: create and save conversions tables from Simulation units to Internacional System.
 <br/>

  .-Vol III: load files and print it in the terminal. Check-out of if all files are created. 
 <br/>

  .-Particle Vol IV: calculate the center of mass (CM) applying loops structures and conditionals that can be adapted to each situation. The CM must be between a lower and upper bound generate for each situation automatically. Create a identification list to vinculate the real number of a particle and the modificated number after load it. Plot a massive number of particles using scatter. With these plots it can be check-out how the techniques to calculate CM are working. Count the number of type of a particle: stars and gas, also dark matter. Finally, it can find particles around a point and track its individual motion.
 <br/>

  .-Tensor Vol V: calculate the inertia tensor, angular momentum and morphology of the large structures. Calculate the eigenvectors and eigenvalues. Check-out the relation gas/star. Plots: Aitoffs proyections, scatter, plots, ... 
