# OrbitalMechanics.jl
This Project is my attempt to use machine learning to optimize ascent trajectories in kerbal space program
Essentially it is supposed to get into orbit using less delta-V then I could have done by hand. 

The model being trained dictates the angle of the (negative) thrust vector vector from the vertical at a given point in time.
It does this by storing a vector of angles and a rate, which it uses a cubic spline to turn into a function from time to angle.

The model is evaluated by running a simulation of a ship going into orbit in three stages, then calculating the deltaV of the ship once it's in orbit

1. Ascent: the ship burns until the apoapsis is above some target altitude (plus some margin). 
2. Coast: the ship cuts the engines and waits until it is at or above the target altitude.
3. Circularize: the ship burns at 90\degree from the virtical until the periapsis is at or above the target altitude

With model controlls the direction of the thrust during the ascent stage. 
