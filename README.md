# AscentOptimization.jl

## Todo

 - [X] write a basic physics simulation 
 - [X] build a model that dictates how the ship will ascend
 	- [ ] improve the circularizing maneuver
 - [X] train the model using reinforcement learning
 - [X] train the model using gradient descent
 - [X] make simulation and training handle unitful quantities
 - [ ] model the atmosphere in the simulation
 - [ ] use actual parameters from Kerbal Space program
 	- [ ] validate on an atmosphere-less body like the mun
 - [ ] write the kos code to control the ship during ascent
 

## Methodology (WIP)
This Project is my attempt to use machine learning to optimize ascent trajectories in kerbal space program
Essentially it is supposed to get into orbit using less delta-V then I could have done by hand. 

The model being trained dictates the angle of the (negative) thrust vector vector from the vertical at a given point in time.
It does this by storing a vector of angles and a rate, which it uses a cubic spline to turn into a function from time to angle.

The model is evaluated by running a simulation of a ship going into orbit in three stages, then calculating the Delta V of the ship once it's in orbit

1. Ascent: the ship burns until the apoapsis is above some target altitude (plus some margin). 
2. Coast: the ship cuts the engines and waits until it is at or above the target altitude.
3. Circularize: the ship burns at 90\degree from the vertical until the periapsis is at or above the target altitude

With model controls the direction of the thrust during the ascent stage. 

## How this project is organized.

This project is roughly organized into two parts

1. A library of miscellaneous code that includes all the physics and machine learning stuff I need

2. A set of notebooks that I use to visually inspect and prototype the code
