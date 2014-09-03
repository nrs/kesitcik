kesitcik
==========
library for the numerical analysis of arbitrary composite beam cross sections
governed by the Euler-Bernoulli beam theory equations

The section is first discretized (using the award winning
program Triangle (c) by Shewchuk) and then analyzed according to
linear flexural strain relations. The main purpose of this project is
to provide better estimations for the strength of composite cross
sections, such as reinforced concrete cross sections.

Installation
------------
Run ``make`` in the main directory.

Alternatively on a Windows system, you can import to Visual Studio or
another IDE, or install Cygwin with necessary building tools.

Notes
-----
* Triangle (c) by Jonathan Shewchuk is used for 2D triangulation.
* Some classes are from the open source finite element library libmesh.
* Feel free to contribute, send pull requests.
