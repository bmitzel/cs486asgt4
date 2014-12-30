 * Programmer: Brian Mitzel
 * Email: bmitzel@csu.fullerton.edu
 * Course: CPSC 486
 * Assignment: 4 - Separating Axes
 * Due Date: 12/17/2014
 
This program demonstrates the separating axes algorithm
for collision detection between 2D polygons. A number of
squares are rendered on the screen with variable starting
locations and velocities. When a square collides with
another square, its velocity vector is reversed. Also,
when a square collides with an edge of the window, it is
reflected.

As the squares move around the screen, they will slowly
fade out to zero intensity (black). However, when one
square collides with another square, it will change color
and appear at full intensity again.

------------------------------------------------------------
    
Features:

The following hotkeys are available:
    1 - slow down updates to once every 500 milliseconds
    2 - speed up updates to once every 16 milliseconds
    d - toggle the printing of the rectangle data
    r - reset all of the rectangles
    ESC or q - quit the program
    
------------------------------------------------------------

Completed:

The following requirements for this assignment have been
completed:
    1. The program contains 22 rectangles that are colored
       using 2 different hues, green and blue, at varying
       levels of intensity.
    2. The motion of the rectangles is linear with no
       acceleration and constant speed.
    3. Collision detection is calculated using the
       separating axes technique.
    4. A collision with another rectangle results in the
       object reversing its velocity.
    5. A collision with a window boundary results in the
       object being reflected.

------------------------------------------------------------

Bonus:

The color intensity of each polygon slowly degrades over
time, causing the polygons to fade out to black.

However, whenever one polygon collides with another, it both
changes color (from blue to green or vice versa) and returns
to full intensity.
       
------------------------------------------------------------
       
Bugs:

There are currently no known bugs or limitations.

------------------------------------------------------------

Building:

Using OSX, Linux, or the FreeBSD VM distributed by Professor
Michael Shafae, open a terminal window to the source
directory and enter the following command at the shell
prompt:

    make

or:

    gmake

------------------------------------------------------------
    
Executing:

From the same directory where you built the executable,
enter the command:

    ./moving_boxes
