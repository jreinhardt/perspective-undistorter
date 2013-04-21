Perspective Undistorter
#######################

This script allows you to create a persectively undistorted image, that, when
viewed from a certain angle, looks like it is not perspectively distorted, even
though it should. This is for example useful, when you want to paint something
on the road, such that passing drivers can read it without effort.

Published under MIT Expat Licence


Prepare input files
===================

The perspective undistorter allows to calculate a undistorted image from two different types of inputs. Either a two dimensional drawing from a SVG file is used, or three dimensional objects from a STL file is used. Some points need to be considered, when preparing input files

1. You need to draw all lines by hand in the end, so keep the complexity low.
2. The SVG input capabilities are very limited, essentially only paths with straight line segments are supported. So you need to convert everything to a path, and convert all the segments to straight line segments.
3. The STL input capabilities are also very limited, all triangles are broken up into lines, so you will potentially end up with very many lines. Try to restrict yourself to very simple, straight geometrical shapes (cubes, pyramids).
3. The perspective undistorter only takes care to correct the geometry, coloring, lighting and shading needs to be done by you.

Prepare filter files
====================

STL files are based on a triangular boundary representation. As a result of this, if you try to work with shapes that are not made of triangles, you will end up with extra lines. One can suppress the output of these helper lines using the --filter and --filter-file command line arguments. To find the indices of the lines you want to suppress, run the perspectice undistorter on a file, look at the resulting svg and look up the indices there. With STL files your initial result will be very crowded, so usually this is a iterative process.


Known Problems
==============

Very limited svg support
