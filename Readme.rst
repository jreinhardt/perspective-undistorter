Perspective undistorter
#######################

This script allows you to create a persectively undistorted image, that, when
viewed from a certain angle, looks like it is not perspectively distorted, even
though it should. This is also called anamorphosis and is for example useful,
when you want to paint something on the road, such that passing drivers can
read it without effort.

Published under MIT Expat Licence


Prepare input files
===================

The perspective undistorter allows to calculate a undistorted image from two
different types of inputs. Either a two dimensional drawing from a SVG file is
used, or three dimensional objects from a STL file is used. Some points need to
be considered, when preparing input files

1. You need to paint all lines by hand in the end, so keep the complexity low.
2. The SVG input capabilities are very limited, essentially only paths with
straight line segments are supported. So you need to convert everything to a
path, and convert all the segments to straight line segments. Choose a coarse
approximation of round curves, otherwise there will be very many lines to draw.
An example for that is hello.svg in the example folder.
3. The STL input capabilities are also very limited, all triangles are broken
up into lines, so you will potentially end up with very many lines. Try to
restrict yourself to very simple, straight geometrical shapes (cubes,
pyramids).
3. The perspective undistorter only takes care to correct the geometry,
coloring, lighting and shading needs to be done by you.

Prepare filter files
====================

STL files are based on a triangular boundary representation. As a result of
this, if you try to work with shapes that are not made of triangles, you will
end up with extra lines. One can suppress these additional lines using the
--filter argument or the [filtered] section of the filter file specified with
the --filter-file command line argument. To find the indices of the lines you
want to suppress, run the perspectice undistorter on a file, look at the
resulting svg and look up the indices there. The --colorful option is very
useful for this. With STL files your initial result will be very crowded, so
usually this is a iterative process. For 3D scenes you might want to suppress
also lines that are hidden by other objects to get a better illusion. If lines
are very close together, it is also sometimes useful to consult the point and
lines list.

One can also mark partially hidden lines by specifying their indices in the
[partial] section of the filter file. These lines show up in a different color
in the overview drawing and in a special section of the line list.

Finally there is also the possibility to add extra lines in the [extra]
section, which sometimes allows to simplify a drawing by replacing several
lines by a single one.

Output
======

Two files are written by the perspective undistorter:

* filename_perspective.svg: A overview file with numbered points and lines,
including the origin and reference point. This is useful for orientation when
painting and to find edges that you want to filter during filter file creation.
* filename_perspective.table: A file containing an overview over all points and
lines. Use this when painting. The distances to origin and reference point are
specified for easy triangulation.

Drawing the thing
=================
The process is much easier if two or more people can work together. 

To be able to accurately determine the position of points for a painting that
might extend over several tens of meters, triangulation is used. In the tables
file all points are specified by their distance to the origin (the point where
you have to stand to enjoy the illusion), and a second reference point. Begin
by choosing these two points with the correct distance and direction. To get
the orientation right, check the overview SVG. Then fix two sufficiently long
tape measures (flexible ones are better than steel tape) at these points.  Now
all other points can be found by triangulation with the distances given in the
table file. The points can be marked by small flags, pieces of tape, sharpie
markings, chalk or other markers, depending on the scale of the drawing and the
character of the ground.

After you got all points right, add the lines. Be careful with partially hidden
lines, draw them at the end. Lines can be marked with string, chalk, spray
chalk, tape or other means. Comments:

* Solid sidewalk chalk: Easy to obtain in toy stores, but usually a limited
range of relatively light colors, does not cover too well.
* Spray chalk: Covers very well, but masking is necessary to get nice lines and
sharply bounded areas. Difficult to use when windy. One can lighten and darken
areas by gently spraying white or black chalk from a larger distance. Is pretty
resilent to rain, usually stays on for a few weeks. Larger paintings require
several cans.

If you succeeded, please take a picture and post it to the wiki. I would love to hear about it.

Known Problems
==============

Very limited svg support
