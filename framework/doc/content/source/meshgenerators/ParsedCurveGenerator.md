# ParsedCurveGenerator

!syntax description /Mesh/ParsedCurveGenerator

## Overview

The `ParsedCurveGenerator` object generates a 3D curve mesh composed of EDGE2 elements which connect the series of points given by $[x(t),~y(t),~z(t)]$, where the range of t is specified by the user.

!equation id=xy_formula
\begin{cases}
  x = x(t)\\
  y = y(t)\\
  z = z(t)
\end{cases},

where, $x(t)$, $y(t)$, and $z(t)$ are all provided in the form of [C++ function parser](http://warp.povusers.org/FunctionParser/) through [!param](/Mesh/ParsedCurveGenerator/x_formulus), [!param](/Mesh/ParsedCurveGenerator/y_formulus), and [!param](/Mesh/ParsedCurveGenerator/z_formulus), respectively. The constants used in these formula can be defined by [!param](/Mesh/ParsedCurveGenerator/constant_names) and [!param](/Mesh/ParsedCurveGenerator/constant_expressions).

!media reactor/meshgenerators/xyz_curve.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=xyz_curve
      caption=A 3D open curve generated by `ParsedCurveGenerator` defined as $x(t)\=\cos t$; $y(t)\=\sin t$; $z(t)\= t$ with $t$ ranging from $0$ to $4\pi$.

Key $t$ values including the starting and ending values of $t$ must be specified by [!param](/Mesh/ParsedCurveGenerator/critical_t_series). Optionally, intermediate $t$ values can be added so that the curve can be divided into several sections. Note that the elements in [!param](/Mesh/ParsedCurveGenerator/critical_t_series) must be unique and change monotonically. Each interval defined by [!param](/Mesh/ParsedCurveGenerator/critical_t_series) must have a corresponding number of segments (i.e., EDGE2 elements), $N_{seg}$, defined by [!param](/Mesh/ParsedCurveGenerator/nums_segments).

## Calculating Segment Division Points

Each section of the curve should ideally have segments of similar length. However, it is challenging to predict the corresponding $t$ values that yield segments with similar length. Hence, oversampling is used to determine the $t$ values which result in consistent segement length. The oversampling factor $N_{os}$ can be defined through [!param](/Mesh/ParsedCurveGenerator/oversample_factor). Assuming that a section of curve is defined by $t_n$ and $t_{n+1}$, the distance between the starting and ending points of this section has the following form,

!equation id=section_distance
d_n = \sqrt{\left[x(t_n)-x(t_{n+1})\right]^2+\left[y(t_n)-y(t_{n+1})\right]^2}

Thus, the oversampling target is to achieve that the maximum interval between neighboring sampling points is lower than a threshold value defined as follows,

!equation id=threshold_distance
d_{os,threshold} = \frac{d_n}{N_{seg}N_{os}}

The oversampling is realized by a binary algorithm, which divides oversized intervals in halves until all the intervals are shorter than $d_{os,threshold}$. Then the oversampled section points are used to determine the actual point locations (i.e., $t$ values).

## Example Syntax

`ParsedCurveGenerator` is capable of generating both open and closed curves.

For open curve generation, the approach is straight forward with the example shown as follows.

!listing test/tests/meshgenerators/parsed_curve_generator/parsed_curve_open.i block=Mesh/pcg
         id=open_curve_input
         caption=the input syntax sample to generate an open logarithmic curve

!media reactor/meshgenerators/open_curve.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=open_curve
      caption=An open logarithmic curve generated by `ParsedCurveGenerator`.

On the other hand, for closed curve generation (defined by [!param](/Mesh/ParsedCurveGenerator/is_closed_loop)), ideally the starting and ending values of $t$ should lead to the same $x(t)$ and $y(t)$ values, as shown below.

!listing test/tests/meshgenerators/parsed_curve_generator/parsed_curve_closed.i block=Mesh/pcg
         id=closed_curve_input
         caption=the input syntax sample to generate a closed half circular and half elliptical curve

!media reactor/meshgenerators/closed_curve.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=closed_curve
      caption=A closed half circular and half elliptical curve generated by `ParsedCurveGenerator`.

If the starting and ending values of $t$ lead to different $x(t)$ and $y(t)$ values, the curve will be "forced" to close by directly connection the starting and ending points (see [forced_closed_curve]).

!media reactor/meshgenerators/forced_closed_curve.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=forced_closed_curve
      caption=a fraction of the curve shown in [closed_curve] forced to be closed.

## Used with Other Mesh Generators

If [!param](/Mesh/ParsedCurveGenerator/z_formulus) is set as zero (default value), the generated curve resides in the XY-plane, and a pair of such `ParsedCurveGenerator` objects can naturally be connected by [`FillBetweenCurvesGenerator`](/FillBetweenCurvesGenerator.md) using [`FillBetweenPointVectorsTools`](/FillBetweenPointVectorsTools.md).

Additionally, closed XY-plane curve meshes generated by `ParsedCurveGenerator` can be used by [`XYDelaunayGenerator`](/XYDelaunayGenerator.md) as either [!param](/Mesh/XYDelaunayGenerator/boundary) or [!param](/Mesh/XYDelaunayGenerator/holes). See example below.

!listing test/tests/meshgenerators/parsed_curve_generator/xy_curve.i block=Mesh
         id=xy_curve_list
         caption=the input syntax sample to use `ParsedCurveGenerator` with [`XYDelaunayGenerator`](/XYDelaunayGenerator.md)

!media reactor/meshgenerators/xy_curve.png
      style=display: block;margin-left:auto;margin-right:auto;width:40%;
      id=xy_curve_fig
      caption=An example mesh generated by using `ParsedCurveGenerator` with [`XYDelaunayGenerator`](/XYDelaunayGenerator.md)

!syntax parameters /Mesh/ParsedCurveGenerator

!syntax inputs /Mesh/ParsedCurveGenerator

!syntax children /Mesh/ParsedCurveGenerator