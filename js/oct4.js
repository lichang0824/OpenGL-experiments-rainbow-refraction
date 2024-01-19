
rooms.oct4 = function() {

description = `<b>Notes for October 4</b>
               <p>
	       Refraction +
	       <br>
	       Ray tracing to 2nd order surfaces`;

code = {

'refraction':`
S.html(\`
<b>Refraction in more detail:</b>
<blockquote>

Given an incoming ray in direction W<sub>1</sub>, a surface normal N,
incoming index of refraction <i>n</sub>1</sub></i> and
outgoing index of refraction <i>n</sub>2</sub></i>,
we can compute the direction of the outgoing ray
W<sub>2</sub> in seven steps, as shown in the diagram:

<ol>
<li>
We note that the component of incoming ray W<sub>1</sub>
that is parallel to normal N is given by:
<font color=blue>
C<sub>1</sub> = N (W<sub>1</sub> &bull; N)
</font>
<p>

<li>
We get the component of W<sub>1</sub>
perpendicular to normal N just by subtracting C<sub>1</sub>
from incoming ray W<sub>1</sub></sub>:
<font color=blue>
S<sub>1</sub> = W<sub>1</sub> - C<sub>1</sub>
</font>
<p>

<li>
Length of S<sub>1</sub> is the sine of incoming
angle &theta;<sub>1</sub>:
<font color=blue>
sin &theta;<sub>1</sub>
=
<big>|</big>S<sub>1</sub><big>|</big>
</font>
<p>

<li>
Snell's law tells us that:
<font color=blue>
<i>n<sub>1</sub></i> sin &theta;<sub>1</sub> =
<i>n<sub>2</sub></i> sin &theta;<sub>2</sub>
</font color=blue>
<p>

<li>
This gives the outgoing angle:
<font color=blue>
&theta;<sub>2</sub> = asin( sin &theta;<sub>1</sub> *
<i>n<sub>1</sub></i> /
<i>n<sub>2</sub></i> )
</font>
<p>

<li>
We can now construct the components of the outgoing ray:
<p>
<font color=blue>
C<sub>2</sub> = C<sub>1</sub> * cos &theta;<sub>2</sub> / cos &theta;<sub>1</sub>
<br>
S<sub>2</sub> = S<sub>1</sub> * sin &theta;<sub>2</sub> / sin &theta;<sub>1</sub>
</font>
<p>

<li>
W<sub>2</sub> is just the sum of those two components:
<font color=blue>
W<sub>2</sub> = C<sub>2</sub> + S<sub>2</sub>
</font>

</ol>

</blockquote>

\`);
setDescription(description + '<p><img src=imgs/refraction_steps.jpg width=470>'
                           + '<p><img src=imgs/glass.jpg width=169>');
`,

'2nd order':`
S.html(\`
<b>Second order surfaces:</b>
<blockquote>

We already know how to ray trace to spheres.
But we would like to ray trace to general second order surfaces,
like the red shapes you can see to the right.
To do this, we need to generalize from unit shapes
to more general shapes.
<p>
We have generalized from a unit sphere:
<blockquote>
x<sup>2</sup> +
y<sup>2</sup> +
z<sup>2</sup>
- 1 &le; 0
</blockquote>
to a general sphere:
<blockquote>
(x-C<sub>x</sub>)<sup>2</sup> +
(y-C<sub>y</sub>)<sup>2</sup> +
(z-C<sub>z</sub>)<sup>2</sup>
- r<sup>2</sup> &le; 0
</blockquote>
<p>
But suppose, for example, we want to ray trace to this set of equations,
which produce the blue cylinder below right:
<blockquote>
x<sup>2</sup> +
z<sup>2</sup>
- 1 &le; 0
&nbsp; <small><b>AND</b></small> &nbsp;
y &ge; -1
&nbsp; <small><b>AND</b></small> &nbsp;
y &le; 1
</blockquote>
We'd also like to ray trace to general cylinders,
including flat pill boxes and long skinny tubes,
in any location and direction.
<p>
To ray trace to all
general second order shapes,
including cones,
paraboloids and hyperboloids,
can employ matrices,
using a technique that was
first developed by Jim Blinn.

</blockquote>
\`);
setDescription(description + '<p><img src=imgs/second_order_surfaces.jpg width=470>'
                           + '<p><img src=imgs/cylinder.jpg width=265>');
`,

'coefficients':`
S.html(\`
<b>Second order surface coefficients:</b>
<blockquote>
Just as we can define any first order surface with 4 constants:
<blockquote>
<font color=red>C<sub>0</sub></font> x +
<font color=red>C<sub>1</sub></font> y +
<font color=red>C<sub>2</sub></font> z +
<font color=red>C<sub>3</sub></font> &nbsp; &le; &nbsp; 0
</blockquote>
we can define any second order surface with 10 constants:
<blockquote>
<font color=red>C<sub>0</sub></font> x<sup>2</sup> +
<font color=red>C<sub>1</sub></font> xy +
<font color=red>C<sub>2</sub></font> xz +
<font color=red>C<sub>3</sub></font> x
<br>
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
+ <font color=red>C<sub>4</sub></font> y<sup>2</sup>
+ <font color=red>C<sub>5</sub></font> yz
+ <font color=red>C<sub>6</sub></font> y
<br>
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
<small>&nbsp;
&nbsp;
&nbsp;</small>
+ <font color=red>C<sub>7</sub></font> z<sup>2</sup>
+ <font color=red>C<sub>8</sub></font> z
<br>
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
<small>&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;
&nbsp;</small>
+ <font color=red>C<sub>9</sub></font> &nbsp; &le; &nbsp; 0
</blockquote>
In the next section, we show how to represent those
10 constants by a 4x4 matrix.
<p>
This will let us do arbitrary linear transformations
on shapes bounded by second order surfaces.
</blockquote>
\`);
setDescription(description);
`,

'matrices':`
S.html(\`
<b>Using matrices to describe second order surfaces:</b>
<blockquote>

We can describe any second order surface as:
&nbsp;
<font color=blue>P<sup>T</sup></font> &bull; <font color=red>Q</font> &bull; <font color=blue>P</font> &le; 0
<p>

where: <font color=blue>P</font> is a point in space and
<br>
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <font color=red>Q</font> is a 4x4 matrix describing quadratic coefficients.
<p>

For example, the volume inside a unit sphere is described by:

<blockquote><table><tr>
<td><big><big><pre>
<font color=blue>x y z 1</font>
<td width=10>
<td><big><big><big>&bull;
<td width=10>
<th><big><big><pre>
<font color=red>1  0  0  0
0  1  0  0
0  0  1  0
0  0  0 -1</font>
<td width=10>
<td><big><big><big>&bull;
<td width=10>
<td><big><big><pre>
<font color=blue>x
y
z
1</font>
<td width=15>
<td><big><big><pre>&le;  0
<td width=15>
</tr></table></blockquote>

which multiplies out to:
&nbsp;
x<sup>2</sup> +
y<sup>2</sup> +
z<sup>2</sup> - 1 &le; 0
<p>

Suppose we have transformed <font color=blue>P</font> by matrix M, so: &nbsp;
<font color=blue>P</font>
&rarr;
<font color=blue>
M &bull; P</font>
<p>
<b>Problem:</b> Find a matrix that preserves the inequality.
<p>
<b>Hint:</b> Rewrite the transformed transpose <font color=blue>(M &bull; P)<sup>T</sup></font> as: &nbsp;
<font color=blue>P<sup>T</sup> &bull; M<sup>T</sup></font>
<p>
<b>Solution:</b> replace <font color=red>Q</font> by: &nbsp;
<font color=red>M<sup>-1<sup><big>T</big></sup></sup> &bull; Q &bull; M<sup>-1</sup></font>

<blockquote>
<font color=blue>P<sup>T</sup> &bull; M<sup>T</sup></font>
&bull;
<font color=red>M<sup>-1<sup><big>T</big></sup></sup> &bull; Q &bull; M<sup>-1</sup></font>
&bull;
<font color=blue>M &bull; P</font>
&le; 0
&nbsp;
</blockquote>

<b>Proof:</b> The above inequality simplifies to:
&nbsp; <font color=blue>P<sup>T</sup></font> &bull; <font color=red>Q</font> &bull; <font color=blue>P</font> &le; 0

</blockquote>
\`);
setDescription(description);
`,

'ray':`
S.html(\`
<b>Ray tracing to a general second order surface:</b>
<blockquote>

Let's say matrix Q contains coefficients of a quadratic shape:

<blockquote>
P<sup>T</sup>
&bull;
Q
&bull;
P
&le; 0
&nbsp; &nbsp; &nbsp; for any point P
</blockquote>
If we want to shoot a ray to the surface of this shape,
we need to evaluate this equation for points P = V+tW along the ray:

<blockquote>
(V+tW)<sup>T</sup>
&bull;
Q
&bull;
(V+tW)
&le; 0
</blockquote>

Because the left and right term each have two parts, there are four terms all together when we multiply this out:
<blockquote>
V<sup>T</sup>&bull;Q&bull;V
&nbsp;+&nbsp;
V<sup>T</sup>&bull;Q&bull;tW
&nbsp;+&nbsp;
tW<sup>T</sup>&bull;Q&bull;V
&nbsp;+&nbsp;
tW<sup>T</sup>&bull;Q&bull;tW
&nbsp; &le; &nbsp; 0
</blockquote>

Rearrange terms to express this as a quadratic equation in t.
Each of the terms shown in color evaluates to a single number:

<blockquote>
<font color=red>
(W<sup>T</sup>&bull;Q&bull;W)
</font>
t<sup>2</sup>

&nbsp;+&nbsp;

<font color=green>
(V<sup>T</sup>&bull;Q&bull;W + W<sup>T</sup>&bull;Q&bull;V)
</font>
t

&nbsp;+&nbsp;

<font color=blue>
V<sup>T</sup>&bull;Q&bull;V
</font>

&nbsp;&le;&nbsp;

0
</blockquote>

Solve this quadratic equation. Either we get two real roots, from which we can get both an entering and exiting point
for where the ray intersects the quadratic surface, or we get zero real roots, which indicates
that the ray has missed the surface.


</blockquote>
\`);
setDescription(description);
`,

'shapes':`
S.html(\`
<b>Varieties of unit quadratic shapes:</b>
<blockquote>

We have described the coefficients of a unit sphere via a 4x4 matrix.
We can also use a 4x4 matrix to describe the coefficients of other primitive shapes,
such as a unit cylindrical tube, cone, slab, paraboloid or hyperboloid,
then use 4x4 matrix transformations to translate, rotate or
scale.
<p>
We can also intersect shapes.
For example, to get a cylinder, transform a tube,
and a slab for the end caps.
Here are the coefficient matrices of some useful unit quadratic shapes:

<table>

<tr>
<td width=200><center><big><br><b>Slab</b><p>z<sup>2</sup> &le; 1<br>&nbsp;
<td width=40>
<td width=200><center><big><br><b>Tube</b><p>x<sup>2</sup>+y<sup>2</sup> &le; 1<br>&nbsp;
<td width=40>
<td width=200><center><big><br><b>Paraboloid</b><p>x<sup>2</sup>+y<sup>2</sup> &le; z<br>&nbsp;

<tr>

<th><big><pre>
 0  0  0  0
 0  0  0  0
 0  0  1  0
 0  0  0 -1
</pre>
<td>
<th><big><pre>
 1  0  0  0
 0  1  0  0
 0  0  0  0
 0  0  0 -1
</pre>
<td>
<th><big><pre>
 1  0  0  0
 0  1  0  0
 0  0  0 -1
 0  0  0  0
</pre>

<tr>
<td>
<small><small><small><small><p>&nbsp;<p>&nbsp;</small></small></small></small>

<tr>
<td width=200><center><big><b>Hyperboloid<br>of one sheet</b><p>x<sup>2</sup>+y<sup>2</sup>-z<sup>2</sup> &le; 1<br>&nbsp;
<td width=40>
<td width=200><center><big><br><b>Cone</b><p>x<sup>2</sup>+y<sup>2</sup>-z<sup>2</sup> &le; 0<br>&nbsp;
<td width=40>
<td width=200><center><big><b>Hyperboloid<br>of two sheets</b><p>x<sup>2</sup>+y<sup>2</sup>-z<sup>2</sup> &le; -1<br>&nbsp;

<tr>

<th><big><pre>
1  0  0  0
0  1  0  0
0  0 -1  0
0  0  0 -1
</pre>
<td>
<th><big><pre>
1  0  0  0
0  1  0  0
0  0 -1  0
0  0  0  0
</pre>
<td>
<th><big><pre>
1  0  0  0
0  1  0  0
0  0 -1  0
0  0  0  1
</pre>

</table>


</blockquote>
\`);
setDescription(description + \`
<p>&nbsp;<p>&nbsp;<p>&nbsp;<p>
<p>&nbsp;<p>&nbsp;<p><big>&nbsp;</big><br>
<img src=imgs/hyperboloid.png width=470>\`
);
`,


'normal':`
S.html(\`
<b>Computing the surface normal:</b>

<blockquote>
You can take the x,y,z partial derivatives
at surface point
<br>
P = [x,y,z], and then normalize the result, as follows:
<p>
When we multiply out the expression <font color=blue>P<sup>T</sup></font> &bull; <font color=red>Q</font> &bull; <font color=blue>P</font>:
<blockquote><table><tr>
<td><big><big><pre><font color=blue>x y z 1</font>
<td width=10>
<td><big><big><big>&bull;
<td width=10>
<td><big><big><pre>
<font color=red>a  b  c  d
e  f  g  h
i  j  k  l
m  n  o  p</font>
</pre>
<td width=10>
<td><big><big><big>&bull;
<td width=10>
<td><big><big><pre>
<font color=blue>x
y
z
1</font>
</tr></table></blockquote>
we get:
<blockquote>
<pre>
<font color=red>a</font> <font color=blue>x<sup>2</sup></font> + <font color=red>(b+e)</font> <font color=blue>xy</font> + <font color=red>(c+i)</font> <font color=blue>xz</font> + <font color=red>(d+m)</font> <font color=blue>x</font> +
         <font color=red>f</font>   <font color=blue>y<sup>2</sup></font> + <font color=red>(g+j)</font> <font color=blue>yz</font> + <font color=red>(h+n)</font> <font color=blue>y</font> +
                    <font color=red>k</font>   <font color=blue>z<sup>2</sup></font> + <font color=red>(l+o)</font> <font color=blue>z</font> +
                              <small><small> </small></small><font color=red>p</font>
</pre>
</blockquote>
If we take the partial derivatives with respect to x,y,z:
<blockquote>
<pre>
f<sub>x</sub> =   <font color=red>2a </font><font color=blue>x</font> + <font color=red>(b+e)</font><font color=blue>y</font> + <font color=red>(c+i)</font><font color=blue>z</font> + <font color=red>(d+m)</font>
f<sub>y</sub> = <font color=red>(b+e)</font><font color=blue>x </font>+   <font color=red>2f </font><font color=blue>y </font>+ <font color=red>(g+j)</font><font color=blue>z</font> + <font color=red>(h+n)</font>
f<sub>z</sub> = <font color=red>(c+i)</font><font color=blue>x </font>+ <font color=red>(g+j)</font><font color=blue>y </font>+   <font color=red>2k </font><font color=blue>z</font> + <font color=red>(l+o)</font>
</pre>
</blockquote>
the surface normal is:
<blockquote>
<code>normalize([f<sub>x</sub>,f<sub>y</sub>,f<sub>z</sub>])</code>
</blockquote>


</blockquote>
\`);
setDescription(description);
`,

'fog':`
S.html(\`
<b>Computing fog:</b>
<blockquote>
Fog consists of tiny water droplets suspended in the air.
If fog is uniform in density, then visibility
falls off exponentially with distance.
This is because
for every unit distance that light travels
through uniform fog, a
constant proportion of the light is scattered.
For example, if 50% of light is scattered
after the light travels a distance D through a uniform fog, then:
<ul>
<li>1/2 of the light will remain
unscattered after distance D;
<li>1/4 of the light will remain unscattered
after distance 2D;
<li>1/8 of the light will remain unscattered
after distance 3D.
</ul>

A simple minimal example of how to compute uniform fog in your fragment shader
is as follows:

<pre>
      <font color=blue>float f = pow(.8, tMin);         // FOG DENSITY
      color = mix(vec3(1.), color, f); // FOG COLOR</font>
</pre>

The above example is missing a number of things. For example,
<code><font color=blue>.8</font></code>
is a constant density. It is much better to create
a slider and let the user vary the fog density.
Also, the color of the fog is fixed to be white
by the use of <code><font color=blue>vec3(1.)</font></code>.
It would be useful to allow the user to use sliders to
vary the fog color.
<p>
Also, you can use the noise function
to vary both fog color and fog density in different parts of the
image, as seen to the right, as well as over time. Animated fog looks every cool. :-)

</blockquote>
\`);
setDescription(description + '<p><img src=imgs/foggy.png width=470>');
`,

};

}

