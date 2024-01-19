
rooms.raytrace5x = function() {

lib3D();

description = `Complex Quadric Equations and Shapes<br>
		<p>
		I tried to implement some shapes with hyperboloids. One of them is the concave tube, with a hyperboloid in the Z direction, and a qSlabZ. However, whenever the Z direction of the hyperboloid is facing the camera, the shape disappears. <br>
		The reasons for this problem is found. Since in our RT algorithm, we are calculating the last tIn and the first tOut, and determine whether the light hits the object using this method. This algorithm works for normal convex objects, because the last tIn will always be before the first tOut. This is not true for concave objects. Those objects will inevitably have some situations where one of the tIn is after a tOut. This is a major problem for concave objects, and RT algorithms need to be improved for concave objects. From a shape standpoint, the problem is also quite intuitive to understand. For all convex objects, a single light ray can only enter and then exit the shape once. This is not true for concave objects though. <br>
		For the concave tube I have here, it's a combination of two quadric surfaces. Normally, the two tIns and tOuts are sequenced like this: tIn, tIn, tOut, tOut. However, when the Z direction is facing the camera, the sequence becomes: tIn, tOut, tIn, tOut. <br>
		</p>
<small>
    <p>  
         <input type=range id=animationSpeed value=30 > animation speed
	<hr>
		 <input type=checkbox id=phongToggle checked> Phong Shading
	<br> <input type=checkbox id=reflectionToggle checked> Reflection
    <br> <input type=range id=hueCount       value=30 
                min=3      max=360                    >
		 <div id=hueCountInfo>&nbsp;</div>
    <br> <input type=range id=refract        value=50 > refract
    <br> <input type=range id=refractDiff    value=25 > color differential
         <div id=iorInfo>&nbsp;</div>
         <div id=refractDiffInfo>&nbsp;</div>
		Refraction indexes are not real-life. <br>
		They should be proportional to wavelength. 
	<hr>
         <input type=range id=fogRed         value=50 > fog red
    <br> <input type=range id=fogGreen       value=50 > fog green
    <br> <input type=range id=fogBlue        value=50 > fog blue
    <br> <input type=range id=fogDensity     value= 0 > fog density
	<hr>
	     <input type=range id=red            value=10 > bg red
    <br> <input type=range id=green          value=10 > bg green
    <br> <input type=range id=blue           value=10 > bg blue
</small>
`;

code = {
'init':`

	// DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL
	
	// ambient, diffuse, specular (with power)
	S.redPlastic    = [.2,.1,.1,0,  .5,.2,.2,0,  2,2,2,20,  0,0,0,0];
	S.greenPlastic  = [.1,.2,.1,0,  .2,.5,.2,0,  2,2,2,20,  0,0,0,0];
	S.bluePlastic   = [.1,.1,.2,0,  .2,.2,.5,0,  2,2,2,20,  0,0,0,0];
	S.whitePlastic  = [.2,.2,.2,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];
	S.clearPlastic  = [.1,.1,.1,0,   1, 1, 1,0,  1,1,1,20,  0,0,0,0];
`,

fragment: `
S.setFragmentShader(\`

	// DECLARE CONSTANTS, UNIFORMS, VARYING VARIABLES

	const int nQ = \` + S.nQ + \`;
	const int nL = \` + S.nL + \`;
	uniform vec3 uBgColor;
	uniform vec3 uLd[nL];
	uniform vec3 uLc[nL];
	uniform mat4 uQ[nQ];
	uniform mat4 uPhong[nQ];
	uniform int  uShape[nQ];
	uniform float uIor;
	uniform float uRefractDiff;
	uniform float uFogDensitySlider;
	uniform vec3 uFogColor;
	uniform float uTime;
	const float hueCount = float(\` + S.hueCount + \`);
	const bool phongToggle = \` + S.phongToggle + \`;
	const bool reflectionToggle = \` + S.reflectionToggle + \`;

	varying vec3 vPos;

	// DEFINE CAMERA FOCAL LENGTH

	float fl = 3.;

/********* PSEUDO-CODE IMPLEMENTATION OF FRAGMENT SHADER **********

Compute surface normal: P, Q => N

   vec3 normalQ(vec3 P, mat4 Q)

      Just like in the course notes.

Trace a ray to a quadric: V, W, Q => [ tIn, tOut ]

   vec2 rayQ(vec3 V, vec3 W, mat4 Q)

      Just like in the course notes:

         First add homogeneous coordinates:

            V1 = vec4(V,1)
            W0 = vec4(W,0)

         Then compute quadratic equation:

            a: W0 . Q*W0
            b: V1 . Q*W0 + W0 . Q*V1
            c: V1 . Q*V1

         Then solve quadratic equation.

         Return both roots as a vec2

Trace a ray to an intersection of quadric surfaces:

   Q1: T=[n,tIn,tOut], n, t.xy => T
   
      tIn  = t.x
      tOut = t.y
      if tIn > 0 and tIn < tOut and tIn < T.y
         T = [n,t]
      return T
   
   Q2: T=[n,tIn,tOut], n, t0.xy, t1.xy => T
   
      tIn  = max(t0.x,t1.x)
      tOut = min(t0.y,t1.y)
      if tIn > 0 and tIn < tOut and tIn < T.y
         i = t0.x==tIn ? 0 : 1
         T = [n+i,t]
      return T
   
   Q3: T=[n,tIn,tOut], n, t0.xy, t1.xy, t2.xy => T
   
      tIn  = max(t0.x,t1.x,t2.x)
      tOut = min(t0.y,t1.y,t2.y)
      if tIn > 0 and tIn < tOut and tIn < T.y
         i = t0.x==tIn ? 0 : t1.x==tIn ? 1 : 2
         T = [n+i,t]
      return T
   
   Q4: T=[n,tIn,tOut], n, t0.xy, t1.xy, t2.xy, t3.xy => T
   
      tIn  = max(t0.x,t1.x,t2.x,t3.x)
      tOut = min(t0.y,t1.y,t2.y,t3.y)
      if tIn > 0 and tIn < tOut and tIn < T.y
         i = t0.x==tIn ? 0 : t1.x==tIn ? 1 : t2.x==tIn ? 2 : 3
         T = [n+i,t]
      return T
   
Trace a ray to the scene:

   vec3 rayScene(vec3 V, vec3 W):

      T = [-1,1000,0]
      loop though all quadrics n
         if shape_type == 1: T = Q1(T, n, ray to Q[n])
         if shape_type == 2: T = Q2(T, n, ray to Q[n], Q[n+1])
         if shape_type == 3: T = Q3(T, n, ray to Q[n], Q[n+1], Q[n+2])
         if shape_type == 4: T = Q4(T, n, ray to Q[n], Q[n+1], Q[n+2], Q[n+3])
      return T

A note on using array subscripts that you computed within a function:

   Array subscripts need to be constant. So if you compute an array subscript
   within a function, then you need to make use of the resulting int value
   in the proper way, as follows:

   WRONG:

      vec3 T = rayScene(V, W);
      int n = int(T.x);
      mat4 Q = uQ[n];                 // This line produces a compiler error.
      ...

   CORRECT:

      vec3 T = rayScene(V, W);
      for (int n = 0 ; n < nQ ; n++)
         if (n == int(T.x)) {
            mat4 Q = uQ[n];          // This line works, because n is constant.
            ...
         }

Shade surface: P, N, phong => color

   vec3 shadeSurface(vec3 P, vec3 N, mat4 phong)

      The same algorithm you use when shading a sphere.

Refract ray: W, N, index_of_refraction => W

   vec3 refractRay(vec3 W, vec3 N, float n)

      Just like in the course notes.

Main loop

   T=[n,tIn,tOut] = rayScene(T, V, W)
   
   n = int(T.x) // REMEMBER, YOU NEED TO MAKE A LOOP HERE, AS SHOWN ABOVE.

   if n >= 0:

      Compute surface point P = V + T.y * W

      Shade with Phong, using:

         N = normalQ(P, Q[n])

      Do reflection:

         compute R

         T=[n,tIn,tOut] = rayScene ( P , R )

         n = int(T.x)

         if n >= 0:

            M = P + T.y * W        // POINT ON SURFACE OF OTHER OBJECT

               because T.y is tIn

            color += shadeSurface(M, normalQ(M,Q[n]), phong[n]) / 2.

      Do refraction:

         (1) SHOOT RAY TO FIND REAR OF THIS OBJECT (USE 2nd ROOT):

         W = refractRay(W, N, index_of_refraction)

         T=[n,tIn,tOut] = rayScene ( P-.01*W , W )

         n = int(T.x)

         P = P + T.z * W            // FIND POINT AT REAR OF OBJECT

            because T.z is tOut

         N = normalQ(P, Q[n])

         (2) SHOOT RAY FROM REAR OF THIS OBJECT TO FIND ANOTHER OBJECT:

         W = refract ray (W, N, 1 / index_of_refraction)

         T=[n,tIn,tOut] = rayScene( P , W )

         n = int(T.x)

         if n >= 0:

            M = P + T.y * W        // POINT ON SURFACE OF OTHER OBJECT

               because T.y is tIn

            color += diffuse_color_of_this_object *
                     shadeSurface(M, normalQ(M,Q[n]), phong[n])

******************************************************************/

	// compute surface normal
	// P is point on surface
	// Q is quadric matrix representing the transformed shape
	// implementing partial derivative from course notes
	vec3 normalQ(vec3 P, mat4 Q) {
		float fx = 2. * Q[0][0] * P.x + 
					(Q[0][1] + Q[1][0]) * P.y + 
					(Q[0][2] + Q[2][0]) * P.z + 
					(Q[0][3] + Q[3][0]);
		float fy = (Q[0][1] + Q[1][0]) * P.x + 
					2. * Q[1][1] * P.y + 
					(Q[1][2] + Q[2][1]) * P.z + 
					(Q[1][3] + Q[3][1]);
		float fz = (Q[0][2] + Q[2][0]) * P.x + 
					(Q[1][2] + Q[2][1]) * P.y + 
					2. * Q[2][2] * P.z + 
					(Q[2][3] + Q[3][2]);
		return normalize(vec3(fx, fy, fz));
	}

	// ray trace to a quadric surface
	// V is ray starting point
	// W is ray direction
	// Q is quadric matrix representing the transformed shape
	vec2 rayQ(vec3 V3, vec3 W3, mat4 Q) {
		vec4 V = vec4(V3, 1);
		vec4 W = vec4(W3, 0);
		// A point P is on surface if P^T dot Q dot P = 0
		// build quadratic equation
		float a = dot(W, Q * W);
		float b = dot(V, Q * W) + dot(W, Q * V);
		float c = dot(V, Q * V);
		// compute two roots
		if (b * b - 4. * a * c >= 0.) {
			// I forgot to check that the root exists at first.
			// The resulting rendering was surprising to say the least...
			// Attention to details...
			float smallRoot = (-b - sqrt(b * b - 4. * a * c)) / (2. * a);
			float bigRoot   = (-b + sqrt(b * b - 4. * a * c)) / (2. * a);
			return vec2(smallRoot, bigRoot);
		} else {
			return vec2(-1.);
		}
	}

	// functions for tracing a ray to an intersection of quadric surfaces
	// T is current guess
	// n is the index of quadric surface hit
	// t* is the two roots of ray trace results
	// finds the last tIn, first tOut
	// and sets surface index accordingly

	vec3 Q1(vec3 T, int n, vec2 t0) {
		float tIn  = t0.x;
		float tOut = t0.y;
		if (tIn > 0. && tIn < tOut && tIn < T.y) {
			T = vec3(n, t0);
		}
		return T;
	}

	vec3 Q2(vec3 T, int n, vec2 t0, vec2 t1) {
		float tIn  = max(t0.x, t1.x);
		float tOut = min(t0.y, t1.y);
		if (tIn > 0. && tIn < tOut && tIn < T.y) {
			int i = t0.x == tIn ? 0 : 1;
			T = vec3(n + i, tIn, tOut);
		}
		return T;
	}

	vec3 Q3(vec3 T, int n, vec2 t0, vec2 t1, vec2 t2) {
		float tIn  = max(t0.x, max(t1.x, t2.x));
		float tOut = min(t0.y, min(t1.y, t2.y));
		if (tIn > 0. && tIn < tOut && tIn < T.y) {
			int i = t0.x == tIn ? 0 : t1.x == tIn ? 1 : 2;
			T = vec3(n + i, tIn, tOut);
		}
		return T;
	}

	vec3 Q4(vec3 T, int n, vec2 t0, vec2 t1, vec2 t2, vec2 t3) {
		float tIn  = max(t0.x, max(t1.x, max(t2.x, t3.x)));
		float tOut = min(t0.y, min(t1.y, min(t2.y, t3.y)));
		if (tIn > 0. && tIn < tOut && tIn < T.y) {
			int i = t0.x == tIn ? 0 : t1.x == tIn ? 1 : t2.x == tIn ? 2 : 3;
			T = vec3(n + i, tIn, tOut);
		}
		return T;
	}

	// ray trace to the scene with all quadric surfaces
	vec3 rayScene(vec3 V, vec3 W) {
		vec3 T = vec3(-1, 10000000., -10000000.);
		for (int n = 0; n < nQ; n++) {
			if (uShape[n] == 1) T = Q1(T, n, rayQ(V, W, uQ[n]));
			if (uShape[n] == 2) T = Q2(T, n, rayQ(V, W, uQ[n]), rayQ(V, W, uQ[n + 1]));
			if (uShape[n] == 3) T = Q3(T, n, rayQ(V, W, uQ[n]), rayQ(V, W, uQ[n + 1]), rayQ(V, W, uQ[n + 2]));
			if (uShape[n] == 4) T = Q4(T, n, rayQ(V, W, uQ[n]), rayQ(V, W, uQ[n + 1]), rayQ(V, W, uQ[n + 2]), rayQ(V, W, uQ[n + 3]));
		}
		return T;
	}

	// use phong shading
	vec3 shadeSurface(vec3 P, vec3 N, mat4 phong) {
		
		// EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX

		vec3  ambient  = phong[0].rgb;
		vec3  diffuse  = phong[1].rgb;
		vec3  specular = phong[2].rgb;
		float p        = phong[2].a;

		// INIT COLOR, APPROXIMATE VECTOR TO EYE

		vec3 c = mix(ambient, uBgColor, .3);
		vec3 E = vec3(0., 0., 1.);

		// LOOP THROUGH LIGHT SOURCES
		
		for (int l = 0; l < nL; l++) {
			// COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE
			vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
			c += uLc[l] * (diffuse * max(0., dot(N, uLd[l])) 
					+ specular * pow(max(0., dot(R, E)), p));
		}
		// c *= 1. + .5 * noise(3. * N); // OPTIONAL SPOTTY TEXTURE
		return c;
	}

	// refraction using Snell's law
	vec3 refractRay(vec3 W1, vec3 N, float n) {
		vec3 C1 = N * dot(W1, N);
		vec3 S1 = W1 - C1;
		float sin1 = length(S1) / length(W1);
		float cos1 = sqrt(1. - sin1 * sin1);
		float theta2 = asin(sin1 / n);
		float sin2 = sin(theta2);
		float cos2 = cos(theta2);
		vec3 C2 = C1 * cos2 / cos1;
		vec3 S2 = S1 * sin2 / sin1;
		return C2 + S2;
	}

	vec3 hueToRgbPercent(float hue) {
		if (hue < 120.) {
			return vec3((120. - hue) / 120., hue / 120., 0.);
		} else if (hue >= 120. && hue < 240.) {
			return vec3(0., (240. - hue) / 120., (hue - 120.) / 120.);
		} else {
			return vec3((hue - 240.) / 120., 0., (360. - hue) / 120.);
		}
	}

	void main() {
		
		// BACKGROUND COLOR IS THE DEFAULT COLOR
		
		vec3 color = uBgColor;

		// DEFINE RAY INTO SCENE FOR THIS PIXEL

		vec3 V = vec3(0., 0., fl);
		vec3 W = normalize(vec3(vPos.xy, -fl));

		vec3 T = rayScene(V, W);
		
		// SEE WHICH QUADRIC SURFACE WE HIT
		for (int n = 0; n < nQ; n++) {
			if (n == int(T.x)) {
				// PHONG SHADING
				vec3 P = V + T.y * W;
				vec3 N = normalQ(P, uQ[n]);
				if (phongToggle) {
					color = shadeSurface(P, N, uPhong[n]);
				}
				
				// REFLECTION
				vec3 R = 2. * dot(N, -W) * N + W;
				vec3 TReflect = rayScene(P, R);
				for (int nReflect = 0; nReflect < nQ; nReflect++) {
					if (nReflect == int(TReflect.x)) {
						vec3 M = P + TReflect.y * R;
						vec3 NReflect = normalQ(M, uQ[nReflect]);
						if (reflectionToggle) {
							color += shadeSurface(M, NReflect, uPhong[nReflect]) / 2.;
						}
					}
				}

				// REFRACTION
				
				/*
				WInObject is light direction after the first refraction
				TInObject is ray trace results for WInObject. 
				nOut      is index of quadric surface that WInObject comes out on
				POut      is point of light coming out from the object
				NOut      is normal of POut at nOut
				WOut      is light direction after the second refraction
				TOut      is ray trace results for WOut
				nHit      is index of quadric surface that WOut hits (on another object)
				PHit      is point of light hitting the other object
				NHit      is normal of PHit at nHit
				*/

				/*
				for (int rgb = 0; rgb < 3; rgb++) {
					float refractDiff = uRefractDiff;
					float iorMultiplier = rgb == 0 ? refractDiff : rgb == 1 ? 1. : 1. / refractDiff;
					float thisIor = uIor * iorMultiplier;
					// 1. shoot ray to find the rear of this object
					vec3 WInObject = refractRay(W, N, thisIor);
					vec3 TInObject = rayScene(P - .01 * WInObject, WInObject);
					for (int nOut = 0; nOut < nQ; nOut++) {
						if (nOut == int(TInObject.x)) {
							vec3 POut = P + TInObject.z * WInObject;
							vec3 NOut = normalQ(POut, uQ[nOut]);
							// 2. shoot ray from the rear of this object to find another object
							vec3 WOut = refractRay(WInObject, NOut, 1. / thisIor);
							vec3 TOut = rayScene(POut, WOut);
							for (int nHit = 0; nHit < nQ; nHit++) {
								if (nHit == int(TOut.x)) {
									vec3 PHit = POut + TOut.y * WOut;
									vec3 NHit = normalQ(PHit, uQ[nHit]);
									vec3 colorComponent = rgb == 0 ? vec3(uPhong[nOut][1].r, 0., 0.) : 
															rgb == 1 ? vec3(0., uPhong[nOut][1].g, 0.) : 
																vec3(0., 0., uPhong[nOut][1].b);
									color += colorComponent * shadeSurface(PHit, NHit, uPhong[nHit]);
								}
							}
						}
					}
				}
				// END OF REFRACTION
				*/
				for (float hue = 0.; hue < 360.; hue += 360. / hueCount) {
					float thisIor = uIor * (uRefractDiff + (1. - uRefractDiff) * 2. * hue / 360.); // modifier centers around 1
//					float thisIor = uIor * (uRefractDiff + (1. - uRefractDiff) * 2. * hue / 360.) - (1. - uRefractDiff); // modifier max out at 1
					// 1. shoot ray to find the rear of this object
					vec3 WInObject = refractRay(W, N, thisIor);
					vec3 TInObject = rayScene(P - .01 * WInObject, WInObject);
					for (int nOut = 0; nOut < nQ; nOut++) {
						if (nOut == int(TInObject.x)) {
							vec3 POut = P + TInObject.z * WInObject;
							vec3 NOut = normalQ(POut, uQ[nOut]);
							// 2. shoot ray from the rear of this object to find another object
							vec3 WOut = refractRay(WInObject, NOut, 1. / thisIor);
							vec3 TOut = rayScene(POut, WOut);
							for (int nHit = 0; nHit < nQ; nHit++) {
								if (nHit == int(TOut.x)) {
									vec3 PHit = POut + TOut.y * WOut;
									vec3 NHit = normalQ(PHit, uQ[nHit]);
									vec3 rgbPercent = hueToRgbPercent(hue);
									vec3 colorComponent = vec3(uPhong[nOut][1].r * rgbPercent.r / hueCount,
															uPhong[nOut][1].g * rgbPercent.g / hueCount,
															uPhong[nOut][1].b * rgbPercent.b / hueCount);
									color += colorComponent * shadeSurface(PHit, NHit, uPhong[nHit]);
								}
							}
						}
					}
				}
				// END OF REFRACTION
			}
		}

		// FOG
		float density = .6 + (100. - uFogDensitySlider) / 100. * .4; // varies from .6 to 1.
		float f = pow(density, int(T.x) == -1 ? 5. : T.y);
		vec3 fogColor = uFogColor * (1. + .5 * noise(vPos * 1.5 + vec3(0., 0., .5 * uTime)));
		color = mix(fogColor, color, f);
		
		gl_FragColor = vec4(sqrt(color), 1.);
	}
\`);
`,
vertex: `
S.setVertexShader(\`

	attribute vec3 aPos;
	varying   vec3 vPos;

	void main() {
	  vPos = aPos;
	  gl_Position = vec4(aPos, 1.);
	}

\`)

`,
render: `

	// USEFUL VECTOR FUNCTIONS

	let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
	let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	let norm = v => Math.sqrt(dot(v,v));
	let normalize = v => { let s = norm(v); return [ v[0]/s, v[1]/s, v[2]/s ]; }
	let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
	let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];

	// SEND LIGHT SOURCE DATA TO GPU
	
	let ldData = [ normalize([1,1,1]),
				  normalize([-1,-1,-1])];
	S.setUniform('3fv', 'uLd', ldData.flat());
	S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);
	/*
	let ldData = [ normalize([0, 0, -1])];
	S.setUniform('3fv', 'uLd', ldData.flat());
	S.setUniform('3fv', 'uLc', [ 1,1,1 ]);
	*/
	// DEFINE NUMBER OF LIGHTS FOR GPU

	S.nL = ldData.length;

	// SEND BACKGROUND COLOR TO GPU

	S.setUniform('3fv', 'uBgColor', [ red.value   / 100,
									 green.value / 100,
									 blue.value  / 100 ]);

	// SEND INDEX OF REFRACTION TO GPU

	let ior = refract.value / 100 + 1;
	S.setUniform('1f', 'uIor', ior);

	let rd = 1 - refractDiff.value / 100 * .8;
	S.setUniform('1f', 'uRefractDiff', rd);

	// SEND HUE COUNT TO GPU
	
	S.hueCount = hueCount.value;

	// SEND TOGGLES TO GPU
	
	S.phongToggle = phongToggle.checked;
	S.reflectionToggle = reflectionToggle.checked;

	// SEND FOG DATA TO GPU
	
	S.setUniform('1f', 'uFogDensitySlider', fogDensity.value);
	S.setUniform('3fv', 'uFogColor', [fogRed.value / 100,
										fogGreen.value / 100,
										fogBlue.value / 100]);

	// SEND ANIMATION TIME AND SPEED TO GPU
	
	S.setUniform('1f', 'uTime', time);

	// DIFFERENT QUADRIC SURFACES

	//                xx        yy         zz           c

	let qSlabX  = [1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,-1]; // x*x - 1 <= 0
	let qSlabY  = [0,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,-1]; // y*y - 1 <= 0
	let qSlabZ  = [0,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,-1]; // z*z - 1 <= 0
	let qSphere = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]; // x*x + y*y + z*z - 1 <= 0
	let qTubeX  = [0,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1]; // y*y + z*z - 1 <= 0
	let qTubeY  = [1,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,-1]; // x*x + z*z - 1 <= 0
	let qTubeZ  = [1,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,-1]; // x*x + y*y - 1 <= 0

	let qPrismBottom = [0,0,0,0, 0,1,0,-Math.sqrt(3)/6, 0,0,0,0, 0,0,0,-1/6];
	let qPrismLeft   = [3,-2*Math.sqrt(3),0,-1, 0,1,0,Math.sqrt(3)/3, 0,0,0,0, 0,0,0,-5/48];
	let qPrismRight  = [3,2*Math.sqrt(3),0,1, 0,1,0,Math.sqrt(3)/3, 0,0,0,0, 0,0,0,-5/48];

	let qHyperboloidX = [-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-.5];
	let qHyperboloidY = [1,0,0,0, 0,-1,0,0, 0,0,1,0, 0,0,0,-.5];
	let qHyperboloidZ = [1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,-.5];

	// SHAPES ARE INTERSECTIONS OF QUADRIC SURFACES

	let shape = [], coefs = [], xform = [], phong = [], M;

	let sphere = (m, M) => {
	  shape.push(1);
	  phong.push(m);
	  xform.push(M);
	  coefs.push(qSphere);
	}

	let tubeX = (m, M) => {
	  shape.push(2, 0);
	  phong.push(m, m);
	  xform.push(M, M);
	  coefs.push(qTubeX, qSlabX);
	}

	let tubeY = (m, M) => {
	  shape.push(2, 0);
	  phong.push(m, m);
	  xform.push(M, M);
	  coefs.push(qTubeY, qSlabY);
	}

	let tubeZ = (m, M) => {
	  shape.push(2, 0);
	  phong.push(m, m);
	  xform.push(M, M);
	  coefs.push(qTubeZ, qSlabZ);
	}

	let cube = (m, M) => {
	  shape.push(3, 0, 0);
	  phong.push(m, m, m);
	  xform.push(M, M, M);
	  coefs.push(qSlabX, qSlabY, qSlabZ);
	}

	let curvedCube = (m, M) => {
	  shape.push(3, 0, 0);
	  phong.push(m, m, m);
	  xform.push(M, M, M);
	  coefs.push(qTubeX, qTubeY, qTubeZ);
	}

	let octahedron = (m, M) => {
	  shape.push(4, 0, 0, 0);
	  phong.push(m, m, m, m);
	  xform.push(M, M, M, M);
	  coefs.push([1, 2, 2, 0,  0, 1, 2, 0,  0,0,1,0,  0,0,0,-1]);
	  coefs.push([1,-2,-2, 0,  0, 1, 2, 0,  0,0,1,0,  0,0,0,-1]);
	  coefs.push([1,-2, 2, 0,  0, 1,-2, 0,  0,0,1,0,  0,0,0,-1]);
	  coefs.push([1, 2,-2, 0,  0, 1,-2, 0,  0,0,1,0,  0,0,0,-1]);
	}

	let prism = (m, M) => {
		shape.push(4, 0, 0, 0);
		phong.push(m, m, m, m);
		xform.push(M, M, M, M);
		coefs.push(qSlabZ, qPrismBottom, qPrismLeft, qPrismRight);
	}

	let concaveTube = (m, M) => {
		shape.push(2, 0);
		phong.push(m, m);
		xform.push(M, M);
		coefs.push(qHyperboloidZ, qSlabZ);
	}

	let star = (m, M) => {
		shape.push(3, 0, 0);
		phong.push(m, m, m);
		xform.push(M, M, M);
		coefs.push(qHyperboloidX, qHyperboloidY, qHyperboloidZ);
	}

	time *= animationSpeed.value / 100; // to see reflection and refraction more clearly
	
	// CREATE MY OWN SCENE
	/*
	prism(S.whitePlastic,
		mScale(.9, .9, .25,
		mRoty(Math.PI / 2,
		mRotz(0,
		mRotx(time,
		matrixTranslate(0, .5, 0))))));
	*/
	/*
	octahedron(S.whitePlastic,
		mScale(.25, .25, .25,
		mRoty(time * 1.1,
		mRotz(time * 1.2,
		mRotx(time * 1.3,
		matrixTranslate(0, 0, 2))))));
	*/
	
	concaveTube(S.clearPlastic,
		mScale(.1, .1, .1,
		mRoty(time * -1.3,
		mRotz(time * -1.1,
		mRotx(time * -1.2,
		matrixTranslate(-.3, 0, 1))))));

	star(S.clearPlastic,
		mScale(.1, .1, .1,
		mRoty(time * -1.1,
		mRotz(time * -1.2,
		mRotx(time * -1.3,
		matrixTranslate(.3, 0, 1))))));
	/*
	curvedCube(S.clearPlastic,
		mScale(.4, .4, .4,
		mRoty(-time * -1.3,
		mRotz(-time * -1.1,
		mRotx(-time * -1.2,
		matrixTranslate(0, 0, -1))))));
	*/
	/*
	cube(S.whitePlastic,
		mScale(.1, 1, 1,
		mRoty(Math.PI / 4,
		matrixTranslate(0, 0, 0))));
	*/
	sphere(S.bluePlastic,
		mScale(.05, .05, .05,
		matrixTranslate(Math.cos(time) * .3, Math.sin(time) * .3, 2)));

	// SEND SCENE DATA TO GPU

	for (let n = 0 ; n < coefs.length ; n++) {
	  let IM = matrixInverse(xform[n]);
	  coefs[n] = matrixMultiply(matrixTranspose(IM), matrixMultiply(coefs[n], IM));
	}
	S.setUniform('1iv', 'uShape', shape);
	S.setUniform('Matrix4fv', 'uQ', false, coefs.flat());
	S.setUniform('Matrix4fv', 'uPhong', false, phong.flat());

	// DEFINE NUMBER OF QUADRIC SURFACES FOR GPU

	S.nQ = coefs.length;

	// RENDER THIS ANIMATION FRAME

	S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);

	// SET ANY HTML INFO

	iorInfo.innerHTML = 'index of refraction = ' + (ior * 100 >> 0) / 100;
	let redComponent = (rd * 100 >> 0) / 100;
	let greenComponent = 1;
	let blueComponent = (1 / rd * 100 >> 0) / 100;
	refractDiffInfo.innerHTML = 'coefficient for rgb components: ' + redComponent + ', ' + greenComponent + ', ' + blueComponent;
	hueCountInfo.innerHTML = 'Hue Count: ' + hueCount.value;
`,
events: `
	;
`
};

}


