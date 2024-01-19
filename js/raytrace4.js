
rooms.raytrace4 = function() {

lib3D();

description = `Raytracing to spheres<br>in a fragment shader
<small>
    <p>
    <b>Background color</b>
    <br> <input type=range id=red   value= 5> red
    <br> <input type=range id=green value=10> green
    <br> <input type=range id=blue  value=50> blue
</small>
<p>
Here's a pyramid I created with shadow implementation.
<br><br> 
The central octahedron rotates and moves up and down, and shadows on the lower pyramid is visible. 
<br><br>
I tried to do shadows for the sphere, but the shading function for spheres and polyhedrons are separate, so I haven't figured out how to do this yet.
<br><br>
This pyramid design is inspired from the game Genshin Impact. Here's an in-game screenshot of the pyramid. 
<br>
<img src='imgs/genshin.png', width=600>
</p>
`;

code = {
'init':`

   // DEFINE NUMBER OF SPHERES AND NUMBER OF LIGHTS

   S.nS = 1;
   S.nL = 2;

   // DEFINE MATERIALS TO BE RENDERED VIA PHONG REFLECTANCE MODEL

   let materials = [
      [.15,.05,.025,0, .3,.1,.05,0, .6,.2,.1,3, 0,0,0,0], // COPPER
      [.25,.15,.025,0, .5,.3,.05,0, 1,.6,.1,6,  0,0,0,0], // GOLD
      [.25,0,0,0,      .5,0,0,0,    2,2,2,20,   0,0,0,0], // PLASTIC
      [.05,.05,.05,0,  .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // LEAD
      [.1,.1,.1,0,     .1,.1,.1,0,  1,1,1,5,    0,0,0,0], // SILVER
   ];

   S.material = [];
   for (let n = 0 ; n < S.nS ; n++)
      S.material.push(materials[n % materials.length]);

   S.polyhedronMaterial = [];
   S.polyhedronMaterial.push(materials[3]);
   S.polyhedronMaterial.push(materials[4]);
   S.polyhedronMaterial.push(materials[1]);

   // INITIALIZE SPHERE POSITIONS AND VELOCITIES

   S.sPos = [];
   S.sVel = [];
   for (let n = 0 ; n < S.nS ; n++) {
      /*
      S.sPos.push([ Math.random() - .5,
                    Math.random() - .5,
                    Math.random() - .5 ]);
	  */
	  S.sPos.push([0,0,0]);
      S.sVel.push([0,0,0]);
   }
`,
fragment: `
S.setFragmentShader(\`

   // DECLARE CONSTANTS, UNIFORMS, VARYING VARIABLES

   const int nS = \` + S.nS + \`;
   const int nL = \` + S.nL + \`;
   uniform float uTime;
   uniform vec3 uBgColor;
   uniform vec4 uS[nS];
   uniform mat4 uSm[nS];
   uniform mat4 uIMCube1;
   uniform mat4 uIMCube2;
   uniform mat4 uIMOct;
   uniform mat4 uCube1Material;
   uniform mat4 uCube2Material;
   uniform mat4 uOctahedronMaterial;
   uniform vec4 uCube[6];
   uniform vec4 uOctahedron[8];
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];

   varying vec3 vPos;

   // DEFINE CAMERA FOCAL LENGTH

   float fl = 3.;

   // TRACE A RAY TO A HALFSPACE

   float rayHalfspace(vec3 V, vec3 W, vec4 H) {
      vec4 V1 = vec4(V, 1.);
      vec4 W0 = vec4(W, 0.);
      return -dot(H, V1) / dot(H, W0);
   }

   // TRACE A RAY TO A TRANSFORMED CUBE

   vec4 rayCube(vec3 V, vec3 W, mat4 IM) {
      vec3 N = vec3(0.);
      float tIn = -1000., tOut = 1000.;
      for (int i = 0 ; i < 6 ; i++) {
         vec4 H = uCube[i] * IM;
	 H /= sqrt(dot(H.xyz, H.xyz));
	 float t = rayHalfspace(V, W, H);
	 if (dot(W, H.xyz) < 0.) {
	    if (t > tIn)
	       N = H.xyz;
	    tIn = max(tIn, t);
	 }
	 else
	    tOut = min(tOut, t);
      }
      return vec4(N, tIn < tOut ? tIn : -1.);
   }

   vec4 rayOctahedron(vec3 V, vec3 W, mat4 IM) {
      vec3 N = vec3(0.);
      float tIn = -1000., tOut = 1000.;
      for (int i = 0 ; i < 8 ; i++) {
         vec4 H = uOctahedron[i] * IM;
	 H /= sqrt(dot(H.xyz, H.xyz));
	 float t = rayHalfspace(V, W, H);
	 if (dot(W, H.xyz) < 0.) {
	    if (t > tIn)
	       N = H.xyz;
	    tIn = max(tIn, t);
	 }
	 else
	    tOut = min(tOut, t);
      }
      return vec4(N, tIn < tOut ? tIn : -1.);
   }

   // TRACE A RAY TO A SPHERE

   float raySphere(vec3 V, vec3 W, vec4 S) {
      V -= S.xyz;
      V += .01 * W;
      float b = dot(V, W);
      float d = b * b - dot(V, V) + S.w * S.w;
      return d < 0. ? -1. : -b - sqrt(d);
   }

   // SHADE A SPHERE AT ONE SURFACE POINT

   vec3 shadeSphere(vec3 P, vec4 S, mat4 M) {

      // EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX

      vec3  ambient  = M[0].rgb;
      vec3  diffuse  = M[1].rgb;
      vec3  specular = M[2].rgb;
      float p        = M[2].a;

      // COMPUTE NORMAL, INIT COLOR, APPROXIMATE VECTOR TO EYE

      vec3 N = normalize(P - S.xyz);
      vec3 c = mix(ambient, uBgColor, .3);
      vec3 E = vec3(0.,0.,1.);

      // LOOP THROUGH LIGHT SOURCES

      for (int l = 0 ; l < nL ; l++) {

         // ADD DIFFUSE AND SPECULAR ONLY IF NOT IN SHADOW

         float t = -1.;
         for (int n = 0 ; n < nS ; n++)
            t = max(t, raySphere(P, uLd[l], uS[n]));

         // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE

         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                         + specular * pow(max(0., dot(R, E)), p));
         }
      }
      //c *= 1. + .5 * noise(3.*N); // OPTIONAL SPOTTY TEXTURE
      return c;
   }

   vec3 shadeHalfspace(vec3 P, vec3 N, mat4 M) {

      // EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX

      vec3  ambient  = M[0].rgb;
      vec3  diffuse  = M[1].rgb;
      vec3  specular = M[2].rgb;
      float p        = M[2].a;

      // COMPUTE INIT COLOR, APPROXIMATE VECTOR TO EYE

      vec3 c = mix(ambient, uBgColor, .3);
      vec3 E = vec3(0.,0.,1.);

      // LOOP THROUGH LIGHT SOURCES

      for (int l = 0 ; l < nL ; l++) {

         // ADD DIFFUSE AND SPECULAR ONLY IF NOT IN SHADOW

         float t = -1.;

		 // COMPUTE t FOR THREE POLYHEDRONS
		 t = max(t, rayOctahedron(P, uLd[l], uIMOct).w);
         t = max(t, rayCube(P, uLd[l], uIMCube1).w);
		 t = max(t, rayCube(P, uLd[l], uIMCube2).w);

         // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE

         if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                         + specular * pow(max(0., dot(R, E)), p));
         }
      }
      //c *= 1. + .5 * noise(3.*N); // OPTIONAL SPOTTY TEXTURE
      return c;
   }

   void main() {

      // BACKGROUND COLOR IS THE DEFAULT COLOR

      vec3 color = uBgColor;

      // DEFINE RAY INTO SCENE FOR THIS PIXEL

      vec3 V = vec3(0.,0.,fl);
      vec3 W = normalize(vec3(vPos.xy, -fl));

      // LOOP THROUGH SPHERES

      float tMin = 10000.;
      for (int n = 0 ; n < nS ; n++) {
         float t = raySphere(V, W, uS[n]);
         if (t > 0. && t < tMin) {

            // IF THIS IS THE NEAREST SPHERE, DO PHONG SHADING

            vec3 P = V + t * W;
            color = shadeSphere(P, uS[n], uSm[n]);
            tMin = t;

            // SEE WHETHER ANY OTHER SPHERE IS VISIBLE VIA REFLECTION

            vec3 N = normalize(P - uS[n].xyz);
            vec3 R = 2. * dot(N, -W) * N + W;
            float rtMin = 10000.;
            vec3 rColor;
            for (int rn = 0 ; rn < nS ; rn++) {
               float rt = raySphere(P, R, uS[rn]);
               if (rt > 0. && rt < rtMin) {
                  rtMin = rt;
                  rColor = shadeSphere(P + rt * R, uS[rn], uSm[rn]);
               }
            }
            if (rtMin < 10000.)
               color += .5 * rColor;
         }
      }

      // THIS SECTION COULD BE REPLACED WITH A LOOP

      // RAY TRACE TO THE OCTAHEDRON
      vec4 Nt = rayOctahedron(V, W, uIMOct);
      if (Nt.w > 0. && Nt.w < tMin) {
	  	 tMin = Nt.w;
         vec3 a = mix(vec3(.1), uBgColor, .3);
         color = a + vec3(max(0., dot(Nt.xyz, vec3(.5))));
		 color = shadeHalfspace(V + Nt.w * W, Nt.xyz, uOctahedronMaterial); 
      }
	  
	  // RAY TRACE TO CUBE1
	  vec4 NtCube1 = rayCube(V, W, uIMCube1);
      if (NtCube1.w > 0. && NtCube1.w < tMin) {
	     tMin = NtCube1.w;
         vec3 a = mix(vec3(.1), uBgColor, .3);
         color = a + vec3(max(0., dot(NtCube1.xyz, vec3(.5))));
		 color = shadeHalfspace(V + NtCube1.w * W, NtCube1.xyz, uCube1Material); 
      }

	  // RAY TRACE TO CUBE2
	  vec4 NtCube2 = rayCube(V, W, uIMCube2);
      if (NtCube2.w > 0. && NtCube2.w < tMin) {
	     tMin = NtCube2.w;
         vec3 a = mix(vec3(.1), uBgColor, .3);
         color = a + vec3(max(0., dot(NtCube2.xyz, vec3(.5))));
		 color = shadeHalfspace(V + NtCube2.w * W, NtCube2.xyz, uCube2Material); 
      }

      // SET PIXEL COLOR

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

   // HANDY DANDY VECTOR LIBRARY

   let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
   let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   let norm = v => Math.sqrt(dot(v,v));
   let normalize = v => {
      let s = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      return [ v[0]/s, v[1]/s, v[2]/s ];
   }
   let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
   let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];

   // DEFINE RADIUS OF EACH SPHERE

   let radius = .07;

   // SEND LIGHT SOURCE DATA TO GPU

   let ldData = [ normalize([1,1,-1]),
                  normalize([-1,-1,1]) ];
   S.setUniform('3fv', 'uLd', ldData.flat());
   S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);

   // SEND ANIMATION TIME TO GPU

   S.setUniform('1f', 'uTime', time);

   // MOVE SPHERES INTO PLACE FOR THIS ANIMATION FRAME

   for (let n = 0 ; n < S.nS ; n++) {
      /*
      for (let i = 0 ; i < 3 ; i++) {
         S.sVel[n][i] += .003 * Math.cos(time + (2+i) * n);
         S.sVel[n][i] += .01 * (Math.random() - .5);
         S.sPos[n][i] += .1 * S.sVel[n][i];
      }
	  */
	  S.sPos[n] = [.3 * Math.cos(time), 0., .3 * Math.sin(time)];
      //S.sPos[n] = scale(normalize(S.sPos[n]), .7);
	  console.log(S.sPos[n]);
   }

   // AVOID SPHERE INTERPENETRATION

   for (let m = 0 ; m < S.nS ; m++)
   for (let n = 0 ; n < S.nS ; n++)
      if (m != n) {
         let D = subtract(S.sPos[m], S.sPos[n]);
         let d = norm(D);
         if (d < 2 * radius) {
            let t = 2 * radius - d;
            for (let i = 0 ; i < 3 ; i++) {
               S.sPos[m][i] += t * D[i] / d;
               S.sPos[n][i] -= t * D[i] / d;
            }
         }
      }

   // SEND SPHERES DATA TO GPU
   let sData = [];
   for (let n = 0 ; n < S.nS ; n++)
      sData.push(S.sPos[n], radius);
   S.setUniform('4fv', 'uS', sData.flat());
   S.setUniform('Matrix4fv', 'uSm', false, S.material.flat());

   // SEND POLYHEDRON MATERIAL TO GPU
   S.setUniform('Matrix4fv', 'uCube1Material', false, S.polyhedronMaterial[0]);
   S.setUniform('Matrix4fv', 'uCube2Material', false, S.polyhedronMaterial[1]);
   S.setUniform('Matrix4fv', 'uOctahedronMaterial', false, S.polyhedronMaterial[2]);

   // SEND BACKGROUND COLOR TO GPU

   S.setUniform('3fv', 'uBgColor', [ red.value  /100,
                                     green.value/100,
                                     blue.value /100 ]);

   // BEGINNINGS OF A MATRIX LIBRARY

   let matrixIdentity = () => [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];
   let matrixTranslate = (x,y,z) => [1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1];
   let matrixScale = (x,y,z) => [x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1];
   let matrixRotx = t => {
      let c = Math.cos(t), s = Math.sin(t);
      return [1,0,0,0, 0,c,s,0, 0,-s,c,0, 0,0,0,1];
   }
   let matrixRoty = t => {
      let c = Math.cos(t), s = Math.sin(t);
      return [c,0,-s,0, 0,1,0,0, s,0,c,0, 0,0,0,1];
   }
   let matrixRotz = t => {
      let c = Math.cos(t), s = Math.sin(t);
      return [c,s,0,0, -s,c,0,0, 0,0,1,0, 0,0,0,1];
   }

   let cM;

   // RENDER THE OCTAHEDRON

   cM = matrixIdentity();
   //cM = matrixMultiply(cM, matrixTranslate(Math.cos(time)/2,Math.sin(time)/2,.5));
   //cM = matrixMultiply(cM, matrixRotx(time));
   //cM = matrixMultiply(cM, matrixRotz(time));
   cM = matrixMultiply(cM, matrixRoty(time)); // rotate around the yaw axis
   cM = matrixMultiply(cM, matrixScale(.2,.15,.15)); // set size
   cM = matrixMultiply(cM, matrixTranslate(0, Math.sin(time) / 5, 0)); // floating movement
   S.setUniform('Matrix4fv', 'uIMOct', false, matrixInverse(cM));

   // RENDER THE LOWER FLAT PYRAMID

   cM = matrixIdentity();
   //cM = matrixMultiply(cM, matrixRoty(time));
   cM = matrixMultiply(cM, matrixRoty(Math.PI / 6)); // rotate to see more faces
   cM = matrixMultiply(cM, matrixTranslate(0., -.3, 0.)); // set location
   cM = matrixMultiply(cM, matrixScale(.3, .3, .3)); // set size
   S.setUniform('Matrix4fv', 'uIMCube1', false, matrixInverse(cM));

   // RENDER THE UPPER FLAT PYRAMID

   cM = matrixIdentity();
   cM = matrixMultiply(cM, matrixRoty(Math.PI / 6)); // rotate to see more faces
   cM = matrixMultiply(cM, matrixRotx(Math.PI)); // flip the pyramid
   cM = matrixMultiply(cM, matrixTranslate(0., -.23, 0.)); // set location
   cM = matrixMultiply(cM, matrixScale(.15, .15, .15)); // set size
   cM = matrixMultiply(cM, matrixTranslate(0, -Math.sin(time) / 2.5, 0)); // floating movement
   S.setUniform('Matrix4fv', 'uIMCube2', false, matrixInverse(cM));

   // MODIFIED CUBE HALFSPACES
   // LOOKS LIKE A FLAT PYRAMID
   S.setUniform('4fv', 'uCube', [
      -1,1,0,-1, 1,1,0,-1,
      0,-1,0,-1, 0,1,0,-.5,
      0,1,-1,-1, 0,1,1,-1,
   ]);

   S.setUniform('4fv', 'uOctahedron', [
      -1,-1,-1,-1, 1,1,1,-1,
      1,-1,-1,-1, -1,1,1,-1,
      -1,1,-1,-1, 1,-1,1,-1,
      -1,-1,1,-1, 1,1,-1,-1,
   ]);

	S.setUniform('4fv', 'uPrism', [
		-1/4, Math.sqrt(3)/12, 0, -1, 
		 1/4, Math.sqrt(3)/12, 0, -1, 
		0, -1, 0, -1, 
		0, 0, 1, -1, 
		0, 0, -1, -1,
	]);

   // RENDER THIS ANIMATION FRAME

   S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
`,
events: `
   ;
`
};

}


