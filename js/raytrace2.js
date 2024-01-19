
rooms.raytrace2 = function() {

lib3D();

description = `Raytracing to spheres<br>in a fragment shader
    <p>  <input type=range id=red   > red
    <br> <input type=range id=green > green
    <br> <input type=range id=blue  > blue
`;

code = {
'init':`
   S.nS = 5;
   S.nL = 2;

   S.sPos = [];
   for (let n = 0 ; n < S.nS ; n++)
      for (let i = 0 ; i < 3 ; i++)
         S.sPos.push(Math.random() - .5);
`,
fragment: `
S.setFragmentShader(\`
   const int nS = \` + S.nS + \`;
   const int nL = \` + S.nL + \`;
   uniform float uTime;
   uniform vec4 uS[nS];
   uniform vec3 uLd[2];
   uniform vec3 uLc[2];
   uniform vec3 uDc;

   varying vec3 vPos;

   float fl = 3.;

   vec3 ld = normalize(vec3(sin(10.*uTime),1.,1.));

   float raySphere(vec3 V, vec3 W, vec4 S) {
      V -= S.xyz;
      float b = dot(V, W);
      float d = b * b - dot(V, V) + S.w * S.w;
      return d < 0. ? -1. : -b - sqrt(d);
   }

   vec3 shadeSphere(vec3 P, vec4 S) {
      vec3 N = normalize(P - S.xyz);
      vec3 c = vec3(.1,.1,.2);
      vec3 E = vec3(0.,0.,1.);
      for (int l = 0 ; l < nL ; l++) {

         float t = -1.;
	 for (int n = 0 ; n < nS ; n++)
	    t = max(t, raySphere(P, uLd[l], uS[n]));

         float weight = step(t, 0.);
         //if (t < 0.) {
            vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
            c += weight * uLc[l] * (uDc * max(0.,dot(N, uLd[l]))
	                   + pow(max(0., dot(R, E)), 20.));
         //}
      }
      //c *= 1. + .5 * noise(3.*N);
      return c;
   }

   void main() {
      vec3 color = vec3(.1,.2,.5);
      vec3 V = vec3(0.,0.,fl);
      vec3 W = normalize(vec3(vPos.xy, -fl));

      float tMin = 10000.;
      for (int n = 0 ; n < nS ; n++) {
         float t = raySphere(V, W, uS[n]);
         if (t > 0. && t < tMin) {
            color = shadeSphere(V + t * W, uS[n]);
	    tMin = t;
         }
      }

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

   let normalize = v => {
      let s = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      return [ v[0]/s, v[1]/s, v[2]/s ];
   }

   let ld0 = normalize([1,1,1]);
   let ld1 = normalize([-1,-1,-1]);

   let ldData = [];
   for (let i = 0 ; i < 3 ; i++)
      ldData.push(ld0[i]);
   for (let i = 0 ; i < 3 ; i++)
      ldData.push(ld1[i]);

   S.setUniform('3fv', 'uLd', ldData);
   S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);
   S.setUniform('3fv', 'uDc', [ red.value/100,
                                green.value/100,
                                blue.value/100 ]);

   S.setUniform('4fv', 'uS', [ 0,0,0,.5 ]);
   S.setUniform('1f', 'uTime', time);

   let sData = [];
   let k = 0;
   for (let n = 0 ; n < S.nS ; n++) {
      for (let i = 0 ; i < 3 ; i++)
         sData.push(S.sPos[k++]);
      sData.push(.3);
   }

   S.setUniform('4fv', 'uS', sData);

   S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, 4);
`,
events: `
   ;
`
};

}


