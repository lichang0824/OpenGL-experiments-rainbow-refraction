
rooms.presentation = function() {

code = {

'Color Components':`
S.html(\`
<b>Separating Colors for Refraction</b><br> 
<textarea rows='20' cols='80' readonly>
// Refract with rainbow
for (hue = 0; hue <= 360; hue += step) {
	// calculate IOR for this hue
	// do RT refraction with this IOR
	// calculate RGB percents for this hue
	color += phong * rgbPercents / hueCount
}

// Original refract
// do RT refraction with general IOR
color = phong
</textarea>
\`);
setDescription('');
`,

'RGB and HSV':`
S.html(\`
<b>Calculating RGB Percents</b><br>
<textarea rows='20' cols='80' readonly>
vec3 hueToRgbPercent(float hue) {
	if (hue < 120.) {
		return vec3((120. - hue) / 120., hue / 120., 0.);
	} else if (hue >= 120. && hue < 240.) {
		return vec3(0., (240. - hue) / 120., (hue - 120.) / 120.);
	} else {
		return vec3((hue - 240.) / 120., 0., (360. - hue) / 120.);
	}
}
</textarea>
\`);
setDescription('<p><img src=imgs/conversion.png width=550><p><img src=imgs/color_wheel.png width=550>');
`,

'Quadric Surfaces':`
S.html(\`
<b>Experimenting with Quadric Surfaces</b>
<textarea rows='20' cols='80' readonly>
let qPrismBottom = [0,0,0,0, 0,1,0,-Math.sqrt(3)/6, 0,0,0,0, 0,0,0,-1/6];
let qPrismLeft   = [3,-2*Math.sqrt(3),0,-1, 0,1,0,Math.sqrt(3)/3, 0,0,0,0, 0,0,0,-5/48];
let qPrismRight  = [3,2*Math.sqrt(3),0,1, 0,1,0,Math.sqrt(3)/3, 0,0,0,0, 0,0,0,-5/48];

let qHyperboloidX = [-1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-.5];
let qHyperboloidY = [1,0,0,0, 0,-1,0,0, 0,0,1,0, 0,0,0,-.5];
let qHyperboloidZ = [1,0,0,0, 0,1,0,0, 0,0,-1,0, 0,0,0,-.5];
</textarea>
\`);
setDescription('<p><img src=imgs/prism.jpg width=200><img src=imgs/hyperboloid_one.png width=200><p><img src=imgs/triangle.png width=550><img src=imgs/qPrismLeft.png width=550>');
`,

'Interactions':`
S.html(\`
Toggle Phong shading<br>
Hue Count Slider<br>
Refract Parameters
\`);
setDescription('');
`,

};

}

