#version 330

uniform sampler2D texmap;
uniform sampler2D texmap1;
uniform sampler2D texmap2;


uniform int texMode;

uniform bool spotlight_mode;
uniform vec4 spotDir;
uniform float spotCosCutOff;

out vec4 colorOut;

struct Materials {
	vec4 diffuse;
	vec4 ambient;
	vec4 specular;
	vec4 emissive;
	float shininess;
	int texCount;
};
uniform Materials mat;


in Data {
	vec3 normal;
	vec3 eye;
	vec3 lightDir;
	vec2 tex_coord;
} DataIn;

void main() {

	vec4 texel, texel1;

	vec4 spec = vec4(0.0);
	float intensity = 0.0f;
	float intSpec = 0.0f;

	vec3 n = normalize(DataIn.normal);
	vec3 l = normalize(DataIn.lightDir);
	vec3 e = normalize(DataIn.eye);
	vec3 sd = normalize(vec3(-spotDir));

	if(spotlight_mode == true)  {     //Scene iluminated by a spotlight
		if(dot(l, sd) > spotCosCutOff)  {				//inside cone?
			intensity = max(dot(n,l), 0.0);
			if (intensity > 0.0) {
				vec3 h = normalize(l + e);
				intSpec = max(dot(h,n), 0.0);
				spec = mat.specular * pow(intSpec, mat.shininess);
			}
		}
	}
	else {				//Scene iluminated by a pointlight
		intensity = max(dot(n,l), 0.0);
		if (intensity > 0.0) {
			vec3 h = normalize(l + e);
			intSpec = max(dot(h,n), 0.0);
			spec = mat.specular * pow(intSpec, mat.shininess);
		}
	}

	if(texMode == 0) // modulate diffuse color with texel color
	{
		texel = texture(texmap2, DataIn.tex_coord);  // texel from lighwood.tga
		colorOut = max(intensity * mat.diffuse * texel + spec,mat.ambient * texel);
	}
	else if (texMode == 2) // diffuse color is replaced by texel color, with specular area or ambient (0.1*texel)
	{
		texel = texture(texmap, DataIn.tex_coord);  // texel from stone.tga
		colorOut = max(intensity*texel + spec, 0.1*texel);
	}
	else if (texMode == 1)  // multitexturing
	{
		texel = texture(texmap2, DataIn.tex_coord);  // texel from lighwood.tga
		texel1 = texture(texmap1, DataIn.tex_coord);  // texel from checker.tga
		colorOut = max(intensity*texel*texel1 + spec, 0.1*texel*texel1);
	}
	else // modulated texture for particle
	{
		texel = texture(texmap, DataIn.tex_coord);  //texel from element flare texture
		if((texel.a == 0.0)  || (mat.diffuse.a == 0.0) ) discard;
		else
			colorOut = mat.diffuse * texel;
	}
}
