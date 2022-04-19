#version 330

uniform mat4 m_pvm;
uniform mat4 m_viewModel;
uniform mat3 m_normal;
uniform mat4 m_Model;
uniform mat4 m_View;

uniform int texMode;

uniform vec4 l_pos;
uniform vec4 l_pos2;
uniform vec4 l_pos3;
uniform vec4 l_pos4;
uniform vec4 l_pos5;
uniform vec4 l_pos6;
uniform vec4 l_posG;

uniform vec3 offsets[89];

uniform bool flag_instance;
uniform bool flag_bumpmap;

in vec4 position;
in vec4 normal;
in vec4 texCoord;
in vec4 tangent;

out float visibility;

out Data {
	vec3 normal;
	vec3 eye;
	vec3 lightDir;
	vec3 lightDir2;
	vec3 lightDir3;
	vec3 lightDir4;
	vec3 lightDir5;
	vec3 lightDir6;
	vec3 lightDirG;
	vec4 position;
	vec2 tex_coord;
	vec4 pos;
	//vec2 sphere_coord;
	vec3 skyboxTexCoord;
} DataOut;

const float density = 0.01f;
const float gradient = 1.5f;

void main () {
	vec3 n, t, b;
	vec3 lightDir, eyeDir;
	vec3 aux1;
	vec4 aux2;

	DataOut.skyboxTexCoord = vec3(m_Model * position);
	DataOut.skyboxTexCoord.x = - DataOut.skyboxTexCoord.x;
	DataOut.tex_coord = texCoord.st;

	vec4 pos = m_viewModel * position;


	
	DataOut.position = pos;
	
	DataOut.lightDir2 = vec3(l_pos2 - pos);
	DataOut.lightDir3 = vec3(l_pos3 - pos);
	DataOut.lightDir4 = vec3(l_pos4 - pos);
	DataOut.lightDir5 = vec3(l_pos5 - pos);
	DataOut.lightDir6 = vec3(l_pos6 - pos);
	DataOut.lightDirG = vec3(l_posG - pos);

	lightDir = vec3(l_pos - pos);
	n = normalize(m_normal * normal.xyz);
	eyeDir = vec3(-pos);
	

	if(flag_bumpmap) {  //convert eye and light vectors to tangent space

		//Calculate components of TBN basis in eye space
		t = normalize(m_normal * tangent.xyz);  
		b = tangent.w * cross(n,t);

		aux1.x = dot(lightDir, t);
		aux1.y = dot(lightDir, b);
		aux1.z = dot(lightDir, n);
		lightDir = normalize(aux1);

		aux1.x = dot(eyeDir, t);
		aux1.y = dot(eyeDir, b);
		aux1.z = dot(eyeDir, n);
		eyeDir = normalize(aux1);
	}

	DataOut.normal = n;
	DataOut.lightDir = lightDir;
	DataOut.eye = eyeDir;

	vec4 poscam = pos;

	
	if (flag_instance){
		vec3 offset = offsets[gl_InstanceID];
		aux2 =  vec4(vec3(position) + offset, 1.0);
	}
	else{
		aux2 = position;
	}

	gl_Position = m_pvm * aux2;

	float distance = length(poscam.xyzw);
	visibility = exp(-pow((distance * density) , gradient));
	visibility = clamp(visibility, 0.0f, 1.0f);
}