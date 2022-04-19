#version 330

out vec4 colorOut;

uniform sampler2D texmap0;
uniform sampler2D texmap1; //when bumpamp this holds the normal
uniform samplerCube cubeMap;

uniform bool spotlight_mode;
uniform bool directional;
uniform vec4 coneDir;
uniform float spotCosCutOff;
uniform bool shadowMode;

uniform int texMode;

struct Materials {
	vec4 diffuse;
	vec4 ambient;
	vec4 specular;
	vec4 emissive;
	float shininess;
	int texCount;
};

uniform Materials mat;
uniform vec4 fogColor;

uniform bool text_flag;
uniform int type_text_flag;
//uniform bool bill_flag;
//uniform bool sky_flag;


in Data {
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
} DataIn;

in float visibility;

void main() {

	vec3 color = vec3(0.0);

	vec4 texel0, texel1;

	vec3 n;

	vec4 spec = vec4(0.0);

	float intensity = 0.0f;
	float intensity2 = 0.0f;
	float intensity3 = 0.0f;
	float intensity4 = 0.0f;
	float intensity5 = 0.0f;
	float intensity6 = 0.0f;
	float intensityG = 0.0f;

	float intSpec = 0.0f;
	float intSpec2 = 0.0f;
	float intSpec3 = 0.0f;
	float intSpec4 = 0.0f;
	float intSpec5 = 0.0f;
	float intSpec6 = 0.0f;

	float att = 0.0;
	float att2 = 0.0;
	float att3 = 0.0;
	float att4 = 0.0;
	float att5 = 0.0;
	float att6 = 0.0;
	float attG = 0.0;
	float spotExp = 80.0;

	if(type_text_flag == 5){
		n = normalize(2.0 * texture(texmap1, DataIn.tex_coord).rgb - 1.0);
	}else n = normalize(DataIn.normal);


	vec3 lightDir = DataIn.lightDir - DataIn.position.xyz;
	vec3 lightDir2 = DataIn.lightDir2 - DataIn.position.xyz;
	vec3 lightDir3 = DataIn.lightDir3 - DataIn.position.xyz;
	vec3 lightDir4 = DataIn.lightDir4 - DataIn.position.xyz;
	vec3 lightDir5 = DataIn.lightDir5 - DataIn.position.xyz;
	vec3 lightDir6 = DataIn.lightDir6 - DataIn.position.xyz;
	vec3 lightDirG = DataIn.lightDirG - DataIn.position.xyz;
	vec3 l = normalize(lightDir);
	vec3 l2 = normalize(lightDir2);
	vec3 l3 = normalize(lightDir3);
	vec3 l4 = normalize(lightDir4);
	vec3 l5 = normalize(lightDir5);
	vec3 l6 = normalize(lightDir6);

	float light_distance = length(lightDir);
	float light_distance2 = length(lightDir2);
	float light_distance3 = length(lightDir3);
	float light_distance4 = length(lightDir4);
	float light_distance5 = length(lightDir5);
	float light_distance6 = length(lightDir6);
	float light_distanceG = length(lightDirG);
	vec3 e = normalize(DataIn.eye-DataIn.position.xyz);
	vec3 sd = normalize(vec3(-coneDir));

	vec3 half_vector = normalize(l+e);
	vec3 half_vector2 = normalize(l2+e);
	vec3 half_vector3 = normalize(l3+e);
	vec3 half_vector4 = normalize(l4+e);
	vec3 half_vector5 = normalize(l5+e);
	vec3 half_vector6 = normalize(l6+e);

	float diffuse = max(0.0, dot(n, l));//intensity diff
	float specular = max(0.0, dot(n, half_vector));//intensity spec

	float diffuse2 = max(0.0, dot(n, l2));//intensity diff
	float specular2 = max(0.0, dot(n, half_vector2));//intensity spec

	float diffuse3 = max(0.0, dot(n, l3));//intensity diff
	float specular3 = max(0.0, dot(n, half_vector3));//intensity spec

	float diffuse4 = max(0.0, dot(n, l4));//intensity diff
	float specular4 = max(0.0, dot(n, half_vector4));//intensity spec

	float diffuse5 = max(0.0, dot(n, l5));//intensity diff
	float specular5 = max(0.0, dot(n, half_vector5));//intensity spec

	float diffuse6 = max(0.0, dot(n, l6));//intensity diff
	float specular6 = max(0.0, dot(n, half_vector6));//intensity spec


	if (diffuse == 0) {
		specular = 0;
	} else {
		specular = pow (specular, mat.shininess);
	}
	if (diffuse2 == 0) {
		specular2 = 0;
	} else {
		specular2 = pow (specular2, mat.shininess);
	}

	if (diffuse3 == 0) {
		specular3 = 0;
	} else {
		specular3 = pow (specular3, mat.shininess);
	}
	if (diffuse4 == 0) {
		specular4 = 0;
	} else {
		specular4 = pow (specular4, mat.shininess);
	}
	if (diffuse5 == 0) {
		specular5 = 0;
	} else {
		specular5 = pow (specular5, mat.shininess);
	}
	if (diffuse6 == 0) {
		specular6 = 0;
	} else {
		specular6 = pow (specular6, mat.shininess);
	}

	if(shadowMode){  //constant color
		colorOut = vec4(0.5, 0.5, 0.5, 1.0);
	}
	else{
		if (text_flag){

			if(type_text_flag == 1) { //multitexturing
				texel0 = texture(texmap0, DataIn.tex_coord);
				texel1 = texture(texmap1, DataIn.tex_coord);

				if (spotlight_mode == true)  { //Scene iluminated by a spotlight
					//		colorOut = vec4(normalize(DataIn.position).xyz,1.0);
					//		return;

					att = 1.0f/(0.001f*light_distance+0.001f*light_distance*light_distance+0.1);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse * mat.diffuse.rgb * texel0.rgb* texel1.rgb *att;
					color += specular * mat.specular.rgb * texel0.rgb * texel1.rgb*att;

					att2 = 1.0f/(0.001f*light_distance2+0.001f*light_distance2*light_distance2+0.1f);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse2 * mat.diffuse.rgb * texel0.rgb* texel1.rgb *att2;
					color += specular2 * mat.specular.rgb * texel0.rgb* texel1.rgb *att2;

					att3 = 1.0f/(0.001f*light_distance3+0.001f*light_distance3*light_distance3+0.1f);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse3 * mat.diffuse.rgb * texel0.rgb * texel1.rgb*att3;
					color += specular3 * mat.specular.rgb * texel0.rgb * texel1.rgb*att3;

					att4 = 1.0f/(0.001f*light_distance4+0.001f*light_distance4*light_distance4+0.1f);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse4 * mat.diffuse.rgb * texel0.rgb * texel1.rgb*att4;
					color += specular4 * mat.specular.rgb * texel0.rgb * texel1.rgb*att4;

					att5 = 1.0f/(0.001f*light_distance5+0.001f*light_distance5*light_distance5+0.1f);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse5 * mat.diffuse.rgb * texel0.rgb * texel1.rgb*att5;
					color += specular5 * mat.specular.rgb * texel0.rgb * texel1.rgb*att5;

					att6 = 1.0f/(0.001f*light_distance6+0.001f*light_distance6*light_distance6+0.1f);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse6 * mat.diffuse.rgb * texel0.rgb * texel1.rgb*att6;
					color += specular6 * mat.specular.rgb * texel0.rgb * texel1.rgb*att6;

					colorOut = vec4(min(max(color, 0.07* texel0.rgb* texel1.rgb), 1.0), mat.diffuse.a);
				}
				else if (directional == true) {
					vec3 lightDirG = DataIn.lightDirG - DataIn.position.xyz;
					vec3 half_vectorG = normalize(lightDirG+e);

					float diffuseG = max(0.0, dot(n, lightDirG));
					float specularG = max(0.0, dot(n, half_vectorG));

					if (diffuseG == 0.0)
					specularG = 0.0;
					else
					specularG = pow(specularG, mat.shininess);

					attG = 1.0f/(0.01f*light_distanceG+0.01f*light_distanceG*light_distanceG+0.1);
					//intensity = max(dot(n,l), 0.0) ;
					color = diffuseG * mat.diffuse.rgb * texel0.rgb* texel1.rgb *attG;
					color += specularG * mat.specular.rgb * texel0.rgb * texel1.rgb*attG;
					colorOut = vec4(min(max(color, 0.07* texel0.rgb* texel1.rgb), 1.0), mat.diffuse.a);
					
					
				}
				else {
					att = 1.0f/(0.001f*light_distance+0.001f*light_distance*light_distance+0.1);
					//intensity = max(dot(n,l), 0.0) ;
					color += diffuse * mat.diffuse.rgb * texel0.rgb* texel1.rgb *att;
					color += specular * mat.specular.rgb * texel0.rgb * texel1.rgb*att;

					colorOut = vec4(min(max(color, 0.07* texel0.rgb* texel1.rgb), 1.0), mat.diffuse.a);

				}
			}
			else if(type_text_flag == 2){ //particules
				texel0 = texture(texmap0, DataIn.tex_coord);
				if ((texel0.a == 0.0) || (mat.diffuse.a == 0.0) ){
					discard;
				}
				else {
					colorOut = mat.diffuse * texel0;
				}
			}
			else if(type_text_flag == 3){ //billboard
				texel0 = texture(texmap0, DataIn.tex_coord);  
			
				if(texel0.a == 0.0){ 
					discard;
				}
				else{
					vec3 h = normalize(l + e);
					vec3 spec = vec3(0.0);
					intensity = max(dot(n,l), 0.0);
					intSpec = max(dot(h,n), 0.0);
					spec = mat.specular.rgb * pow(intSpec, mat.shininess);
					colorOut = vec4(max(intensity*texel0.rgb + spec, 0.1*texel0.rgb), texel0.a);
				}
			}
			else if(type_text_flag == 4){ //skybox
				//colorOut = texture(cubeMap, DataIn.skyboxTexCoord);
			}
			else if(type_text_flag == 5 || type_text_flag == 6){ 
				float intensity = max(dot(n,l), 0.0);
				if (intensity > 0.0) {
					vec3 h = normalize(l + e);
					float intSpec = max(dot(h,n), 0.0);
					spec = mat.specular * pow(intSpec, mat.shininess);
				}
				//n = normalize(2.0 * texture(texmap1, DataIn.tex_coord).rgb - 1.0);
				texel0 = texture(texmap0, DataIn.tex_coord);  // texel from stone.tga
				colorOut = vec4((max(intensity*texel0 + spec, 0.2*texel0)).rgb, 1.0);

			}
			
				
		}
		else {
			if (spotlight_mode == true)  { //Scene iluminated by a spotlight
				//		colorOut = vec4(normalize(DataIn.position).xyz,1.0);
				//		return;

				att = 1.0f/(0.001f*light_distance+0.001f*light_distance*light_distance+0.1);
				//intensity = max(dot(n,l), 0.0) ;
				color += diffuse * mat.diffuse.rgb *att;
				color += specular * mat.specular.rgb *att;

				att2 = 1.0f/(0.001f*light_distance2+0.001f*light_distance2*light_distance2+0.1f);
				//intensity = max(dot(n,l), 0.0) ;
				color += diffuse2 * mat.diffuse.rgb  *att2;
				color += specular2 * mat.specular.rgb  *att2;

				att3 = 1.0f/(0.001f*light_distance3+0.001f*light_distance3*light_distance3+0.1f);
				//intensity = max(dot(n,l), 0.0) ;
				color += diffuse3 * mat.diffuse.rgb*att3;
				color += specular3 * mat.specular.rgb *att3;

				att4 = 1.0f/(0.001f*light_distance4+0.001f*light_distance4*light_distance4+0.1f);
				//intensity = max(dot(n,l), 0.0) ;
				color += diffuse4 * mat.diffuse.rgb *att4;
				color += specular4 * mat.specular.rgb *att4;

				att5 = 1.0f/(0.001f*light_distance5+0.001f*light_distance5*light_distance5+0.1f);
				//intensity = max(dot(n,l), 0.0) ;
				color += diffuse5 * mat.diffuse.rgb *att5;
				color += specular5 * mat.specular.rgb * att5;

				att6 = 1.0f/(0.001f*light_distance6+0.001f*light_distance6*light_distance6+0.1f);
				//intensity = max(dot(n,l), 0.0) ;
				color += diffuse6 * mat.diffuse.rgb *att6;
				color += specular6 * mat.specular.rgb *att6;

				colorOut = vec4(color, mat.diffuse.a);
			}
			else if (directional == true) {
				vec3 lightDirG = DataIn.lightDirG - DataIn.position.xyz;
				vec3 half_vectorG = normalize(lightDirG+e);

				float diffuseG = max(0.0, dot(n, lightDirG));
				float specularG = max(0.0, dot(n, half_vectorG));

				if (diffuseG == 0.0)
					specularG = 0.0;
				else
					specularG = pow(specularG, mat.shininess);

				color = diffuseG * mat.diffuse.rgb;
				color += specularG * mat.specular.rgb;
				colorOut = vec4(color, mat.diffuse.a);
			}
			if(texMode==3) {
				texel0 = texture(texmap0, DataIn.tex_coord);  //texel from element flare texture
				if((texel0.a == 0.0)  || (mat.diffuse.a == 0.0) ) 
					discard;
				else
					colorOut = mat.diffuse * texel0;
					
			}
			else {
			//intensity = max(dot(n,l), 0.0) ;

			color += diffuse * mat.diffuse.rgb;
			color += specular * mat.specular.rgb;
			colorOut = vec4(color, mat.diffuse.a);
			}
			
		}

	}
		//colorOut = mix(fogColor , colorOut, visibility); //fog

}
