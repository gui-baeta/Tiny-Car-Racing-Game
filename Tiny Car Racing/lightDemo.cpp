//
// AVT demo light 
// based on demos from GLSL Core Tutorial in Lighthouse3D.com   
//
// This demo was built for learning purposes only.
// Some code could be severely optimised, but I tried to
// keep as simple and clear as possible.
//
// The code comes with no warranties, use it at your own risk.
// You may use it, or parts of it, wherever you want.
//

#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <ratio>
#include <chrono>
#include <filesystem>
#include <cmath>

// include GLEW to access OpenGL 3.3 functions
#include <GL/glew.h>


#include <IL/il.h>


// GLUT is the toolkit to interface with the OS
#include <GL/freeglut.h>

// Use Very Simple Libs
#include "VSShaderlib.h"
#include "AVTmathLib.h"
#include "VertexAttrDef.h"
#include "geometry.h"
#include "Scene.h"
#include "avtFreeType.h"


using namespace std;

#define CAPTION "Micro Machines Halloween"

bool oranges_kill = false; // TODO change this back to true

#define frand()			((float)rand()/RAND_MAX)

int WindowHandle = 0;
int WinX = 640, WinY = 480;


unsigned int FrameCount = 0;

VSShaderLib shader, shaderText;

//Vector with meshes
vector<struct MyMesh> myMeshes;

//External array storage defined in AVTmathLib.cpp

/// The storage for matrices
extern float mMatrix[COUNT_MATRICES][16];
extern float mCompMatrix[COUNT_COMPUTED_MATRICES][16];

/// The normal matrix
extern float mNormal3x3[9];


inline double clamp(const double x, const double min, const double max) {
	return (x < min ? min : (x > max ? max : x));
}

inline int clampi(const int x, const int min, const int max) {
	return (x < min ? min : (x > max ? max : x));
}


GLint pvm_uniformId;
GLint vm_uniformId;
GLint model_uniformId;
GLint normal_uniformId;
GLint view_uniformId;
GLint shadowMode_uniformId;
GLint lPos_uniformId;
GLint lPos2_uniformId;
GLint lPos3_uniformId;
GLint lPos4_uniformId;
GLint lPos5_uniformId;
GLint lPos6_uniformId;
GLint lPosG_uniformId;
GLint tex_loc, tex_loc1, tex_loc2;
GLint texMode_uniformId;
float lightScreenPos[3]; 
GLuint FlareTextureArray[5];
FLARE_DEF AVTflare; 

GLuint TextureArray[3];

float alphaAux = 39.0f, betaAux = 21.0f;
float rAux = 10.0f;

// Camera Spherical Coordinates
float _alpha = 39.0f, _beta = 21.0f;
float r = 10.0f;

float _alphaStatic = 39.0f, _betaStatic = 51.0f;
float rStatic = 10.0f;
	
// Camera Position
float camX, camY, camZ;

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = 39.0f, beta = 51.0f;

// Frame counting and FPS computation
long myTime,timebase = 0,frame = 0;
char s[32];
float lightPos[4] = { 0.0f, 30.0f, 0.0f, 1.0f };

float lightPos2[4] = { 50.0f, 30.0f, 15.0f, 1.0f };
float lightPos3[4] = { 70.0f, 30.0f, 20.0f, 1.0f };
float lightPos4[4] = { -60.0f, 30.0f, -12.0f, 1.0f };
float lightPos5[4] = { -12.0f, 30.0f, -60.0f, 1.0f };
float lightPos6[4] = { 100.0f, 30.0f, 0.0f, 1.0f };
float lightPosG[4] = { 0.0f, 30.0f, 0.0f, 1.0f };
float lightPosR[4] = {4.0f, 7.0f, 2.0f, 1.0f};
bool spotlight_mode = false;
bool directional = false;
bool flareEffect = false;

const int WIDTH = 1366;
const int HEIGHT = 768;

int amount = 89; // Cheerios
int max_particles = 1500; 
bool snow = false;

std::vector<std::vector<float>> cheercoor(amount);

std::vector<struct Particle> particles(max_particles);

extern inline float rand_uf(float top);

unsigned long long n_frames = 0;

double last_time = 0.0f;

double t_delta = 0;

bool _first_time = true;
bool _reset = false;
bool _to_reset = false;

Scene scene;

std::chrono::high_resolution_clock::time_point part_timer;
std::chrono::high_resolution_clock::time_point part_timer_end;


bool car_forwards = false;
bool car_backwards = false;
bool car_left = false;
bool car_right = false;

char font_name[] = "fonts/arial.ttf";


const char* VERT = "./shaders/pointlight.vert";
const char* FRAG = "./shaders/pointlight.frag";
const char* ASSET1 = "./stone.tga";
const char* ASSET2 = "./checker.png";
const char* ASSET3 = "./lightwood.tga";

#define STARTING_LIFES 5;
int lifes_left = STARTING_LIFES;
std::string str_life = "HP:";
std::string str_points = "Points:";
double points_n = 0;

bool _game_over = false;
bool _pause = false;

// Normalize a vec3
std::vector<float> normalize_vec(std::vector<float> a) {
	float mag = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

	a[0] /= mag;
	a[1] /= mag;
	a[2] /= mag;

	return a;
}

void roadSetUp() {
	for (size_t i = 0; i < amount; i++) {
		if (i > 0 && i < 45) {
			cheercoor.at(i) = vec_f({ static_cast<float>(18 * cos(i) * (3 + cos(4 * i))), 0.0f,
								  static_cast<float>(18 * sin(i) * (3 + cos(4 * i))) });
		}
		if (i >= 45 && i < amount) {
			cheercoor.at(i) = vec_f({ static_cast<float>(13 * cos(i) * (3 + cos(4 * i))), 0.0f,
								  static_cast<float>(13 * sin(i) * (3 + cos(4 * i))) });
		}
		if (i == 0) {
			cheercoor.at(i) = vec_f({ 0.0f, -100.0f, 0.0f });
		}
	}
}

void iniParticles() {
	GLfloat v, theta, phi;

	for (auto&& particle : particles) {
		v = 0.8 * frand() + 0.2;
		phi = frand() * PI;
		theta = 2.0 * frand() * PI;

		particle.position = vec_f({0.0f, 10.0f, 0.0f});
		particle.vx = v * cos(theta) * sin(phi);
		particle.vy = v * cos(phi);
		particle.vz = v * sin(theta) * sin(phi);
		particle.ax = 0.1f; /* simular um pouco de vento */
		particle.ay = -0.15f; 
		particle.az = 0.0f;

		/* tom amarelado que vai ser multiplicado pela textura que varia entre branco e preto */


		particle.r = 0.882f;
		particle.g = 0.552f;
		particle.b = 0.211f;

		particle.life = rand_uf(3.0);		/* vida inicial */
		particle.fade = 0.0025f;	    
	}

	scene.get_entities()->at(EntityType::SNOW).get()->get_objects()->at(0).set_mesh_diffuse_alpha(1.0f);
}

void updateParticles() {
	int i;
	float h;

	h = 0.125f;
	//h = 0.033;
	if (snow) {

		for (auto&& particle : particles) {
			auto x = particle.position.at(0) + (h * particle.vx);
			auto y = particle.position.at(1) + (h * particle.vy);
			auto z = particle.position.at(2) + (h * particle.vz);
			particle.position = vec_f({x , y, z});
			particle.vx += (h*particle.ax);
			particle.vy += (h*particle.ay);
			particle.vz += (h*particle.az);
			//particle.life -= particle.fade;
		}
	}
}

void render_skybox(VSShaderLib* shader, Object* obj) {
	GLuint loc;
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(shader->getProgramIndex());

	//loc = glGetUniformLocation(shader->getProgramIndex(), "bill_flag");
    //glUniform1i(loc, false);
	loc = glGetUniformLocation(shader->getProgramIndex(), "text_flag");
    glUniform1i(loc, true);
	//loc = glGetUniformLocation(shader->getProgramIndex(), "particle_flag");
    //glUniform1i(loc, false);
	loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
    glUniform1i(loc, 4);


	glActiveTexture(GL_TEXTURE0);
	std::vector<GLuint>* _textureArray = obj->get_textureArray();
    glBindTexture(GL_TEXTURE_CUBE_MAP, _textureArray->at(0));

	loc = glGetUniformLocation(shader->getProgramIndex(), "texMode");
    glUniform1i(loc, 3);
    loc = glGetUniformLocation(shader->getProgramIndex(), "cubeMap");
    glUniform1i(loc, 0);

	loadIdentity(MatrixTypes::VIEW);
	loadIdentity(MatrixTypes::MODEL);

	glDepthMask(GL_FALSE); 
	glFrontFace(GL_CW);

	pushMatrix(MatrixTypes::MODEL);
	pushMatrix(MatrixTypes::VIEW);

	mMatrix[MatrixTypes::VIEW][12] = 0.0f;
	mMatrix[MatrixTypes::VIEW][13] = 0.0f;
	mMatrix[MatrixTypes::VIEW][14] = 0.0f;

	obj->transform();

	glUniformMatrix4fv(model_uniformId, 1, GL_FALSE, mMatrix[MatrixTypes::MODEL]); 
	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::PROJ_VIEW_MODEL]);

	glBindVertexArray(obj->get_mesh()->vao);
	glDrawElements(obj->get_mesh()->type, obj->get_mesh()->numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	popMatrix(MatrixTypes::MODEL);
	popMatrix(MatrixTypes::VIEW);
	
	glFrontFace(GL_CCW); // restore counter clockwise vertex order to mean the front
	glDepthMask(GL_TRUE);

	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	glutSwapBuffers();

}

void draw_mirror(VSShaderLib* shader, Object* obj){
	GLuint loc;

	loc = glGetUniformLocation(shader->getProgramIndex(), "text_flag");
    glUniform1i(loc, false);
	//loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
    //glUniform1i(loc, 1);

	//std::vector<GLuint>* _textureArray = obj->get_textureArray();
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, _textureArray->at(0));
	//glActiveTexture(GL_TEXTURE1);
	//glBindTexture(GL_TEXTURE_2D, _textureArray->at(1));

	// loc = glGetUniformLocation(shader->getProgramIndex(), "texMode");
	// glUniform1i(loc, 1);
	//loc = glGetUniformLocation(shader->getProgramIndex(), "texmap0");
	//glUniform1i(loc, 0);
	//loc = glGetUniformLocation(shader->getProgramIndex(), "texmap1");
	//glUniform1i(loc, 1);

	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.ambient");
	glUniform4fv(loc, 1, obj->get_mesh()->mat.ambient);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.diffuse");
	glUniform4fv(loc, 1, obj->get_mesh()->mat.diffuse);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.specular");
	glUniform4fv(loc, 1, obj->get_mesh()->mat.specular);
	loc = glGetUniformLocation(shader->getProgramIndex(), "mat.shininess");
	glUniform1f(loc, obj->get_mesh()->mat.shininess);

	loc = glGetUniformLocation(shader->getProgramIndex(), "fogColor");
	vec_f fog_color({ 0.33f, 0.33f, 0.33f, 0.5f });
    glUniform4fv(loc, 1, fog_color.data());

	obj->transform();

	//loadIdentity(MatrixTypes::MODEL);

	pushMatrix(MatrixTypes::MODEL);

	computeDerivedMatrix(PROJ_VIEW_MODEL);
	glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::VIEW_MODEL]);
	glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::PROJ_VIEW_MODEL]);
	computeNormalMatrix3x3();
	glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

	// loc = glGetUniformLocation(shader->getProgramIndex(), "texMode");
    // glUniform1i(loc, 2);
	glBindVertexArray(obj->get_mesh()->vao);
	glDrawElements(obj->get_mesh()->type, obj->get_mesh()->numIndexes, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	popMatrix(MatrixTypes::MODEL);
}

/*std::vector<std::vector<float>> calculate_corners(Object &obj) {
	float x_max = (obj.get_coords().at(0) + (obj.get_size().at(0) / 2)) * cos(obj.get_y_angle() * PI / 180.0);
	float x_min = (obj.get_coords().at(0) - (obj.get_size().at(0) / 2)) * cos(obj.get_y_angle() * PI / 180.0);
	float z_max = (obj.get_coords().at(2) + (obj.get_size().at(2) / 2)) * cos((obj.get_y_angle()) * PI / 180.0);
	float z_min = (obj.get_coords().at(2) - (obj.get_size().at(2) / 2)) * cos((obj.get_y_angle()) * PI / 180.0);

	std::vector<float> top_left = {x_max, obj.get_coords().at(1), z_min};
	std::vector<float> top_right = {x_max, obj.get_coords().at(1), z_max};
	std::vector<float> bottom_left = {x_min, obj.get_coords().at(1), z_min};
	std::vector<float> bottom_right = {x_min, obj.get_coords().at(1), z_max};

	return {top_left, top_right, bottom_left, bottom_right};
}*/

//Collision car and butter
bool check_collision_AA(Object first, Object second) {
	float x_max1 = (first.get_coords().at(0) + (first.get_size().at(0) / 2)) - (first.get_size().at(2) / 2) * std::abs(sin(first.get_y_angle() * PI / 180.0));
	float x_min1 = (first.get_coords().at(0) - (first.get_size().at(0) / 2));
	float z_max1 = (first.get_coords().at(2) + (first.get_size().at(2) / 2)) + (first.get_size().at(0) / 2) * std::abs(sin(first.get_y_angle() * PI / 180.0));
	float z_min1 = (first.get_coords().at(2) - (first.get_size().at(2) / 2));

	float x_max2 = (second.get_coords().at(0) + (second.get_size().at(0) / 2)) - (second.get_size().at(2) / 2) * std::abs(sin(second.get_y_angle() * PI / 180.0));
	float x_min2 = (second.get_coords().at(0) - (second.get_size().at(0) / 2));
	float z_max2 = (second.get_coords().at(2) + (second.get_size().at(2) / 2)) + (second.get_size().at(0) / 2) * std::abs(sin(second.get_y_angle() * PI / 180.0));
	float z_min2 = (second.get_coords().at(2) - (second.get_size().at(2) / 2));

	return (x_min1 <= x_max2 && x_max1 >= x_min2) && (z_min1 <= z_max2 && z_max1 >= z_min2);
}

//Collision car and cheerios
bool check_collision_CH(Object first, std::vector<float> center, float radius) {
	float x_max1 = (first.get_coords().at(0) + (first.get_size().at(0) / 2)) -
		(first.get_size().at(2) / 2) * std::abs(sin(first.get_y_angle() * PI / 180.0));
	float x_min1 = (first.get_coords().at(0) - (first.get_size().at(0) / 2));
	float z_max1 = (first.get_coords().at(2) + (first.get_size().at(2) / 2)) +
		(first.get_size().at(0) / 2) * std::abs(sin(first.get_y_angle() * PI / 180.0));
	float z_min1 = (first.get_coords().at(2) - (first.get_size().at(2) / 2));

	float x_max2 = (center.at(0) + radius);
	float x_min2 = (center.at(0) - radius);
	float z_max2 = (center.at(2) + radius);
	float z_min2 = (center.at(2) - radius);

	return (x_min1 <= x_max2 && x_max1 >= x_min2) && (z_min1 <= z_max2 && z_max1 >= z_min2);
}

//Collision car and oranges
bool check_collision_SA(Object box, Object sphere, float radius) {
	float x_max = (box.get_coords().at(0) + (box.get_size().at(0) / 2)) -
		(box.get_size().at(2) / 2) * std::abs(sin(box.get_y_angle() * PI / 180.0));
	float x_min = (box.get_coords().at(0) - (box.get_size().at(0) / 2));
	float z_max = (box.get_coords().at(2) + (box.get_size().at(2) / 2)) +
		(box.get_size().at(0) / 2) * std::abs(sin(box.get_y_angle() * PI / 180.0));
	float z_min = (box.get_coords().at(2) - (box.get_size().at(2) / 2));


	float x = fmax(x_min, fmin(sphere.get_coords().at(0), x_max));
	float z = fmax(z_min, fmin(sphere.get_coords().at(2), z_max));

	// this is the same as isPointInsideSphere
	float distance = sqrt((x - sphere.get_coords().at(0)) * (x - sphere.get_coords().at(0))
		+ (z - sphere.get_coords().at(2)) * (z - sphere.get_coords().at(2)));

	return distance < radius;
}

/*void SAT_test(std::vector<float> normal, std::vector<std::vector<float>> corners, float *minAlong, float *maxAlong) {
	*minAlong = 10000, *maxAlong = -10000;
	for (auto &&corner: corners) {
		float dot_val = dotProduct(corner.data(), normal.data());
		if (dot_val < *minAlong) *minAlong = dot_val;
		if (dot_val > *maxAlong) *maxAlong = dot_val;
	}
}

std::vector<std::vector<float>> calculate_normals(std::vector<std::vector<float>> corners) {
	// Calculate normal of every side

	std::vector<float> top_left = corners.at(0);
	std::vector<float> top_right = corners.at(1);
	std::vector<float> bottom_left = corners.at(2);
	std::vector<float> bottom_right = corners.at(3);
	//float y = top_left.at(1);

	std::vector<float> top = {top_right.at(0) - top_left.at(0), 0, top_right.at(2) - top_left.at(2)};
	std::vector<float> bottom = {bottom_right.at(0) - bottom_left.at(0), 0, bottom_right.at(2) - bottom_left.at(2)};
	std::vector<float> left = {top_left.at(0) - bottom_left.at(0), 0, top_left.at(2) - bottom_left.at(2)};
	std::vector<float> right = {top_right.at(0) - bottom_right.at(0), 0, top_right.at(2) - bottom_right.at(2)};

	std::vector<float> n_top = {-top.at(2), 0, top.at(0)};
	std::vector<float> n_bottom = {-bottom.at(2), 0, bottom.at(0)};
	std::vector<float> n_left = {-left.at(2), 0, left.at(0)};
	std::vector<float> n_right = {-right.at(2), 0, right.at(0)};

	return {normalize_vec(n_top), normalize_vec(n_bottom), normalize_vec(n_left), normalize_vec(n_right)};

}



inline bool isBetweenOrdered(float val, float lowerBound, float upperBound) {
	return (lowerBound <= val) && (val <= upperBound);
}

bool overlaps(float min1, float max1, float min2, float max2) {
	return isBetweenOrdered(min2, min1, max1) || isBetweenOrdered(min1, min2, max2);
}

bool check_intersection(Object first, Object second, std::vector<std::vector<float>> normals_1,
						std::vector<std::vector<float>> normals_2, std::vector<std::vector<float>> corners_1,
						std::vector<std::vector<float>> corners_2) {
	for (auto &&normal: normals_1) {
		float shape1_min, shape1_max, shape2_min, shape2_max;
		SAT_test(normal, corners_1, &shape1_min, &shape1_max);
		SAT_test(normal, corners_2, &shape2_min, &shape2_max);
		if (!overlaps(shape1_min, shape1_max, shape2_min, shape2_max)) {
			return false;
		}
	}

	for (auto &&normal: normals_2) {
		float shape1_min, shape1_max, shape2_min, shape2_max;
		SAT_test(normal, corners_1, &shape1_min, &shape1_max);
		SAT_test(normal, corners_2, &shape2_min, &shape2_max);
		if (!overlaps(shape1_min, shape1_max, shape2_min, shape2_max)) {
			return false;
		}
	}

	return true;
}*/

/*bool check_collision(Object first, Object second) {
	auto first_corners = calculate_corners(first);
	auto second_corners = calculate_corners(second);
	auto first_normals = calculate_normals(first_corners);
	auto second_normals = calculate_normals(second_corners);

//    for (auto &comp : first_corners) {
//        std::cout << comp;
//    }

	return check_intersection(first, second, first_normals, second_normals, first_corners, second_corners);
}

bool CheckCollision(Object one, Object two) { // AABB - AABB
	std::cout << "car: " << one.get_coords().at(0) + (one.get_size().at(0) / 2) << std::endl;
	std::cout << "butter_bottom: " << two.get_coords().at(0) - (two.get_size().at(0) / 2) << std::endl;
	std::cout << "butter_bottom: " << two.get_coords().at(0) + (two.get_size().at(0) / 2) << std::endl;

	// collision x-axis?
	bool collisionX = ((two.get_coords().at(0) - (two.get_size().at(0) / 2)) <=
					   (one.get_coords().at(0) + (one.get_size().at(0) / 2))) &&
					  ((one.get_coords().at(0) + (one.get_size().at(0) / 2)) <=
					   (two.get_coords().at(0) + (two.get_size().at(0) / 2)));

	// collision y-axis?
	//bool collisionY = one.get_coords().at(1) - (one.get_size().at(1)/ 2) >= two.get_coords().at(1) + (two.get_size().at(1) /2) &&
	two.get_coords().at(1) - (two.get_size().at(1) / 2) >= one.get_coords().at(1) + (one.get_size().at(1) / 2);

	//collision z-axis?
	bool collisionZ = ((two.get_coords().at(2) - (two.get_size().at(2) / 2)) <=
					   (one.get_coords().at(2) + (one.get_size().at(2) / 2))) &&
					  ((one.get_coords().at(2) + (one.get_size().at(2) / 2)) <=
					   (two.get_coords().at(2) + (two.get_size().at(2) / 2)));

	return collisionX && collisionZ;
}*/

void reset_cam() {
	Car* car = dynamic_cast<Car*>(scene.get_entities()->at(EntityType::CAR).get());
	alphaAux = -car->get_angle() - 90;

	_alpha = alphaAux;
}




void timer(int value)
{
	std::ostringstream oss;
	oss << CAPTION << " (" << FrameCount << " FPS) : [" << WIDTH << "x" << HEIGHT << "]";
	std::string s = oss.str();
	glutSetWindow(WindowHandle);
	glutSetWindowTitle(s.c_str());
    FrameCount = 0;
    glutTimerFunc(1000, timer, 0);
}

void refresh(int value)
{
	//PUT YOUR CODE HERE
}

// ------------------------------------------------------------
//
// Reshape Callback Function
//

void changeSize(int w, int h) {

	float ratio;
	// Prevent a divide by zero, when window is too short
	if(h == 0)
		h = 1;
	// set the viewport to be the entire window
	glViewport(0, 0, w, h);
	// set the projection matrix
	ratio = (1.0f * w) / h;
	loadIdentity(PROJECTION);
	perspective(53.13f, ratio, 0.1f, 1000.0f);
}


void render_text() {
	glDisable(GL_DEPTH_TEST);
	//the glyph contains background colors and non-transparent for the actual character pixels. So we use the blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	int m_viewport[4];
	glGetIntegerv(GL_VIEWPORT, m_viewport);
	//viewer looking down at  negative z direction
	pushMatrix(MODEL);
	loadIdentity(MODEL);
	pushMatrix(PROJECTION);
	loadIdentity(PROJECTION);
	pushMatrix(VIEW);
	loadIdentity(VIEW);
	ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);

	auto _life = str_life;
	for (int n = 1; n <= lifes_left; n++) {
		_life += " <3";
	}

	if (points_n < 0.0f) {
		points_n = 0.0f;
	}

	// std::cout << points_n << std::endl;
	auto _points = str_points + " " + std::to_string(int(ceil(points_n)));

	if (_pause) {
		auto border = vec_f({0.20f, 0.107f, 0.58f});
		auto frground = vec_f({0.234f, 0.70f, 0.48});
		RenderText(shaderText, "Game is paused!", 125.0f, WinY / 2 + 150, 3.00f, border.at(0), border.at(1), border.at(2));
		RenderText(shaderText, "Finish your soup later!", 446.0f, WinY / 2 - 20, 1.0f, border.at(0), border.at(1), border.at(2));
		RenderText(shaderText, "Press 's' and get back to stacking cheese", 247.0f, WinY / 2 - 70, 1.0f, border.at(0), border.at(1), border.at(2));

		RenderText(shaderText, "Game is paused!", 130.0f, WinY / 2 + 150, 3.0f, frground.at(0), frground.at(1), frground.at(2));
		RenderText(shaderText, "Finish your soup later!", 450.0f, WinY / 2 - 20, 1.0f, frground.at(0), frground.at(1), frground.at(2));
		RenderText(shaderText, "Press 's' and get back to stacking cheese", 250.0f, WinY / 2 - 70, 1.0f, frground.at(0), frground.at(1), frground.at(2));
	} else if (_game_over) {
		RenderText(shaderText, "GAME OVER!", 225.0f, WinY / 2 + 250.0f, 3.03f, 0.71f, 0.71f, 0.71f);
		RenderText(shaderText, "Maybe next time :P", 105.0f, WinY / 2 - 50.0f, 3.03f, 0.71f, 0.71f, 0.71f);

		RenderText(shaderText, "GAME OVER!", 230.0f, WinY / 2 + 250.0f, 3.0f, 0.71f, 0.24f, 0.14f);
		RenderText(shaderText, "Maybe next time :P", 110.0f, WinY / 2 - 50.0f, 3.0f, 0.71f, 0.24f, 0.14f);
	} else {
		RenderText(shaderText, _life, 15.0f, 25.0f, 1.0f, 0.5f, 0.8f, 0.2f);
		RenderText(shaderText, _points, 15.0f, 720.0f, 1.0f, 0.3f, 0.7f, 0.9f);
	}

	popMatrix(PROJECTION);
	popMatrix(VIEW);
	popMatrix(MODEL);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
}


// ------------------------------------------------------------
//
// Render stufff
//

void render_flare(FLARE_DEF *flare, int lx, int ly, int *m_viewport) {  //lx, ly represent the projected position of light on viewport

	int     dx, dy;          // Screen coordinates of "destination"
	int     px, py;          // Screen coordinates of flare element
	int		cx, cy;
	float    maxflaredist, flaredist, flaremaxsize, flarescale, scaleDistance;
	int     width, height, alpha;    // Piece parameters;
	int     i;
	float	diffuse[4];

	GLint loc;

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	int screenMaxCoordX = m_viewport[0] + m_viewport[2] - 1;
	int screenMaxCoordY = m_viewport[1] + m_viewport[3] - 1;

	//viewport center
	cx = m_viewport[0] + (int)(0.5f * (float)m_viewport[2]) - 1;
	cy = m_viewport[1] + (int)(0.5f * (float)m_viewport[3]) - 1;

	// Compute how far off-center the flare source is.
	maxflaredist = sqrt(cx*cx + cy * cy);
	flaredist = sqrt((lx - cx)*(lx - cx) + (ly - cy)*(ly - cy));
	scaleDistance = (maxflaredist - flaredist) / maxflaredist;
	flaremaxsize = (int)(m_viewport[2] * flare->fMaxSize);
	flarescale = (int)(m_viewport[2] * flare->fScale);

	

	// Destination is opposite side of centre from source
	dx = clampi(cx + (cx - lx), m_viewport[0], screenMaxCoordX);
	dy = clampi(cy + (cy - ly), m_viewport[1], screenMaxCoordY);

	// Render each element. To be used Texture Unit 0

	glUniform1i(texMode_uniformId, 3); // draw modulated textured particles 
	glUniform1i(tex_loc, 0);  //use TU 0

	for (i = 0; i < flare->nPieces; ++i)
	{
		
		// Position is interpolated along line between start and destination.
		px = (int)((1.0f - flare->element[i].fDistance)*lx + flare->element[i].fDistance*dx);
		py = (int)((1.0f - flare->element[i].fDistance)*ly + flare->element[i].fDistance*dy);
		px = clampi(px, m_viewport[0], screenMaxCoordX);
		py = clampi(py, m_viewport[1], screenMaxCoordY);

		// Piece size are 0 to 1; flare size is proportion of screen width; scale by flaredist/maxflaredist.
		width = (int)(scaleDistance*flarescale*flare->element[i].fSize);
		
		// Width gets clamped, to allows the off-axis flaresto keep a good size without letting the elements get big when centered.
		if (width > flaremaxsize)  width = flaremaxsize;
		
		height = (int)((float)m_viewport[3] / (float)m_viewport[2] * (float)width);
		memcpy(diffuse, flare->element[i].matDiffuse, 4 * sizeof(float));
		diffuse[3] *= scaleDistance;   //scale the alpha channel
		
		int objId = 0;
		if (width > 1)
		{

			// send the material - diffuse color modulated with texture
			loc = glGetUniformLocation(shader.getProgramIndex(), "mat.diffuse");
			glUniform4fv(loc, 1, diffuse);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, FlareTextureArray[flare->element[i].textureId]);

			pushMatrix(MODEL);
			translate(MODEL, (float)(px - width * 0.0f), (float)(py - height * 0.0f), 0.0f);
			scale(MODEL, (float)width, (float)height, 1);
			computeDerivedMatrix(PROJ_VIEW_MODEL);
			glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[VIEW_MODEL]);
			glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[PROJ_VIEW_MODEL]);
			computeNormalMatrix3x3();
			glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

			glBindVertexArray(myMeshes[0].vao);
			glDrawElements(myMeshes[0].type, myMeshes[0].numIndexes, GL_UNSIGNED_INT, 0);
			glBindVertexArray(0);
			popMatrix(MODEL);
		}
	}
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
}


void renderScene() {
	if (_pause) {
		t_delta = 0;
	}
	if (_first_time) {
		roadSetUp();
		_first_time = false;
	}
	if (_reset) {
		roadSetUp();
		_reset = false;
	}
	if (_to_reset) {
		_to_reset = false;
		scene.reset();
	}

	std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

	FrameCount++;
	/*glClearColor(0.2f, 0.2f, 0.2f, 1.0f);*/
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	lookAt(camX, camY, camZ, 0,0,0, 0,1,0);
	// load identity matrices
	loadIdentity(VIEW);
	loadIdentity(MODEL);

	scene.cam_update();
	// use our shader
	glUseProgram(shader.getProgramIndex());

	if (!_pause) {
		scene.input();
		scene.update();
	}

	float res[4];
	float mat[16];
	GLfloat plane_floor[4] = { 0,1,0,0 };

	glEnable(GL_DEPTH_TEST);

	Object table = scene.get_entities()->at(EntityType::TABLE)->get_objects()->at(0);

	if(camY > 0.0f) {  //camera in front of the floor so render reflections and shadows. Inner product between the viewing direction and the normal of the ground
	
		glEnable(GL_STENCIL_TEST);
		glStencilFunc(GL_NEVER, 0x1, 0x1);
		glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);

		
		draw_mirror(&shader, &table);

		glStencilFunc(GL_EQUAL, 0x1, 0x1);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

		lightPosG[1] *= (-1.0f);  //mirror the position of light
		multMatrixPoint(VIEW, lightPosG, res);

		glUniform4fv(lPosG_uniformId, 1, res);
		pushMatrix(MatrixTypes::MODEL);
		scale(MatrixTypes::MODEL, 1.0f, -1.0f, 1.0f);
		glCullFace(GL_FRONT);
		scene.draw(&shader);
		glCullFace(GL_BACK);
		popMatrix(MatrixTypes::MODEL);

		lightPosG[1] *= (-1.0f);  //reset the light position
		multMatrixPoint(MatrixTypes::VIEW, lightPosG, res);
		glUniform4fv(lPosG_uniformId, 1, res);

		
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);// Blend specular Ground with reflected geometry
		draw_mirror(&shader, &table);
	
		//Render the Shadows
		 glUniform1i(shadowMode_uniformId, 1);  //Render with constant color
		 shadow_matrix(mat, plane_floor, lightPos);

		glDisable(GL_DEPTH_TEST);

		glBlendFunc(GL_DST_COLOR, GL_ZERO);
		glStencilOp(GL_KEEP, GL_KEEP, GL_ZERO);

		pushMatrix(MatrixTypes::MODEL);
		multMatrix(MatrixTypes::MODEL, mat);

		scene.draw(&shader);
		popMatrix(MatrixTypes::MODEL);

		glDisable(GL_STENCIL_TEST);
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);

		//render the geometry
		glUniform1i(shadowMode_uniformId, 0);

		scene.draw(&shader);

	} else{

		glUniform1i(shadowMode_uniformId, 0);
		scene.draw(&shader);
		draw_mirror(&shader, &table);
		//scene.draw(&shader);
	}
	


	//send the light position in eye coordinates

		//glUniform4fv(lPos_uniformId, 1, lightPos); //efeito capacete do mineiro, ou seja lighPos foi definido em eye coord 

		//float res[4];
		//multMatrixPoint(VIEW, lightPos,res);   //lightPos definido em World Coord so is converted to eye space
		//glUniform4fv(lPos_uniformId, 1, res);



	glutSwapBuffers();

	std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
	t_delta = time_span.count();
}

// ------------------------------------------------------------
//
// Events from the Keyboard
//

void KeyboardUpHandler(unsigned char key, int xx, int yy) {
	if (_pause) {
		return;
	}
	switch (key) {
		case 'q':
		case GLUT_KEY_UP:
			car_forwards = false;
			break;
		case 'a':
		case GLUT_KEY_DOWN:
			car_backwards = false;
			break;
		case 'o':
		case GLUT_KEY_LEFT:
			car_left = false;
			break;
		case 'p':
		case GLUT_KEY_RIGHT:
			car_right = false;
			break;
	}
}

void special_keys_callback(int key, int xx, int yy) {
	switch (key) {
		case (GLUT_KEY_F5):
			type++;
			if (type == 5) type = 0;
			break;

		case (GLUT_KEY_UP):
			car_forwards = true;
			break;
		case (GLUT_KEY_DOWN):
			car_backwards = true;
			break;
		case (GLUT_KEY_LEFT):
			car_left = true;
			break;
		case (GLUT_KEY_RIGHT):
			car_right = true;
			break;
		default:
			break;
	}
}
void special_keys_up_callback(int key, int xx, int yy) {
	switch (key) {
		case (GLUT_KEY_UP):
			car_forwards = false;
			break;
		case (GLUT_KEY_DOWN):
			car_backwards = false;
			break;
		case (GLUT_KEY_LEFT):
			car_left = false;
			break;
		case (GLUT_KEY_RIGHT):
			car_right = false;
			break;
		default:
			break;
	}
}


void processKeys(unsigned char key, int xx, int yy) {
	if (_pause) {
		if (key == 's') {
			_pause = false;
		}
		return;
	}
	switch(key) {
		case 'q':
			car_forwards = true;
			break;
		case 'a':
			car_backwards = true;
			break;
		case 'o':
			car_left = true;
			break;
		case 'p':
			car_right = true;
			break;
		case 27:
			glutLeaveMainLoop();
			break;

		default:
			scene.key_input(key, xx, yy);
	}
}


// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	if (_pause) {
		return;
	}
	// start tracking the mouse
	if (state == GLUT_DOWN)  {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	//stop tracking the mouse
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			_alpha -= (xx - startX);
			_beta += (yy - startY);
		}
		else if (tracking == 2) {
			r += (yy - startY) * 0.01f;
			if (r < 0.1f)
				r = 0.1f;
		}
		tracking = 0;
	}

	// TODO maybe able to delete this. When deleted, camera orbits faster
	camX = r * sin(_alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camZ = r * cos(_alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camY = r * sin(_beta * 3.14f / 180.0f);
}

// Track mouse motion while buttons are pressed

void processMouseMotion(int xx, int yy)
{
	if (_pause) {
		return;
	}

	int deltaX, deltaY;

	deltaX =  - xx + startX;
	deltaY =    yy - startY;

	Car* car = dynamic_cast<Car*>(scene.get_entities()->at(EntityType::CAR).get());
	alphaAux = -car->get_angle() - 90;

	_beta = _beta < -90 ? -89 : _beta;
	_beta = _beta > 90 ? 89 : _beta;

	// left mouse button: move camera
	if (tracking == 1) {


		alphaAux = _alpha + deltaX;
		betaAux = _beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;
		rAux = r;
	}
	// right mouse button: zoom
	else if (tracking == 2) {

		alphaAux = _alpha;
		betaAux = _beta;
		rAux = r + (deltaY * 0.01f);
		if (rAux < 0.1f)
			rAux = 0.1f;
	}

	betaAux = betaAux < -90 ? -89 : betaAux;
	betaAux = betaAux > 90 ? 89 : betaAux;

	betaAux = betaAux == -90 ? -89 : betaAux;
	betaAux = betaAux == 90 ? 89 : betaAux;

//  uncomment this if not using an idle or refresh func
//	glutPostRedisplay();
}


void mouseWheel(int wheel, int direction, int x, int y) {
	if (_pause) {
		return;
	}

	r -= direction * 0.1f;
	if (r < 0.1f)
		r = 0.1f;

	rAux = r;

//  uncomment this if not using an idle or refresh func
//	glutPostRedisplay();
}

// --------------------------------------------------------
//
// Shader Stuff
//


GLuint setupShaders() {

	// Shader for models
	shader.init();
	shader.loadShader(VSShaderLib::VERTEX_SHADER, VERT);
	shader.loadShader(VSShaderLib::FRAGMENT_SHADER, FRAG);

	// set semantics for the shader variables
	glBindFragDataLocation(shader.getProgramIndex(), 0,"colorOut");
	glBindAttribLocation(shader.getProgramIndex(), VERTEX_COORD_ATTRIB, "position");
	glBindAttribLocation(shader.getProgramIndex(), NORMAL_ATTRIB, "normal");
	glBindAttribLocation(shader.getProgramIndex(), TEXTURE_COORD_ATTRIB, "texCoord");

	glLinkProgram(shader.getProgramIndex());

	texMode_uniformId = glGetUniformLocation(shader.getProgramIndex(), "texMode"); // different modes of texturing
	pvm_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_pvm");
	vm_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_viewModel");
	model_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_Model");
	normal_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_normal");
	view_uniformId = glGetUniformLocation(shader.getProgramIndex(), "m_View");
	shadowMode_uniformId = glGetUniformLocation(shader.getProgramIndex(), "shadowMode");
	lPos_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_pos");
	lPos2_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_pos2");
	lPos3_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_pos3");
	lPos4_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_pos4");
	lPos5_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_pos5");
	lPos6_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_pos6");
	lPosG_uniformId = glGetUniformLocation(shader.getProgramIndex(), "l_posG");
	
	printf("InfoLog for Per Fragment Phong Lightning Shader\n%s\n\n", shader.getAllInfoLogs().c_str());
	
	shaderText.init();
	shaderText.loadShader(VSShaderLib::VERTEX_SHADER, "shaders/text.vert");
	shaderText.loadShader(VSShaderLib::FRAGMENT_SHADER, "shaders/text.frag");
	glLinkProgram(shaderText.getProgramIndex());
	printf("InfoLog for Text Rendering Shader\n%s\n\n", shaderText.getAllInfoLogs().c_str());
	return(shader.isProgramLinked() && shaderText.isProgramLinked());
}

// ------------------------------------------------------------
//
// Model loading and OpenGL setup
//

void init()
{
	// set the camera position based on its spherical coordinates
	camX = r * sin(_alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camZ = r * cos(_alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camY = r * sin(_beta * 3.14f / 180.0f);

	MyMesh amesh;
	/* Initialization of DevIL */
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	freeType_init(font_name);

	glGenTextures(3, TextureArray);
	Texture2D_Loader(TextureArray, ASSET1, 0);
	Texture2D_Loader(TextureArray, ASSET2, 1);
	Texture2D_Loader(TextureArray, ASSET3, 2);

	//Flare elements textures
	glGenTextures(5, FlareTextureArray);
	Texture2D_Loader(FlareTextureArray, "crcl.tga", 0);
	Texture2D_Loader(FlareTextureArray, "flar.tga", 1);
	Texture2D_Loader(FlareTextureArray, "hxgn.tga", 2);
	Texture2D_Loader(FlareTextureArray, "ring.tga", 3);
	Texture2D_Loader(FlareTextureArray, "sun.tga", 4);

	scene.create_objs();
	Car* car = dynamic_cast<Car*>(scene.get_entities()->at(EntityType::CAR).get());
	alphaAux = -car->get_angle() - 90;

	//Load flare from file
	amesh = createQuad(1, 1);
	myMeshes.push_back(amesh);

	//Load flare from file
	loadFlareFile(&AVTflare, "flare.txt");


	// some GL settings
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
	//glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearColor(255.0f, 255.0f, 255.0f, 1.0f);

}


// ------------------------------------------------------------
//
// Main function
//


int main(int argc, char **argv) {
	lifes_left = STARTING_LIFES;

//  GLUT initialization
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA|GLUT_STENCIL|GLUT_MULTISAMPLE);

	glutInitContextVersion (3, 3);
	glutInitContextProfile (GLUT_CORE_PROFILE );
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE | GLUT_DEBUG);

	glutInitWindowPosition(100,100);
	glutInitWindowSize(WIDTH, HEIGHT);
	WindowHandle = glutCreateWindow(CAPTION);

	//glutSetKeyRepeat(GLUT_KEY_REPEAT_OFF);


//  Callback Registration
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);

	glutTimerFunc(0, timer, 0);
	glutIdleFunc(renderScene);  // Use it for maximum performance
	//glutTimerFunc(0, refresh, 0);    //use it to to get 60 FPS whatever

//	Mouse and Keyboard Callbacks
	glutKeyboardFunc(processKeys);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutKeyboardUpFunc(KeyboardUpHandler);
	glutSpecialFunc(special_keys_callback);
	glutSpecialUpFunc(special_keys_up_callback);
	glutMouseWheelFunc(mouseWheel);

//	return from main loop
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

//	Init GLEW
	glewExperimental = GL_TRUE;
	GLenum res = glewInit();
	if (res != GLEW_OK) {
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		exit(EXIT_FAILURE);
	}

	printf ("Vendor: %s\n", glGetString (GL_VENDOR));
	printf ("Renderer: %s\n", glGetString (GL_RENDERER));
	printf ("Version: %s\n", glGetString (GL_VERSION));
	printf ("GLSL: %s\n", glGetString (GL_SHADING_LANGUAGE_VERSION));

	if (!setupShaders())
		return(1);

	init();

	//  GLUT main loop
	glutMainLoop();

	return(0);

}

unsigned int getTextureId(char *name) {
	int i;

	for (i = 0; i < NTEXTURES; ++i)
	{
		if (strncmp(name, flareTextureNames[i], strlen(name)) == 0)
			return i;
	}
	return -1;
}


void    loadFlareFile(FLARE_DEF *flare, char *filename)
{
	int     n = 0;
	FILE* f;
	char    buf[256];
	int fields;

	memset(flare, 0, sizeof(FLARE_DEF));

	f = fopen(filename, "r");
	if (f)
	{
		fgets(buf, sizeof(buf), f);
		sscanf(buf, "%f %f", &flare->fScale, &flare->fMaxSize);

		while (!feof(f))
		{

			printf("Flare file opening \n");
			char            name[8] = { '\0', };
			double          dDist = 0.0, dSize = 0.0;
			float			color[4];
			int				id;

			fgets(buf, sizeof(buf), f);
			fields = sscanf(buf, "%8s %lf %lf ( %f %f %f %f )", name, &dDist, &dSize, &color[0], &color[1], &color[2], &color[3]);
			if (fields == 7)
			{
				for (int i = 0; i < 4; ++i) color[i] = clamp(color[i] / 255.0f, 0.0f, 1.0f);
				id = getTextureId(name);
				if (id < 0) printf("Texture name not recognized\n");
				else
					flare->element[n].textureId = id;
				flare->element[n].fDistance = (float)dDist;
				flare->element[n].fSize = (float)dSize;
				memcpy(flare->element[n].matDiffuse, color, 4 * sizeof(float));
				++n;
			}
		}

		flare->nPieces = n;
		fclose(f);
	}
	else printf("Flare file opening error\n");
}

