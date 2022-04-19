#pragma once

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <vector>
#include <memory>
#include <cstdint>
#include <utility>
#include <cstring>
#include <chrono>
#include <ctime>
#include <ratio>
#include <ranges>
#include <cmath>
#include <chrono>
#include <random>
#include <iostream>
#include <limits>
#include <cctype>
#include "vsShaderLib.h"
#include "AVTmathLib.h"
#include "geometry.h"
#include "Texture_Loader.h"
#include "flare.h"
#include "l3dBillboard.h"


#define PI 3.14159265
#define DISTANCE_TO_CAR 5.0

extern bool _reset;

extern bool oranges_kill;

extern inline double clamp(const double x, const double min, const double max);

extern inline int clampi(const int x, const int min, const int max);

extern GLint pvm_uniformId;
extern GLint vm_uniformId;
extern GLint normal_uniformId;
extern GLint view_uniformId;
extern GLint lPos_uniformId;
extern GLint lPos2_uniformId;
extern GLint lPos3_uniformId;
extern GLint lPos4_uniformId;
extern GLint lPos5_uniformId;
extern GLint lPos6_uniformId;
extern GLint lPosG_uniformId;

extern GLint texMode_uniformId;

extern float lightScreenPos[3]; 
extern GLuint TextureArray[3];
extern GLuint FlareTextureArray[5];
extern FLARE_DEF AVTflare;


/// The storage for matrices
extern float mMatrix[COUNT_MATRICES][16];
extern float mCompMatrix[COUNT_COMPUTED_MATRICES][16];

/// The normal matrix
extern float mNormal3x3[9];

extern struct OrthoProperties orthoProperties;
extern float _alphaStatic, _betaStatic;
extern float rStatic;
extern float alphaAux, betaAux;
extern float rAux;
typedef enum {
    TOP_ORTHO, TOP_PERS, FOLLOW
} Camera;
extern float _alpha, _beta;

extern float camX, camY, camZ;

extern const int WIDTH;
extern const int HEIGHT;
extern float r;

extern double t_delta;

typedef std::vector<float> vec_f;

//typedef std::chrono::high_resolution_clock clock;

using u32 = uint_least32_t;
using engine = std::mt19937;

extern int amount; // Cheerios
extern int max_particles; //snow
extern bool snow;
int dead_num_particles = 0;
int type;

bool bumpmap;

extern float lightPos[4];

extern float lightPos2[4];
extern float lightPos3[4];
extern float lightPos4[4];
extern float lightPos5[4];
extern float lightPos6[4];
extern float lightPosG[4];

extern bool spotlight_mode;
extern bool directional;
extern bool flareEffect;
extern std::vector<std::vector<float>> cheercoor;

extern int lifes_left;
extern double points_n;

extern bool _game_over;

extern bool _pause;

extern bool _to_reset;


struct Particle {
    float	life;		// life
    float	fade;		// fade
    std::vector<float> position;
    float r, g , b;   //color
    GLfloat vx, vy, vz; // velocity
    GLfloat ax, ay, az; // accelaration
};

extern std::vector<struct Particle> particles;

class Entity;
class Object;
class Orange;

void roadSetUp();
bool check_collision_AA(Object first, Object second);
bool check_collision_SA(Object box, Object sphere, float radius);
bool check_collision_CH(Object first, std::vector<float> center, float radius);

void render_text();
void iniParticles();
void updateParticles();

void reset_cam();
void changeSize(int w, int h);
void render_skybox(VSShaderLib* shader, Object* obj);


extern std::chrono::high_resolution_clock::time_point part_timer;
extern std::chrono::high_resolution_clock::time_point part_timer_end;

inline float rand_f(float top) {
    std::random_device os_seed;
    const u32 seed = os_seed();

    engine generator(seed);
    std::uniform_int_distribution<u32> distribute(1.0, top);
    float rand_f = distribute(generator) / 1.1;
    rand_f = distribute(generator) > top / 2 ? -rand_f : rand_f;
    return rand_f;
}

inline float rand_uf(float top) {
    std::random_device os_seed;
    const u32 seed = os_seed();

    engine generator(seed);
    std::uniform_int_distribution<u32> distribute(1.0, top);
    float rand_f = distribute(generator) / 1.1;
    return rand_f;
}


/// Objects index
enum EntityType {
    TABLE,
    ORANGES,
    BUTTER,
    //SKYBOX,
    TREES,
    SNOW,
    ROAD,
    CAR,
};

enum TransType {
    Scale,
    Translate,
    RotateX,
    RotateY,
    RotateZ
};

enum Direction {
    FORWARD,
    BACKWARDS
};

extern bool car_forwards;
extern bool car_backwards;
extern bool car_left;
extern bool car_right;


class Transformation {
    TransType _type;
    std::vector<float> _comps;
    /// Default is that the transformation is to be modified, meaning is the one to be used in the dynamic transformation
    bool _modifiable;

public:
    Transformation(TransType type, std::vector<float> comps, bool modifiable = true) : _type(type), _comps(comps),
        _modifiable(modifiable) {}

    void set(std::vector<float> comps) {
        _comps = comps;
    }

    TransType get_type() {
        return _type;
    }

    std::vector<float>* get_comps() {
        return &_comps;
    }

    bool is_modifiable() {
        return _modifiable;
    }

    void apply() {
        switch (_type) {
        case TransType::Scale:
            scale(MatrixTypes::MODEL,
                _comps.at(0),
                _comps.at(1),
                _comps.at(2));
            break;
        case TransType::Translate:
            translate(MatrixTypes::MODEL,
                _comps.at(0),
                _comps.at(1),
                _comps.at(2));
            break;
        case TransType::RotateX:
        case TransType::RotateY:
        case TransType::RotateZ:
            rotate(MatrixTypes::MODEL,
                _comps.at(0),
                _comps.at(1),
                _comps.at(2),
                _comps.at(3));
            break;
        default:
            throw std::runtime_error("Transformation type provided not supported!!");
            break;
        }
    }
};

typedef std::vector<Transformation> vec_t;

class Object {
protected:
    MyMesh _mesh;

    std::vector<Transformation> _transfs;

    std::vector<float> _initial_offset;
    std::vector<float> _offset;

    std::vector<float> _size;

    std::vector<GLuint> _textureArray;

    float _y_angle = 0;

    bool _drawable = true;

public:
    Object() {}

    Object(MyMesh mesh, std::vector<Transformation> transfs,
        std::vector<float> initial_offset = vec_f({ 0.0f, 0.0f, 0.0f }),
        std::vector<float> size = vec_f({ 0.0f, 0.0f, 0.0f }), std::vector<GLuint> textureArray = {}) :
        _mesh(mesh), _transfs(transfs), _initial_offset(initial_offset), _size(size), _textureArray(textureArray) {

        _offset = { 0.0f, 0.0f, 0.0f };
    }

    std::vector<Transformation>* get_transfs() {
        return &this->_transfs;
    }

    void set_transf_at(size_t ix, std::vector<float> transf) {
        _transfs.at(ix).set(transf);
    }

    std::vector<float>* get_offset() {
        return &this->_offset;
    }

    std::vector<float> get_ini_offset() {
        return _initial_offset;
    }

    MyMesh* get_mesh(){
        return &this->_mesh;
    }

    void set_y_angle(float y_angle) {
        _y_angle = y_angle;
    }

    float get_y_angle() {
        return _y_angle;
    }

    void set_mesh_diffuse_alpha(float alpha) {
        this->_mesh.mat.diffuse[3] = alpha;
    }

    std::vector<GLuint>* get_textureArray() {
        return &this->_textureArray;
    }

    bool is_drawable() {
        return _drawable;
    }

    void set_drawability(bool drawable) {
        _drawable = drawable;
    }


    std::vector<float> get_size() {
        return this->_size;
    }


    std::vector<float> get_coords() {
        vec_f coords(3);
        size_t comp = 0;
        for (auto&& coord : coords) {
            coord = _initial_offset.at(comp) + _offset.at(comp);

            comp++;
        }

        return coords;
    }

    void set_ini_offset(std::vector<float> ini_offset) {
        _initial_offset = ini_offset;
    }

    void set_offset(std::vector<float> offset) {
        _offset = offset;
    }

    void draw(VSShaderLib* shader, size_t entity, double _time_span) {
        GLint loc;

        float res[4];
        multMatrixPoint(VIEW, lightPos, res);
        glUniform4fv(lPos_uniformId, 1, lightPos);
        multMatrixPoint(VIEW, lightPos2, res);
        glUniform4fv(lPos2_uniformId, 1, lightPos2);
        multMatrixPoint(VIEW, lightPos3, res);
        glUniform4fv(lPos3_uniformId, 1, lightPos3);
        multMatrixPoint(VIEW, lightPos4, res);
        glUniform4fv(lPos4_uniformId, 1, lightPos4);
        multMatrixPoint(VIEW, lightPos5, res);
        glUniform4fv(lPos5_uniformId, 1, lightPos5);
        multMatrixPoint(VIEW, lightPos6, res);
        glUniform4fv(lPos6_uniformId, 1, lightPos6);
        multMatrixPoint(VIEW, lightPosG, res);
        glUniform4fv(lPosG_uniformId, 1, lightPosG);

        //if(snow) glDepthMask(GL_FALSE);

        loc = glGetUniformLocation(shader->getProgramIndex(), "sky_flag");
        glUniform1i(loc, false);

        
        if (spotlight_mode) {

            loc = glGetUniformLocation(shader->getProgramIndex(), "spotlight_mode");
            glUniform1i(loc, 1);
        }
        else {
            loc = glGetUniformLocation(shader->getProgramIndex(), "spotlight_mode");
            glUniform1i(loc, 0);
        }
        if (directional) {

            loc = glGetUniformLocation(shader->getProgramIndex(), "directional");
            glUniform1i(loc, 1);
        }
        else {
            loc = glGetUniformLocation(shader->getProgramIndex(), "directional");
            glUniform1i(loc, 0);
        }


        loc = glGetUniformLocation(shader->getProgramIndex(), "coneDir");

        float coneDir[4] = { 0.0f, 0.0f, -1.0f, 0.0f };  //already in eye coordinates
        glUniform4fv(loc, 1, coneDir);
        loc = glGetUniformLocation(shader->getProgramIndex(), "spotCosCutOff");
        glUniform1f(loc, 0.93f);

        if (entity != EntityType::TREES) {
        //loc = glGetUniformLocation(shader->getProgramIndex(), "bill_flag");
        //glUniform1i(loc, false);
        
            if (!_textureArray.empty()) {
                loc = glGetUniformLocation(shader->getProgramIndex(), "text_flag");
                glUniform1i(loc, true);
                
                if (_textureArray.size() == 2){
                    if(entity!=EntityType::BUTTER) {
                        loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
                        glUniform1i(loc, 1);
                    }else if (entity == EntityType::BUTTER){
                        if (bumpmap){
                            loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
                            glUniform1i(loc, 5);
                           
                            loc = glGetUniformLocation(shader->getProgramIndex(), "flag_bumpmap");
                            glUniform1i(loc, true);                    
                        }else{
                            loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
                            glUniform1i(loc, 6);
                            glUniform1i(texMode_uniformId, 6);
                            loc = glGetUniformLocation(shader->getProgramIndex(), "flag_bumpmap");
                            glUniform1i(loc, false);
                        }
                    }
                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GL_TEXTURE_2D, this->_textureArray[0]);
                    glActiveTexture(GL_TEXTURE1);
                    glBindTexture(GL_TEXTURE_2D, this->_textureArray[1]);

                    loc = glGetUniformLocation(shader->getProgramIndex(), "texMode");
                    glUniform1i(loc, 1);
                    loc = glGetUniformLocation(shader->getProgramIndex(), "texmap0");
                    glUniform1i(loc, 0);
                    loc = glGetUniformLocation(shader->getProgramIndex(), "texmap1");
                    glUniform1i(loc, 1);

                }else if (_textureArray.size() == 1 && entity==EntityType::SNOW) {
                    loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
                    glUniform1i(loc, 2);
                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GL_TEXTURE_2D, this->_textureArray[0]);

                    loc = glGetUniformLocation(shader->getProgramIndex(), "texMode");
                    glUniform1i(loc, 1);
                    loc = glGetUniformLocation(shader->getProgramIndex(), "texmap0");
                    glUniform1i(loc, 0);
                }

            }
            else {
                loc = glGetUniformLocation(shader->getProgramIndex(), "text_flag");
                glUniform1f(loc, false);
                loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
                glUniform1i(loc, 0);
                //loc = glGetUniformLocation(shader->getProgramIndex(), "billboard_flag");
                //glUniform1i(loc, false);
            }

            loc = glGetUniformLocation(shader->getProgramIndex(), "mat.diffuse");
            glUniform4fv(loc, 1, _mesh.mat.diffuse);
            loc = glGetUniformLocation(shader->getProgramIndex(), "mat.ambient");
            glUniform4fv(loc, 1, _mesh.mat.ambient);
            loc = glGetUniformLocation(shader->getProgramIndex(), "mat.specular");
            glUniform4fv(loc, 1, _mesh.mat.specular);
            loc = glGetUniformLocation(shader->getProgramIndex(), "mat.shininess");
            glUniform1f(loc, _mesh.mat.shininess);

            loc = glGetUniformLocation(shader->getProgramIndex(), "fogColor");
            vec_f fog_color({ 0.33f, 0.33f, 0.33f, 0.5f });
            glUniform4fv(loc, 1, fog_color.data());

            pushMatrix(MatrixTypes::MODEL);


            this->transform();


            // send matrices to OGL

            glUniformMatrix4fv(view_uniformId, 1, GL_FALSE, mMatrix[MatrixTypes::VIEW]);

            if (entity == EntityType::SNOW){
                computeDerivedMatrix(VIEW_MODEL);
                BillboardCheatSphericalBegin(); 
                computeDerivedMatrix_PVM();
                
            }else computeDerivedMatrix(ComputedMatrixTypes::PROJ_VIEW_MODEL);

            glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::VIEW_MODEL]);
            glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::PROJ_VIEW_MODEL]);
            computeNormalMatrix3x3();
            glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

            // Render mesh
            glBindVertexArray(_mesh.vao);

            if (!shader->isProgramValid()) {
                printf("Program Not Valid!\n");
                exit(1);
            }


            if (entity == EntityType::ROAD) { // If it's a Cheerio...
                loc = glGetUniformLocation(shader->getProgramIndex(), "flag_instance");
                glUniform1i(loc, true);

                for (int i = 0; i < amount; i++) {

                    loc = glGetUniformLocation(shader->getProgramIndex(),
                        std::string("offsets[" + std::to_string(i) + "]").c_str());

                    glUniform3fv(loc, amount, cheercoor.at(i).data());
                }
                

                glDrawElementsInstanced(_mesh.type, _mesh.numIndexes, GL_UNSIGNED_INT, 0, amount);

            } else if (entity == EntityType::SNOW && snow) {
                updateParticles();

                loc = glGetUniformLocation(shader->getProgramIndex(), "flag_instance");
                glUniform1i(loc, true); 

                for (int i = 0; i < max_particles; i++) {
                    /*_mesh.mat.diffuse[0] = particles.at(i).r;
                    _mesh.mat.diffuse[1] = particles.at(i).g;
                    _mesh.mat.diffuse[2] = particles.at(i).b;*/

                    if (0.5f < _time_span) {
                        _mesh.mat.diffuse[3] -= 0.0004 * t_delta;
                    }
                    if (_mesh.mat.diffuse[3] < 0) {
                        _mesh.mat.diffuse[3] = 0;
                    }
                    loc = glGetUniformLocation(shader->getProgramIndex(), "mat.diffuse");
                    glUniform4fv(loc, 1, _mesh.mat.diffuse);

                    if (_time_span < 3.0f) {
                        loc = glGetUniformLocation(shader->getProgramIndex(),
                            std::string("offsets[" + std::to_string(i) + "]").c_str());

                        //const GLfloat* value = particles.at(i).position.data();

                        glUniform3fv(loc, max_particles, particles.at(i).position.data());
                    } else dead_num_particles = max_particles;
                }

                glDrawElementsInstanced(_mesh.type, _mesh.numIndexes, GL_UNSIGNED_INT, 0, max_particles);

            } else { // If it's not snow
                loc = glGetUniformLocation(shader->getProgramIndex(), "flag_instance");
                glUniform1i(loc, false);

                glDrawElements(_mesh.type, _mesh.numIndexes, GL_UNSIGNED_INT, 0);
            }


            if (dead_num_particles == max_particles) {
                snow = false;
                //glDepthMask(GL_TRUE);
                dead_num_particles = 0;
                printf("All particles dead\n");
            }

        } else if (entity == EntityType::TREES){
            loc = glGetUniformLocation(shader->getProgramIndex(), "text_flag");
            glUniform1i(loc, true);
            //loc = glGetUniformLocation(shader->getProgramIndex(), "particle_flag");
            //glUniform1i(loc, false);
            loc = glGetUniformLocation(shader->getProgramIndex(), "type_text_flag");
            glUniform1i(loc, 3);
            
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, this->_textureArray[0]);

            loc = glGetUniformLocation(shader->getProgramIndex(), "texMode");
            glUniform1i(loc, 1);
            loc = glGetUniformLocation(shader->getProgramIndex(), "texmap0");
            glUniform1i(loc, 0);

            //translate(MatrixTypes::MODEL, 15.0f, 0.0f, 15.0f);

            //float cam[2] = {camX, camZ};
            //float pos[2] = {15.0f, 15.0f};

            //if (type == 2)
			//	l3dBillboardSphericalBegin(cam, pos);
			//else if (type == 3)
			//	l3dBillboardCylindricalBegin(cam, pos);

            loc = glGetUniformLocation(shader->getProgramIndex(), "mat.specular");
			glUniform4fv(loc, 1, _mesh.mat.specular);
			loc = glGetUniformLocation(shader->getProgramIndex(), "mat.shininess");
			glUniform1f(loc, _mesh.mat.shininess);

            loc = glGetUniformLocation(shader->getProgramIndex(), "fogColor");
            vec_f fog_color({ 0.33f, 0.33f, 0.33f, 0.5f });
            glUniform4fv(loc, 1, fog_color.data());

            pushMatrix(MatrixTypes::MODEL);


            this->transform();

            //if (type == 0 || type == 1)	{     //Cheating matrix reset billboard techniques
			computeDerivedMatrix(VIEW_MODEL);
			
			//reset VIEW_MODEL
			//	if(type==0) BillboardCheatSphericalBegin();   
		    BillboardCheatCylindricalBegin();

			computeDerivedMatrix_PVM(); // calculate PROJ_VIEW_MODEL
			//}
			//else computeDerivedMatrix(ComputedMatrixTypes::PROJ_VIEW_MODEL);


            
            // send matrices to OGL
            //computeDerivedMatrix(ComputedMatrixTypes::PROJ_VIEW_MODEL);

            glUniformMatrix4fv(vm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::VIEW_MODEL]);
            glUniformMatrix4fv(pvm_uniformId, 1, GL_FALSE, mCompMatrix[ComputedMatrixTypes::PROJ_VIEW_MODEL]);
            computeNormalMatrix3x3();
            glUniformMatrix3fv(normal_uniformId, 1, GL_FALSE, mNormal3x3);

            // Render mesh
            glBindVertexArray(_mesh.vao);

            
            if (!shader->isProgramValid()) {
                printf("Program Not Valid!\n");
                exit(1);
            }

            glDrawElements(_mesh.type, _mesh.numIndexes, GL_UNSIGNED_INT, 0);

        }

        glBindVertexArray(0);
        popMatrix(MatrixTypes::MODEL);


        glBindTexture(GL_TEXTURE_2D, 0);
    }

    void transform() {
        for (int ix = _transfs.size() - 1; ix >= 0; ix--) {
            _transfs.at(ix).apply();
        }
    }

    void update() {
        for (auto&& transf : _transfs) {
            if (transf.is_modifiable()) {
                switch (transf.get_type()) {
                case TransType::Translate:
                    transf.get_comps()->at(0) = _offset.at(0);
                    transf.get_comps()->at(1) = _offset.at(1);
                    transf.get_comps()->at(2) = _offset.at(2);
                    break;
                }
            }
        }
    }

    void apply_speed(float speed, float y_angle) {
        size_t comp = 0;
        vec_f vel({ static_cast<float>(speed * t_delta * cos(y_angle * PI / 180.0)), 0.0f,
                   static_cast<float>(speed * t_delta * sin(y_angle * PI / 180.0)) });

        for (auto vel_comp : vel) {
            _offset.at(comp) += vel_comp;
            comp += 1;
        }

        for (auto&& transf : _transfs) {
            if (transf.get_type() == TransType::RotateY && transf.is_modifiable()) {
                transf.get_comps()->at(0) = -y_angle;
            }
        }

        _y_angle = y_angle;
    }

    void apply_rotation(float rot_amount) {
        for (auto&& transf : _transfs) {
            if (transf.is_modifiable() && transf.get_type() == TransType::RotateZ) {
                transf.get_comps()->at(0) -= rot_amount;
            }
        }
    }
};

class Entity {
protected:
    std::vector<Object> _objects;
public:
    Entity() {}

    Entity(std::vector<Object> objects) : _objects(objects) {}

    std::vector<Object>* get_objects() {
        return &_objects;
    }

    virtual void input() {}

    virtual void reset(bool hard_reset) {}

    virtual void update() {}
};

class Car : public Entity {
    float _speed = 0;
    float _angle = 0;
    float _accel = 0;

    bool _rolling;
    Direction _direction;

    float _max_accel = 1.5f;
    float _friction = -2.8f;
    float _max_speed = 17.0f;
public:
    Car(std::vector<Object> objects) {
        _objects = objects;
    }

    void reset(bool hard_reset) {
        for (auto&& obj : this->_objects) {
            obj.set_offset(vec_f({ 0.0f, 0.0f, 0.0f }));
            obj.set_y_angle(0.0f);
            obj.set_transf_at(obj.get_transfs()->size() - 1, vec_f({ 0.0f, 0.0f, 0.0f }));
            _speed = 0.0f;
        }
        this->set_y_angle(0.0f);
    }

    std::vector<float> get_coords() {
        return _objects.at(0).get_coords();
    }

    float get_speed() {
        return _speed;
    }

    void set_speed(float speed) {
        _speed = speed;
    }

    void set_y_angle(float y_angle) {
        _angle = y_angle;
    }

    float get_angle() {
        return _angle;
    }

    void check_off_table() {
        auto coords = this->get_objects()->at(0).get_coords();
        if (coords.at(0) < 75.0f && coords.at(2) > -75.0f &&
            coords.at(0) > -75.0f && coords.at(2) < 75.0f) {

            // Car is on the table, all good
        }
        else {
            points_n = points_n - 100;

            _to_reset = true;;
        }
    }

    void input() override {
        int state;
        _rolling = true;


        if (_speed == 0) { // Facing forward if stopped
            _direction = Direction::FORWARD;
        }
        if (_direction == Direction::FORWARD) { // Friction affects in the opposite direction of movement
            _accel = _friction;
        }
        else {
            _accel = -_friction;
        }
        if (car_forwards) { //car goes forward
            reset_cam();
            _accel = _max_accel * 4;
            _rolling = false;

            if (_speed > 0.5f) {
                _direction = Direction::FORWARD;
            }

            if (_direction == Direction::BACKWARDS) {
                _accel = _max_accel * 10;
            }
        }
        if (car_backwards) { // car goes backwards
            reset_cam();
            _accel = -_max_accel * 1.5;
            _rolling = false;

            if (_speed < 0.5f && _speed > -0.5f) {
                _direction = Direction::BACKWARDS;
            }

            if (_direction == Direction::FORWARD) {
                _accel = -_max_accel * 10;
            }
        }
        if (car_left) { // car goes left
            if (_direction == Direction::FORWARD) {
                _angle -= 60 * t_delta;
            }
            else {
                _angle += 60 * t_delta;
            }
            reset_cam();
        }
        if (car_right) { //car goes right
            if (_direction == Direction::FORWARD) {
                _angle += 60 * t_delta;
            }
            else {
                _angle -= 60 * t_delta;
            }
            reset_cam();
        }

        _speed = _speed + _accel * t_delta;
        if (_speed < 0.5f && _speed > -0.5f && _rolling) {
            _speed = 0;
        }

        if (_speed > _max_speed) { // Limit max speed (be it forwards or backwards)
            _speed = _max_speed;
        }
        else if (_speed < (-_max_speed / 3)) {
            _speed = -_max_speed / 3;
        }


        if (_speed < 0.5f && _speed > -0.5f && _rolling) {
            points_n = points_n - 0.5f * t_delta;
        } else if (_direction == Direction::FORWARD && _rolling) {
            points_n = points_n + 0.5f * t_delta;
        } else if (_direction == Direction::FORWARD && !_rolling) {
            points_n = points_n + 1.0f * t_delta;
        } else if (_direction == Direction::BACKWARDS && _rolling) {
            points_n = points_n - 1.0f * t_delta;
        } else if (_direction == Direction::BACKWARDS && !_rolling) {
            points_n = points_n - 0.5f * t_delta;
        }
    }

    void update() override {
        this->check_off_table();
        for (auto&& object : *this->get_objects()) {
            object.apply_speed(_speed, _angle);
            object.update();
        }
    }
};

class Road : public Entity {
public:
    Road(std::vector<Object> objects) {
        _objects = objects;
    }
};


class Butter : public Entity {
    std::vector<float> _speeds;

    float _friction = -2.8f;

public:
    Butter(std::vector<Object> objects) {
        _objects = objects;

        _speeds.resize(objects.size());
        for (auto&& speed : _speeds) {
            speed = 0.0f;
        }
    }

    void reset(bool hard_reset) {
        if (hard_reset) {
            for (auto&& object : this->_objects) {
                object.get_transfs()->at(object.get_transfs()->size() - 1).get_comps()->at(0) = 0.0f;
                object.get_transfs()->at(object.get_transfs()->size() - 1).get_comps()->at(1) = 0.0f;
                object.get_transfs()->at(object.get_transfs()->size() - 1).get_comps()->at(2) = 0.0f;
                object.set_offset(vec_f({ 0.0f, 0.0f, 0.0f }));
            }
            for (auto&& speed : _speeds) {
                speed = 0.0f;
            }
            this->update();
        }
    }

    float get_speed(size_t ix) {
        return _speeds.at(ix);
    }

    void set_speed(size_t ix, float speed) {
        _speeds.at(ix) = speed;
    }

    void update() override {
        size_t butter_ix = 0;
        for (auto&& object : *this->get_objects()) {
            _speeds.at(butter_ix) = _speeds.at(butter_ix) + _friction * t_delta;
            float speed = _speeds.at(butter_ix);
            if (speed < 0.5f) {
                _speeds.at(butter_ix) = speed = 0;
            }
            object.apply_speed(speed, object.get_y_angle());
            object.update();

            butter_ix++;
        }
    }

    void update_at(size_t butter_ix, float y_angle) {
        this->get_objects()->at(butter_ix).apply_speed(_speeds.at(butter_ix), y_angle);
    }
};

class Orange : public Entity {
    std::vector<float> _speeds;
    std::vector<float> _y_angles;
    std::vector<float> _time;

    double _rel_game_time = 0;


public:
    Orange(std::vector<Object> objects) {
        _objects = objects;
        _speeds = vec_f(objects.size(), 0.0f);
        _y_angles = vec_f(objects.size(), 0.0f);
        _time = vec_f(objects.size(), 0.0f);

        for (auto&& speed : _speeds) {
            speed = rand_f(15.0);
        }
        for (auto&& y_angle : _y_angles) {
            y_angle = rand_f(360);
        }
    }

    void reset(bool hard_reset) {
        if (hard_reset) {
            for (auto&& speed : _speeds) {
                speed = rand_f(15.0);
            }
            for (auto&& y_angle : _y_angles) {
                y_angle = rand_f(360);
            }
        }
    }

    void check_off_table() {
        size_t orange_ix = 0;
        for (auto&& orange : *this->get_objects()) {
            vec_f coords = orange.get_coords();
            if (coords.at(0) < 75.0f && coords.at(2) > -75.0f &&
                coords.at(0) > -75.0f && coords.at(2) < 75.0f) {

                // Orange is on the table, all good
            }
            else {
                orange.set_offset(vec_f({ rand_f(70), 0.0f, rand_f(70) }));
                //                _speeds.at(orange_ix) = rand_f(3.0);
                _y_angles.at(orange_ix) = rand_f(360);

                orange.set_drawability(false);
                _time.at(orange_ix) = rand_uf(2.0f);
            }
            orange_ix++;
        }
    }

    void update() override {
        _rel_game_time += t_delta;

        this->check_off_table();
        size_t orange_ix = 0;
        for (auto&& orange : *this->get_objects()) {
            if (orange.is_drawable()) {
                if (!((int)round(_rel_game_time) % 20)) {
                    _speeds.at(orange_ix) *= 1.5f;
                }
                float speed = _speeds.at(orange_ix);
                float y_angle = _y_angles.at(orange_ix);
                orange.apply_speed(speed, y_angle);

                float radius = 1.5f;
                float rot_amount = speed / (radius * PI);
                orange.apply_rotation(rot_amount);

                orange.update();
            }
            orange_ix++;
        }
        if (!((int)round(_rel_game_time) % 20)) {
            _rel_game_time += 1;
        }

        orange_ix = 0;
        for (auto&& time : _time) {
            time -= t_delta;
            if (time <= 0.0f) {
                _objects.at(orange_ix).set_drawability(true);
                time = 0.0f;
            }
            orange_ix++;
        }
    };
};


class Scene {
    std::vector<std::unique_ptr<Entity>> entities;

    Camera _cam = Camera::FOLLOW;

public:
    Scene() {

    }

    std::vector<std::unique_ptr<Entity>>* get_entities() {
        return &entities;
    }

    void reset() {
        bool hard_reset = false;
        lifes_left = lifes_left - 1;

        if (lifes_left <= 0) {
            _game_over = true;
        }

        if (_reset) {
            hard_reset = true;
            lifes_left = 5;
            points_n = 0;
            _game_over = false;
            _reset = false;
        }

        for (auto&& entity : this->entities) {
            entity->reset(hard_reset);
        }
    }

    void cam_update() {
        std::vector<float> car_coords = entities.at(EntityType::CAR)->get_objects()->at(0).get_coords();
        float ratio;
        int h = HEIGHT;
        int w = WIDTH;
        if (h == 0)
            h = 1;
        ratio = (1.0f * w) / h;

        switch (_cam) {
        case Camera::FOLLOW:
            changeSize(WIDTH, HEIGHT);
            // set the camera position based on its spherical coordinates
            camX = rAux * sin(alphaAux * PI / 180.0f) * cos(betaAux * PI / 180.0f) + car_coords.at(0);
            camZ = rAux * cos(alphaAux * PI / 180.0f) * cos(betaAux * PI / 180.0f) + car_coords.at(2);
            camY = rAux * sin(betaAux * PI / 180.0f);
            lookAt(camX, camY, camZ, car_coords.at(0), car_coords.at(1), car_coords.at(2), 0, 1, 0);
            changeSize(WIDTH, HEIGHT);
            break;
        case Camera::TOP_ORTHO:
            camX = rStatic * sin(PI / 180.0f);
            camZ = rStatic * cos(PI / 180.0f);
            camY = rStatic * sin(PI / 180.0f) + 90;
            loadIdentity(MatrixTypes::PROJECTION);
            ortho(-(160 * ratio) / 2, (160 * ratio) / 2,
                -160 / 2, 160 / 2,
                1, 1000);
            lookAt(camX, camY, camZ, 0, 0, 0, 0, 0, -1);
            break;
        case Camera::TOP_PERS:
            changeSize(WIDTH, HEIGHT);

            camX = rStatic * sin(PI / 180.0f);
            camZ = rStatic * cos(PI / 180.0f);
            camY = rStatic * sin(PI / 180.0f) + 160;
            lookAt(camX, camY, camZ, 0, 0, 0, 0, 0, -1);
            changeSize(WIDTH, HEIGHT);
            break;
        }
    }

    void key_input(int key, int xx, int yy) {
        if (key == '1') {
            _cam = Camera::TOP_ORTHO;
        }
        else if (key == '2') {
            _cam = Camera::TOP_PERS;
        }
        else if (key == '3') {
            _cam = Camera::FOLLOW;
        }
        else if (key == 'h') { // TODO Car headlights
            puts("H");
        }
        else if (key == 'e') {
            puts("E");
		    iniParticles();
            part_timer = std::chrono::high_resolution_clock::now();
		    snow = true;
        }
        else if (key == 'r') {
            puts("R");
            _reset = true;
            this->reset();
        }
        else if (key == 's') {
            puts("S");
            _pause = !_pause;
        }
        else if (key == 'b') {
            puts("B");
            bumpmap = !bumpmap;
        }
        else if (key == 'f') {
            puts("F - Flare");
            if (spotlight_mode) flareEffect = false;
			else
				if (flareEffect) flareEffect = false;
				else flareEffect = true;
        }
        else if (key == 'n') {
            puts("N");
            if (!directional) {
                directional = true;

                printf("Point light disabled. Spot light enabled");
            }
            else {
                directional = false;
                printf("Spot light disabled. Point light enabled");
            }

        } else if (key == 'c') { // TODO 6 candle lights
            puts("C");
            if (!spotlight_mode) {
                spotlight_mode = true;

                printf("Point light disabled. Spot light enabled");
            }
            else {
                spotlight_mode = false;
                printf("Spot light disabled. Point light enabled");
            }
        }
    }

    void push_all_transform(Transformation transf) {
        for (auto&& entity : entities) {
            for (auto&& object : *entity->get_objects()) {
                object.get_transfs()->push_back(transf);
            }
        }
    }
    
    void pop_all_transform() {
        for (auto&& entity : entities) {
            for (auto&& object : *entity->get_objects()) {
                object.get_transfs()->pop_back();
            }
        }
    }


    void input() {
        for (auto&& entity : entities) {
            entity->input();
        }
    }

    void do_collision_butter(Car* car, Butter* butters) {
        Object body = car->get_objects()->at(0);
        size_t butter_ix = 0;
        for (auto&& butter : *butters->get_objects()) {
            if (check_collision_AA(body, butter)) {
                butters->set_speed(butter_ix, std::abs(car->get_speed()) < 2.0f ? 5.0f : std::abs(car->get_speed()));
                butters->update_at(butter_ix, car->get_speed() < 0.0f ? car->get_angle() + 180 : car->get_angle());
                car->set_speed(0.0f);
                car->update();

                points_n = points_n - 9.0f;
            }
            butter_ix++;
        }
    }

    void do_collision_cheer(Car* car, Entity* entity) {
        Object body = car->get_objects()->at(0);

        for (auto&& cheerio : cheercoor) {
            if (check_collision_CH(body, cheerio, 0.1f)) {
                vec_f vel{ static_cast<float>(600.0f * t_delta * cos(car->get_angle() * PI / 180.0)), 0.0f,
                          static_cast<float>(600.0f * t_delta * sin(car->get_angle() * PI / 180.0)) };

                size_t ix = 0;
                for (auto&& comp : vel) {
                    cheerio.at(ix) += comp;
                    ix += 1;
                }
                car->set_speed(0.0f);
                car->update();

                points_n = points_n - 5.0f;
            }
        }
    }

    void do_collision_orange(Car* car, Orange* oranges) {
        Object body = car->get_objects()->at(0);
        size_t orange_ix = 0;
        for (auto&& orange : *oranges->get_objects()) {
            if (oranges_kill && orange.is_drawable() && check_collision_SA(body, orange, 1.5f)) {
                car->set_speed(0.0f);
                this->reset();

                points_n = points_n/4.0f;
            }
            orange_ix++;
        }
    }

    void draw(VSShaderLib* shader) {
        size_t entity_ix = 0;
        Car* car = dynamic_cast<Car*>(entities.at(EntityType::CAR).get());
        Butter* butter = dynamic_cast<Butter*>(entities.at(EntityType::BUTTER).get());
        Entity* road = entities.at(EntityType::ROAD).get();
        Orange* orange = dynamic_cast<Orange*>(entities.at(EntityType::ORANGES).get());


        part_timer_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(part_timer_end - part_timer);
        auto _time_span = time_span.count();

        // Draw objects
        for (auto&& entity : entities) {
            if (entity_ix != EntityType::TABLE) {
                for (auto&& object : *entity->get_objects()) {
                    if (object.is_drawable()) {
                        glEnable(GL_BLEND);
                        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        object.draw(shader, entity_ix, _time_span);
                        glDisable(GL_BLEND);
                    }
                }
            }
            entity_ix++;
        }

        //render_skybox(shader, &entities.at(EntityType::SKYBOX).get()->get_objects()->at(0));

        // TODO cheerios


        // TODO 


        if (flareEffect && !spotlight_mode) {

		    int flarePos[2];
		    int m_viewport[4];
		    glGetIntegerv(GL_VIEWPORT, m_viewport);

		    pushMatrix(MODEL);
		    loadIdentity(MODEL);
		    computeDerivedMatrix(PROJ_VIEW_MODEL);  //pvm to be applied to lightPost. pvm is used in project function
		
		    if (!project(lightPos, lightScreenPos, m_viewport))
			    printf("Error in getting projected light in screen\n");  //Calculate the window Coordinates of the light position: the projected position of light on viewport
		    flarePos[0] = clampi((int)lightScreenPos[0], m_viewport[0], m_viewport[0] + m_viewport[2] - 1);
		    flarePos[1] = clampi((int)lightScreenPos[1], m_viewport[1], m_viewport[1] + m_viewport[3] - 1);
		    popMatrix(MODEL);

		       //viewer looking down at  negative z direction
		    pushMatrix(PROJECTION);
		    loadIdentity(PROJECTION);
		    pushMatrix(VIEW);
		    loadIdentity(VIEW);
		    ortho(m_viewport[0], m_viewport[0] + m_viewport[2] - 1, m_viewport[1], m_viewport[1] + m_viewport[3] - 1, -1, 1);
		    render_flare(&AVTflare, flarePos[0], flarePos[1], m_viewport); 
            popMatrix(PROJECTION);
		    popMatrix(VIEW);
	    }

        render_text();

        do_collision_butter(car, butter);
        do_collision_cheer(car, road);
        do_collision_orange(car, orange);

    }

    void update() {
        for (auto&& entity : entities) {
            entity->update();
        }
    }

    void forge_mesh(MyMesh* mesh, std::vector<float> _amb, std::vector<float> _diff, std::vector<float> _spec,
        std::vector<float> _emissive, float _shininess, int _texcount) {
        memcpy(mesh->mat.ambient, _amb.data(), 4 * sizeof(float));
        memcpy(mesh->mat.diffuse, _diff.data(), 4 * sizeof(float));
        memcpy(mesh->mat.specular, _spec.data(), 4 * sizeof(float));
        memcpy(mesh->mat.emissive, _emissive.data(), 4 * sizeof(float));
        mesh->mat.shininess = _shininess;
        mesh->mat.texCount = _texcount;
    }

    std::vector<Object> create_table() {
        std::vector<Object> table;
        std::vector<GLuint> TextureArray(2);


        glGenTextures(2, TextureArray.data());
        //Texture2D_Loader(TextureArray.data(), "./lightwood.tga", 0);
        //Texture2D_Loader(TextureArray.data(), "./spiderwebs.jpeg", 1);

        vec_f _amb({ 0.1f, 0.1f, 0.1f, 1.0f });
        //vec_f _diff({ 0.99f, 0.99f, 0.99f, 0.10f });
        vec_f _spec({ 0.9f, 0.9f, 0.9f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f});
        float _shininess = 100.0f;
        int _texcount = 0;

        vec_f color({ 1.00f, 0.00f, 0.00f, 0.1f });

        MyMesh mesh;
        mesh = createCube(); //table top
        this->forge_mesh(&mesh, _amb, color, _spec,
            _emissive, _shininess, _texcount);

        vec_t transfs;
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -0.5f, -1.0f, -0.5f }), false));
        transfs.push_back(Transformation(TransType::Scale, vec_f({ 150.0f, 1.5f, 150.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));

        table.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));

        transfs.clear();
        // legs
        // 1
       /*  mesh = createCylinder(15.0f, 1.0f, 20);
        this->forge_mesh(&mesh, _amb, color, _spec,
            _emissive, _shininess, _texcount);

        transfs.push_back(Transfo rmation(TransType::Scale, vec_f({ 1.0f, 5.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 74.0f, -39.0f, 74.0f }), false));
        table.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.0f, 0.0f, 0.0f }), TextureArray));

        transfs.clear();
        // 2
        mesh = createCylinder(15.0f, 1.0f, 20);
        this->forge_mesh(&mesh, _amb, color, _spec,
            _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 1.0f, 5.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 74.0f, -39.0f, -74.0f }), false));
        table.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.0f, 0.0f, 0.0f }), TextureArray));

        transfs.clear();
        // 3
        mesh = createCylinder(15.0f, 1.0f, 20);
        this->forge_mesh(&mesh, _amb, color, _spec,

            _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 1.0f, 5.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -74.0f, -39.0f, -74.0f }), false));
        table.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.0f, 0.0f, 0.0f }), TextureArray));

        transfs.clear();

        // 4
        mesh = createCylinder(15.0f, 1.0f, 20);
        this->forge_mesh(&mesh, _amb, color, _spec,
            _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 1.0f, 5.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -74.0f, -39.0f, 74.0f }), false));
        table.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.0f, 0.0f, 0.0f }), TextureArray));
 */
        //transfs.clear();
        return table;
    }

    std::vector<Object> create_car() {
        std::vector<Object> car;

        vec_f _amb({ 0.2f, 0.15f, 0.1f, 1.0f });
        vec_f _diff({ 0.8f, 0.6f, 0.4f, 1.0f });
        vec_f _spec({ 0.8f, 0.8f, 0.8f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });
        float _shininess = 100.0f;
        int _texcount = 0;

        vec_f color1({ 0.8f, 1.0f, 1.0f, 0.2f });
        MyMesh mesh;

        vec_t transfs;

        mesh = createCube(); //body
        this->forge_mesh(&mesh, _amb, color1, _spec,
            _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 4.0f, 0.8f, 2.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -2.0f, 0.2f, -1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));
        car.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.2f, 62.0f }),
            vec_f({ 4.0f, 0.8f, 2.0f })));

        transfs.clear();

        mesh = createCube(); //windows
        this->forge_mesh(&mesh, _amb, color1, _spec,
            _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 1.0f, 0.5f, 2.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -1.0f, 1.0f, -1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));
        car.push_back(Object(mesh, transfs, vec_f({ -0.5f, 1.0f, 61.0f }), vec_f({ 1.0f, 0.5f, 2.0f }), {}));

        transfs.clear();


        mesh = createCylinder(2.0f, 0.5f, 3); //windowshield
        this->forge_mesh(&mesh, _amb, color1, _spec,
            _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 1.0f, 1.0f, 3.0f }), false));
        transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 90, 0.0f, 0.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 90, 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f })));

        car.push_back(
            Object(mesh, transfs, vec_f({ 0.0f, 1.0f,62.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));


        transfs.clear();

        mesh = createCylinder(2.0f, 0.5f, 3); //windowshield
        this->forge_mesh(&mesh, _amb, color1, _spec,
            _emissive, _shininess, _texcount);

        //transfs.push_back(Transformation(TransType::Scale, vec_f({{1.0f, 1.0f, 3.0f}})));
        transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 90, 0.0f, 0.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 90, 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -1.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f })));

        car.push_back(Object(mesh, transfs, vec_f({ -1.0f, 1.0f, 62.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));

        transfs.clear();

        vec_f amb2({ 0.2f, 0.2f, 0.2f, 1.0f });
        vec_f color2({ 1.0f, 1.0f, 1.0f, 1.0f });
        vec_f spec2({ 0.9f, 0.9f, 0.9f, 1.0f });
        vec_f emissive2({ 0.3f, 0.3f, 0.3f, 1.0f });

        //wheel1
        mesh = createTorus(0.10f, 0.4f, 20, 20);
        this->forge_mesh(&mesh, amb2, color2, spec2,
            emissive2, _shininess, _texcount);
        transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 90.0f, 0.0f, 0.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 90.0f, 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -1.5f, 0.43f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));
        car.push_back(Object(mesh, transfs, vec_f({ -1.5f, 0.43f, 63.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));

        transfs.clear();

        //wheel2
        mesh = createTorus(0.10f, 0.4f, 20, 20);
        this->forge_mesh(&mesh, amb2, color2, spec2,
            emissive2, _shininess, _texcount);
        transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 90, 0.0f, 0.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 90, 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -1.5f, 0.43f, -1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));

        car.push_back(Object(mesh, transfs, vec_f({ -1.5f, 0.43f, 61.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));

        transfs.clear();

        //wheel3
        mesh = createTorus(0.10f, 0.4f, 20, 20);
        this->forge_mesh(&mesh, amb2, color2, spec2,
            emissive2, _shininess, _texcount);
        transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 90, 0.0f, 0.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 90, 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 1.5f, 0.43f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));

        car.push_back(Object(mesh, transfs, vec_f({ 1.5f, 0.43f, 63.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));

        transfs.clear();

        //wheel4
        mesh = createTorus(0.10f, 0.4f, 20, 20);
        this->forge_mesh(&mesh, amb2, color2, spec2,
            emissive2, _shininess, _texcount);
        transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 90, 0.0f, 0.0f, 1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 90, 0.0f, 1.0f, 0.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 1.5f, 0.43f, -1.0f }), false));
        transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 62.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));


        car.push_back(Object(mesh, transfs, vec_f({ 1.5f, 0.43f, 61.0f }), vec_f({ 0.0f, 0.0f, 0.0f })));

        transfs.clear();

        return car;
    }


    std::vector<Object> create_road() {
        std::vector<Object> road;
        size_t amount = 3;
        MyMesh mesh;
        vec_t transfs;

        vec_f _amb({ 0.2f, 0.15f, 0.1f, 1.0f });
        vec_f _diff({ 0.8f, 0.6f, 0.4f, 1.0f });
        vec_f _spec({ 0.8f, 0.8f, 0.8f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });
        float _shininess = 100.0f;
        int _texcount = 0;

        vec_f color({ 0.0f, 1.0f, 0.00f, 1.0f });

        mesh = createTorus(0.7f, 0.4f, 20, 20);
        this->forge_mesh(&mesh, _amb, color, _spec,
            _emissive, _shininess, _texcount);


        road.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.8f, 0.8f, 0.8f })));
        transfs.clear();


        return road;
    }

    std::vector<Object> create_oranges() {
        std::vector<Object> oranges;
        std::vector<Transformation> transfs;
        MyMesh mesh;
        std::vector<GLuint> TextureArray(2);


        glGenTextures(2, TextureArray.data());
        Texture2D_Loader(TextureArray.data(), "./jack_o_lantern.png", 0);
        Texture2D_Loader(TextureArray.data(), "./orange.jpeg", 1);

        vec_f _amb({ 0.2f, 0.15f, 0.1f, 1.0f });
        vec_f _diff({ 0.8f, 0.6f, 0.4f, 1.0f });
        vec_f _spec({ 0.8f, 0.8f, 0.8f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });
        float _shininess = 100.0f;
        int _texcount = 2;

        vec_f _color({ 1.0f, 0.5f, 0.0f, 1.0f });

        for (size_t n = 1; n <= 15; n++) {
            float x = rand_f(70), z = rand_f(70);

            mesh = createSphere(1.5f, 40);
            this->forge_mesh(&mesh, _amb, _color, _spec,
                _emissive, _shininess, _texcount);

            transfs.push_back(Transformation(TransType::RotateZ, vec_f({ 0.0f, 0.0f, 0.0f, 1.0f })));
            transfs.push_back(Transformation(TransType::RotateY, vec_f({ 0.0f, 0.0f, 1.0f, 0.0f })));
            transfs.push_back(Transformation(TransType::Translate, vec_f({ x, 1.5f, z }), false));
            transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));

            oranges.push_back(
                Object(mesh, transfs, vec_f({ x, 0.0f, z }), vec_f({ 0.0f, 0.0f, 0.0f }), TextureArray));
            transfs.clear();
        }

        return oranges;
    }

    std::vector<Object> create_butter() {
        std::vector<Object> butter;
        std::vector<Transformation> transfs;
        MyMesh mesh;

        vec_f _amb({ 0.2f, 0.15f, 0.1f, 1.0f });
        vec_f _diff({ 0.8f, 0.6f, 0.4f, 1.0f });
        vec_f _spec({ 0.8f, 0.8f, 0.8f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });
        float _shininess = 100.0f;
        int _texcount = 2;

        vec_f color({ 1.0f, 1.0f, 0.8f, 1.0f });
        std::vector<GLuint> TextureArray(2);

        glGenTextures(2, TextureArray.data());
        Texture2D_Loader(TextureArray.data(), "./stone.tga", 0);
        Texture2D_Loader(TextureArray.data(), "./normal.tga", 1);

        for (size_t n = 1; n <= 30; ++n) {
            float x = rand_f(70), z = rand_f(70), y_angle = rand_f(90);

            mesh = createCube();
            this->forge_mesh(&mesh, _amb, color, _spec,
                _emissive, _shininess, _texcount);
            transfs.push_back(Transformation(TransType::Translate, vec_f({ -0.5f, 0.0f, -0.5f }), false));
            transfs.push_back(Transformation(TransType::Scale, vec_f({ 4.0f, 1.3f, 2.0f }), false));
            transfs.push_back(Transformation(TransType::RotateY, vec_f({ y_angle, 0.0f, 1.0f, 0.0f }), false));
            transfs.push_back(Transformation(TransType::Translate, vec_f({ x, 0.0f, z }), false));
            transfs.push_back(Transformation(TransType::Translate, vec_f({ 0.0f, 0.0f, 0.0f })));

            //butter.push_back(Object(mesh, transfs, vec_f({ x, 0.0f, z }), vec_f({ 4.0f, 1.3f, 2.0f }), TextureArray));
            transfs.clear();
        }


        return butter;
    }

    std::vector<Object> create_snow() {
        std::vector<Object> snow_objs;
        MyMesh mesh;
        std::vector<Transformation> transfs;

        vec_f _amb({ 0.2f, 0.15f, 0.1f, 1.0f });
        vec_f _diff({ 0.8f, 0.6f, 0.4f, 1.0f });
        vec_f _spec({ 0.8f, 0.8f, 0.8f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });
        float _shininess = 100.0f;
        int _texcount = 1;

       // vec_f color({ 0.882f, 0.552f, 0.211f, 1.0f });

        std::vector<GLuint> TextureArray(1);
        
        glGenTextures(1, TextureArray.data());
        Texture2D_Loader(TextureArray.data(), "./snowflake.jpeg", 0);


        mesh = createQuad(2, 2);
        this->forge_mesh(&mesh, _amb, _diff, _spec,
           _emissive, _shininess, _texcount);

        

        snow_objs.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.8f, 0.8f, 0.8f }), TextureArray));
        transfs.clear();

        return snow_objs;
    }

    std::vector<Object> create_trees() {
        std::vector<Object> obj;
        MyMesh mesh;
        std::vector<Transformation> transfs;

        vec_f _amb({ 0.2f, 0.15f, 0.1f, 1.0f });
        vec_f _diff({ 0.8f, 0.6f, 0.4f, 1.0f });
        vec_f _spec({ 0.2f, 0.2f, 0.2f, 1.0f });
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });
        float _shininess = 10.0f;
        int _texcount = 1;

        std::vector<GLuint> TextureArray(1);
        
        glGenTextures(1, TextureArray.data());
        Texture2D_Loader(TextureArray.data(), "./tree.tga", 0);


        mesh = createQuad(6, 6);
        this->forge_mesh(&mesh, _amb, _diff, _spec,
           _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Translate, vec_f({ 7.0f, 3.0f, 7.0f }), false));
        obj.push_back(Object(mesh, transfs, vec_f({ 0.0f, 3.0f, 0.0f }), vec_f({ 0.8f, 0.8f, 0.8f }), TextureArray));
        transfs.clear();

        return obj;
    }

    std::vector<Object> create_skybox() {
        std::vector<Object> obj;
        MyMesh mesh;
        std::vector<Transformation> transfs;

        vec_f _amb({0.2f, 0.15f, 0.1f, 1.0f});
        vec_f _diff({0.8f, 0.6f, 0.4f, 1.0f});
        vec_f _spec({0.8f, 0.8f, 0.8f, 1.0f});
        vec_f _emissive({ 0.0f, 0.0f, 0.0f, 1.0f });

        float _shininess = 10.0f;
        int _texcount = 1;

        std::vector<GLuint> TextureArray(1);
        
        glGenTextures(1, TextureArray.data());

        const char *filenames[] = { "posx.jpg", "negx.jpg", "posy.jpg", "negy.jpg", "posz.jpg", "negz.jpg" };

	    TextureCubeMap_Loader(TextureArray.data(), filenames, 0);


        mesh = createCube();
        this->forge_mesh(&mesh, _amb, _diff, _spec,
           _emissive, _shininess, _texcount);

        transfs.push_back(Transformation(TransType::Scale, vec_f({ 200.0f, 200.0f, 200.0f }), false));
        transfs.push_back(Transformation(TransType::Translate, vec_f({ -0.5f, -0.5f, -0.5f }), false));
        //obj.push_back(Object(mesh, transfs, vec_f({ 0.0f, 0.0f, 0.0f }), vec_f({ 0.8f, 0.8f, 0.8f }), TextureArray));
        transfs.clear();

        return obj;
    }


    void create_objs() {
        entities.emplace_back(new Entity(create_table()));
        entities.emplace_back(new Orange(create_oranges()));
        entities.emplace_back(new Butter(create_butter()));
        //entities.emplace_back(new Entity(create_skybox()));
        entities.emplace_back(new Entity(create_trees()));
        entities.emplace_back(new Entity(create_snow()));
        entities.emplace_back(new Entity(create_road()));
        entities.emplace_back(new Car(create_car()));
    }
};
