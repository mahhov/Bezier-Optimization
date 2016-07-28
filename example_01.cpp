#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

using namespace std;

#define PI 3.14159265 // Should be used from mathlib
inline float sqr(float x) { return x*x; }

/*
WARNING: this is my first c++ program, so there are probably a lot of "bad" coding.
also, this program creates a lot of objects via the "new" keword, but deallocates very few of them
i.e, it has a very huge memory waste, but modern computers should reclaim the memory once the program is closed
*/

class Viewport;
class Patch;
class Camera;
class Point;
class Model;

void drawUniform();
void drawAdaptive();
void drawControl();
void initScene();
void myReshape(int, int);
void myDisplay();
void myFrameMove();
void drawControl();
void drawAdaptive();
void drawUniform();
void initCommands(int, char*[]);
void readInputFile();
void myInput(unsigned char, int, int);
void mySpecialInput(int, int, int);
void myMouseInput(int, int, int, int);
void myMotionInput(int, int);

class Viewport {
public:
	int w, h; // width and height
};

class Point {
public:
	float x, y, z;

	Point(float x, float y, float z, bool normalize = false) {
		if (normalize) {
			float mag = sqrt(x*x+y*y+z*z);
			this->x = x / mag;
			this->y = y / mag;
			this->z = z / mag;
		} else {
			this->x = x;
			this->y = y;
			this->z = z;
		}
	};

	Point* merge(Point* p, float u) {
		float v = 1 - u;
		float x = this->x * v + u * p->x;
		float y = this->y * v + u * p->y;
		float z = this->z * v + u * p->z;
		return new Point(x, y, z);
	}

	Point** merge(Point* p0, Point* p1, Point* p2, Point* p3, float u) {
		float v = 1 - u;
		Point* a = p0->merge(p1, u);
		Point* b = p1->merge(p2, u);
		Point* c = p2->merge(p3, u);
		a = a->merge(b, u);
		b = b->merge(c, u);

		Point* p = a->merge(b, u);
		b = new Point(3 * (b->x - a->x), 3 * (b->y - a->y), 3 * (b->z - a->z));

		delete a;
		delete c;

		Point** d = new Point*[2];
		d[0] = p;
		d[1] = b;
		return d;
	};

	Point* cross(Point* p, bool normal){
		return new Point(y*p->z - z*p->y, z*p->x - x*p->z, x*p->y - y*p->x, normal);	
	};

	Point* normal(Point* p1, Point* p2) {
		float dx1 = p1->x - x;
		float dx2 = p2->x - x;
		float dy1 = p1->y - y;
		float dy2 = p2->y - y;
		float dz1 = p1->z - z;
		float dz2 = p2->z - z;
		return new Point(dy1*dz2 - dz1*dy2, dz1*dx2 - dx1*dz2, dx1*dy2 - dy1*dx2);	
	};

	float distance2(Point* p) {
		float dx = p->x - x;
		float dy = p->y - y;
		float dz = p->z - z;
		return dx*dx + dy*dy + dz*dz;
	};
};

class Patch {
public:
	Point *p[4][4];
	Point ***c;
	Point ***n;
	vector<Point***> tri;

	// uniform
	Patch(Point *p[][4], int subdivision, float step) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				this->p[i][j] = p[i][j];
			}
		}
		uniformSubdivide(subdivision, step);
	};

	// adaptive
	Patch(Point *p[][4], float err2) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				this->p[i][j] = p[i][j];
			}
		}
		adaptiveSubdivide(err2);
	};

	void uniformSubdivide(int subdivision, float step) {
		Point dummy(0, 0, 0);
		c = new Point**[subdivision + 1];
		n = new Point**[subdivision + 1];
		float u, v;
		for (int i = 0; i <= subdivision; i++) {
			u = i * step;
			c[i] = new Point*[subdivision + 1];
			n[i] = new Point*[subdivision + 1];
			Point* a = dummy.merge(p[0][0], p[1][0], p[2][0], p[3][0], u)[0];
			Point* b = dummy.merge(p[0][1], p[1][1], p[2][1], p[3][1], u)[0];
			Point* c = dummy.merge(p[0][2], p[1][2], p[2][2], p[3][2], u)[0];
			Point* d = dummy.merge(p[0][3], p[1][3], p[2][3], p[3][3], u)[0];
			for (int j = 0; j <= subdivision; j++) {
				v = j * step;
				Point ** cn = dummy.merge(a, b, c, d, v);
				this->c[i][j] = cn[0];
				this->n[i][j] = cn[1];
			}
		}

		for (int i = 0; i <= subdivision; i++) {
			v = i * step;
			Point* a = dummy.merge(p[0][0], p[0][1], p[0][2], p[0][3], v)[0];
			Point* b = dummy.merge(p[1][0], p[1][1], p[1][2], p[1][3], v)[0];
			Point* c = dummy.merge(p[2][0], p[2][1], p[2][2], p[2][3], v)[0];
			Point* d = dummy.merge(p[3][0], p[3][1], p[3][2], p[3][3], v)[0];
			for (int j = 0; j <= subdivision; j++) {
				u = j * step;
				Point ** cn = dummy.merge(a, b, c, d, u);
				this->n[j][i] = this->n[j][i]->cross(cn[1], false);
			}
		}
	}

	Point** compute(float u, float v) {
		Point **r1, **r2;

		Point dummy(0, 0, 0);
		Point* a = dummy.merge(p[0][0], p[1][0], p[2][0], p[3][0], u)[0];
		Point* b = dummy.merge(p[0][1], p[1][1], p[2][1], p[3][1], u)[0];
		Point* c = dummy.merge(p[0][2], p[1][2], p[2][2], p[3][2], u)[0];
		Point* d = dummy.merge(p[0][3], p[1][3], p[2][3], p[3][3], u)[0];
		r1 = dummy.merge(a, b, c, d, v);

		a = dummy.merge(p[0][0], p[0][1], p[0][2], p[0][3], v)[0];
		b = dummy.merge(p[1][0], p[1][1], p[1][2], p[1][3], v)[0];
		c = dummy.merge(p[2][0], p[2][1], p[2][2], p[2][3], v)[0];
		d = dummy.merge(p[3][0], p[3][1], p[3][2], p[3][3], v)[0];
		r2 = dummy.merge(a, b, c, d, u);

		/*float e = .00001;
		if (sqr(r1[0]->x - r2[0]->x) > e || sqr(r1[0]->y - r2[0]->y) > e || sqr(r1[0]->z - r2[0]->z) > e) {
		printf("error: ");
		printf("%f=%f, ", r1[0]->x, r2[0]->x);
		printf("%f=%f, ", r1[0]->y, r2[0]->y);
		printf("%f=%f\n", r1[0]->z, r2[0]->z);
		}*/

		r1[1] = r1[1]->cross(r2[1], false);

		return r1;
	};

	void addTri(float u1, float v1, float u2, float v2, float u3, float v3, float err2) {
		float u12 = (u1 + u2) / 2;
		float v12 = (v1 + v2) / 2;
		float u13 = (u1 + u3) / 2;
		float v13 = (v1 + v3) / 2;
		float u23 = (u2 + u3) / 2;
		float v23 = (v2 + v3) / 2;

		Point** p1 = compute(u1, v1);
		Point** p2 = compute(u2, v2);
		Point** p3 = compute(u3, v3);

		Point** r1 = compute(u23, v23);
		Point* w1 = p2[0]->merge(p3[0], .5f);
		float e1 = r1[0]->distance2(w1);

		Point** r2 = compute(u13, v13);
		Point* w2 = p1[0]->merge(p3[0], .5f);
		float e2 = r2[0]->distance2(w2);

		Point** r3 = compute(u12, v12);
		Point* w3 = p1[0]->merge(p2[0], .5f);
		float e3 = r3[0]->distance2(w3);

		addTri(u1, v1, u2, v2, u3, v3, err2, p1, p2, p3, r1, r2, r3, w1, w2, w3, e1, e2, e3);
	};

	void addTri(float u1, float v1, float u2, float v2, float u3, float v3, float err2, Point** p1, Point** p2, Point** p3, Point** r1, Point** r2, Point** r3, Point* w1, Point* w2, Point* w3, float e1, float e2, float e3) {
		float u12 = (u1 + u2) / 2;
		float v12 = (v1 + v2) / 2;
		float u13 = (u1 + u3) / 2;
		float v13 = (v1 + v3) / 2;
		float u23 = (u2 + u3) / 2;
		float v23 = (v2 + v3) / 2;

		if (r1 == NULL) {
			r1 = compute(u23, v23);
			w1 = p2[0]->merge(p3[0], .5f);
			e1 = r1[0]->distance2(w1);
		}

		if (r2 == NULL) {
			r2 = compute(u13, v13);
			w2 = p1[0]->merge(p3[0], .5f);
			e2 = r2[0]->distance2(w2);
		}

		if (r3 == NULL) {
			r3 = compute(u12, v12);
			w3 = p1[0]->merge(p2[0], .5f);
			e3 = r3[0]->distance2(w3);
		}

		int er = 0;
		if (e1 > err2)
			er += 1;
		if (e2 > err2)
			er += 2;
		if (e3 > err2)
			er += 4;

		switch (er) {
		case 0: // no errors
			{ 
				Point*** tri = new Point**[3];
				tri[0] = p1;
				tri[1] = p2;
				tri[2] = p3;
				this->tri.push_back(tri);
			}
			break;
		case 1: // error 1
			addTri(u1, v1, u2, v2, u23, v23, err2, p1, p2, r1, NULL, NULL, r3, NULL, NULL, w3, 0, 0, e3);
			addTri(u23, v23, u3, v3, u1, v1, err2, r1, p3, p1, r2, NULL, NULL, w2, NULL, NULL, e2, 0, 0);
			break;
		case 2: // error 2
			addTri(u1, v1, u2, v2, u13, v13, err2, p1, p2, r2, NULL, NULL, r3, NULL, NULL, w3, 0, 0, e3);
			addTri(u13, v13, u2, v2, u3, v3, err2, r2, p2, p3, r1, NULL, NULL, w1, NULL, NULL, e1, 0, 0);
			break;
		case 3: // error 1 and 2
			addTri(u13, v13, u1, v1, u2, v2, err2, r2, p1, p2, r3, NULL, NULL, w3, NULL, NULL, e3, 0, 0);
			addTri(u3, v3, u13, v13, u23, v23, err2, p3, r2, r1, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u23, v23, u13, v13, u2, v2, err2, r1, r2, p2, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			break;
		case 4: // error 3
			addTri(u3, v3, u1, v1, u12, v12, err2, p3, p1, r3, NULL, NULL, r2, NULL, NULL, w2, 0, 0, e2);
			addTri(u12, v12, u2, v2, u3, v3, err2, r3, p2, p3, r1, NULL, NULL, w1, NULL, NULL, e1, 0, 0);
			break;
		case 5: // error 1 and 3
			addTri(u1, v1, u12, v12, u23, v23, err2, p1, r3, r1, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u12, v12, u2, v2, u23, v23, err2, r3, p2, r1, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u23, v23, u3, v3, u1, v1, err2, r1, p3, p1, r2, NULL, NULL, w2, NULL, NULL, e2, 0, 0);
			break;
		case 6: // error 2 and 3
			addTri(u13, v13, u1, v1, u12, v12, err2, r2, p1, r3, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u12, v12, u2, v2, u3, v3, err2, r3, p2, p3, r1, NULL, NULL, w1, NULL, NULL, e1, 0, 0);
			addTri(u3, v3, u13, v13, u12, v12, err2, p3, r2, r3, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			break;
		case 7: // error 1, 2, and 3
			addTri(u1, v1, u12, v12, u13, v13, err2, p1, r3, r2, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u2, v2, u23, v23, u12, v12, err2, p2, r1, r3, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u3, v3, u13, v13, u23, v23, err2, p3, r2, r1, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			addTri(u12, v12, u23, v23, u13, v13, err2, r3, r1, r2, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0);
			break;
		}
	};

	// this is a (very slightly) slower verision of addTri, but less likely to have typo's
	void addTriSafe(float u1, float v1, float u2, float v2, float u3, float v3, float err2) {
		int er = 0;

		float u12 = (u1 + u2) / 2;
		float v12 = (v1 + v2) / 2;
		float u13 = (u1 + u3) / 2;
		float v13 = (v1 + v3) / 2;
		float u23 = (u2 + u3) / 2;
		float v23 = (v2 + v3) / 2;

		Point** p1 = compute(u1, v1);
		Point** p2 = compute(u2, v2);
		Point** p3 = compute(u3, v3);

		Point** r1 = compute(u23, v23);
		Point* w1 = p2[0]->merge(p3[0], .5f);
		float e1 = r1[0]->distance2(w1);
		if (e1 > err2)
			er += 1;

		Point** r2 = compute(u13, v13);
		Point* w2 = p1[0]->merge(p3[0], .5f);
		float e2 = r2[0]->distance2(w2);
		if (e2 > err2)
			er += 2;

		Point** r3 = compute(u12, v12);
		Point* w3 = p1[0]->merge(p2[0], .5f);
		float e3 = r3[0]->distance2(w3);
		if (e3 > err2)
			er += 4;

		switch (er) {
		case 0: // no errors
			{ 
				Point*** tri = new Point**[3];
				tri[0] = p1;
				tri[1] = p2;
				tri[2] = p3;
				this->tri.push_back(tri);
			}
			break;
		case 1: // error 1
			addTri(u1, v1, u2, v2, u23, v23, err2);
			addTri(u23, v23, u3, v3, u1, v1, err2);
			break;
		case 2: // error 2
			addTri(u1, v1, u2, v2, u13, v13, err2);
			addTri(u13, v13, u2, v2, u3, v3, err2);
			break;
		case 3: // error 1 and 2
			addTri(u13, v13, u1, v1, u2, v2, err2);
			addTri(u3, v3, u13, v13, u23, v23, err2);
			addTri(u23, v23, u13, v13, u2, v2, err2);
			break;
		case 4: // error 3
			addTri(u3, v3, u1, v1, u12, v12, err2);
			addTri(u12, v12, u2, v2, u3, v3, err2);
			break;
		case 5: // error 1 and 3
			addTri(u1, v1, u12, v12, u23, v23, err2);
			addTri(u12, v12, u2, v2, u23, v23, err2);
			addTri(u23, v23, u3, v3, u1, v1, err2);
			break;
		case 6: // error 2 and 3
			addTri(u13, v13, u1, v1, u12, v12, err2);
			addTri(u12, v12, u2, v2, u3, v3, err2);
			addTri(u3, v3, u13, v13, u12, v12, err2);
			break;
		case 7: // error 1, 2, and 3
			addTri(u1, v1, u12, v12, u13, v13, err2);
			addTri(u2, v2, u23, v23, u12, v12, err2);
			addTri(u3, v3, u13, v13, u23, v23, err2);
			addTri(u12, v12, u23, v23, u13, v13, err2);
			break;
		}
	};

	void adaptiveSubdivide(float err2) {
		addTri(0, 0, 1, 1, 0, 1, err2);
		addTri(0, 0, 1, 0, 1, 1, err2);

		// alternatively:
		// addTri[Safe](1, 0, 1, 1, 0, 1, err2);
		// addTri[Safe](1, 0, 0, 1, 0, 0, err2);
	};
};

class Model {
public:
	float x, y, z;
	float theta, phi;
	float s;
};

class Camera {
public:
	float ox , oy, oz;
	float tx, ty, tz;
	float r, flatr;
	float theta, thetac, thetas;
	float phi, phic, phis;

	bool wireframe, smooth, control;

	Camera(float tx, float ty, float tz, float ox, float oy, float oz){
		this->tx = tx;
		this->ty = ty;
		this->tz = tz;
		this->ox = ox;
		this->oy = oy;
		this->oz = oz;
		float dx = ox - tx;
		float dy = oy - ty;
		float dz = oz - tz;
		r = sqrt(dx*dx + dy*dy + dz*dz);
		flatr = sqrt(dx*dx + dz*dz);
		theta = atan2(dx, dz);
		thetac = cos(theta);
		thetas = sin(theta);
		phi = atan2(dy, flatr);
		phic = cos(phi);
		phis = sin(phi);
	}
};

Viewport viewport;
char *inputFile;
int subdivision;
float step, err2; 
bool adaptive;
int patchCount;
Patch **patch;
Camera camera(0, 0, 0, 10, 10, 10);
Model model;
int mouse, mx, my;
int displayListIndex;

void initScene(){
	//This tells glut to use a double-buffered window with red, green, and blue channels 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_DEPTH

	// Initalize theviewport size
	viewport.w = 800;
	viewport.h = 800;

	// The size and position of the window
	glutInitWindowSize(viewport.w, viewport.h);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Bezoar Stones");

	glEnable(GL_DEPTH_TEST);

	// lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHT0);
	GLfloat pos0[] = {10000, 5000, 4000, 1.0f};
	GLfloat amb0[] = {.0f, .0f, .0f, 1.0f};
	GLfloat dif0[] = {1, 1, 1, 1.0f};
	GLfloat spec0[] = {.0f, .0f, .0f, 1.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, pos0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, amb0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, dif0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, spec0);

	glEnable(GL_LIGHT1);
	GLfloat pos1[] = {-30000, 100000, 40000, 1.0f};
	GLfloat amb1[] = {.0f, .2f, .2f, 1.0f};
	GLfloat dif1[] = {.4f, .8f, .4f, 1.0f};
	GLfloat spec1[] = {.0f, .8f, .0f, 1.0f};
	glLightfv(GL_LIGHT1, GL_POSITION, pos1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, amb1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, dif1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, spec1);

	// default flat shading
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);
}

void myReshape(int w, int h) {
	int s = min(w, h);
	
	viewport.w = s;
	viewport.h = s;

	glViewport(0, 0, viewport.w, viewport.h);
}

void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// camera
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, 0.1, 100);
	gluLookAt(camera.ox, camera.oy, camera.oz, // look from
		camera.tx, camera.ty, camera.tz, // look at
		0, 1 , 0); // up axis

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(model.x, model.y, model.z);
	glRotatef(model.theta, 0.0, 1.0, 0.0 );
	glRotatef(model.phi, 1, 0.0, 0);

	// Start drawing
	if (adaptive)
		drawAdaptive();
	else
		drawUniform();
	if (camera.control)
		drawControl();

	glFlush();
	glutSwapBuffers();
}

void myFrameMove() {
	//nothing here for now
#ifdef _WIN32
	Sleep(10); //give ~10ms back to OS (so as not to waste the CPU)
#endif
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

void drawControl() {
	glDisable(GL_LIGHTING);

	glColor3f(0.0, 0.0, 1.0);

	for (int i = 0; i < patchCount; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				glBegin(GL_LINE_LOOP);

				glVertex3f(patch[i]->p[j][k]->x, patch[i]->p[j][k]->y, patch[i]->p[j][k]->z);
				glVertex3f(patch[i]->p[j+1][k]->x, patch[i]->p[j+1][k]->y, patch[i]->p[j+1][k]->z);
				glVertex3f(patch[i]->p[j+1][k+1]->x, patch[i]->p[j+1][k+1]->y, patch[i]->p[j+1][k+1]->z);
				glVertex3f(patch[i]->p[j][k+1]->x, patch[i]->p[j][k+1]->y, patch[i]->p[j][k+1]->z);

				glEnd();
			}
		}
	}

}

void drawAdaptive() {	
	if (displayListIndex == 0) {
		displayListIndex = glGenLists(1);
		glNewList(displayListIndex, GL_COMPILE_AND_EXECUTE);

		glEnable(GL_LIGHTING);

		glBegin(GL_TRIANGLES);

		for (int i = 0; i < patchCount; i++) {
			for (int j = 0; j < patch[i]->tri.size(); j++) {
				Point*** tri = patch[i]->tri.at(j);

				glNormal3f(tri[0][1]->x, tri[0][1]->y, tri[0][1]->z);
				glVertex3f(tri[0][0]->x, tri[0][0]->y, tri[0][0]->z);

				glNormal3f(tri[1][1]->x, tri[1][1]->y, tri[1][1]->z);
				glVertex3f(tri[1][0]->x, tri[1][0]->y, tri[1][0]->z);

				glNormal3f(tri[2][1]->x, tri[2][1]->y, tri[2][1]->z);
				glVertex3f(tri[2][0]->x, tri[2][0]->y, tri[2][0]->z);
			}
		}

		glEnd();

		glEndList();
	} else
		glCallList(displayListIndex);
}

void drawUniform() {	
	if (displayListIndex == 0) {
		displayListIndex = glGenLists(1);
		glNewList(displayListIndex, GL_COMPILE_AND_EXECUTE);

		glEnable(GL_LIGHTING);

		glBegin(GL_QUADS);

		for (int i = 0; i < patchCount; i++) {
			for (int j = 0; j < subdivision; j++) {
				for (int k = 0; k < subdivision; k++) {
					glNormal3f(patch[i]->n[j][k]->x, patch[i]->n[j][k]->y, patch[i]->n[j][k]->z);
					glVertex3f(patch[i]->c[j][k]->x, patch[i]->c[j][k]->y, patch[i]->c[j][k]->z);

					glNormal3f(patch[i]->n[j+1][k]->x, patch[i]->n[j+1][k]->y, patch[i]->n[j+1][k]->z);
					glVertex3f(patch[i]->c[j+1][k]->x, patch[i]->c[j+1][k]->y, patch[i]->c[j+1][k]->z);

					glNormal3f(patch[i]->n[j+1][k+1]->x, patch[i]->n[j+1][k+1]->y, patch[i]->n[j+1][k+1]->z);
					glVertex3f(patch[i]->c[j+1][k+1]->x, patch[i]->c[j+1][k+1]->y, patch[i]->c[j+1][k+1]->z);

					glNormal3f(patch[i]->n[j][k+1]->x, patch[i]->n[j][k+1]->y, patch[i]->n[j][k+1]->z);
					glVertex3f(patch[i]->c[j][k+1]->x, patch[i]->c[j][k+1]->y, patch[i]->c[j][k+1]->z);
				}
			}
		}

		glEnd();

		glEndList();
	} else
		glCallList(displayListIndex);
}

void initCommands(int argc, char *argv[]) {
	if (argc == 4) {		
		inputFile = argv[1];
		err2 = atof(argv[2]);
		subdivision = (int) ceilf(1 / err2);
		step = 1.0f / subdivision;
		err2 *= err2;
		adaptive = (strcmp(argv[3], "-a") == 0);
		printf("file: %s, subdivision: %d, step: %f, adaptive: %d \n", inputFile, subdivision, step, adaptive);
	} else if (argc == 3) {
		inputFile = argv[1];
		subdivision = (int) ceilf(1 / atof(argv[2]));
		step = 1.0f / subdivision;
		adaptive = false;
		printf("file: %s, subdivision: %d, step: %f, adaptive: %d \n", inputFile, subdivision, step, adaptive);
	} else {
		printf("bad parmameters\n");
		exit(0);
	}
}

void readInputFile() {
	//long t1 = GetTickCount();

	string data;
	Point *p[4][4];

	ifstream file;
	file.open (inputFile);
	file >> data;

	patchCount = atof(data.c_str());
	patch = new Patch*[patchCount];

	for (int i  = 0; i < patchCount; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				file >> data;
				float x = atof(data.c_str());
				file >> data;
				float y = atof(data.c_str());
				file >> data;
				float z = atof(data.c_str());
				p[j][k] = new Point(x, y, z);
			}
		}
		if (adaptive)
			patch[i] = new Patch(p, err2);
		else
			patch[i] = new Patch(p, subdivision, step);
	}

	//long t2 = GetTickCount();

	printf("patches: %d", patchCount);

	//printf("time: %d\n", t2 - t1);

	file.close();
}

void myInput(unsigned char key, int x, int y) {
	if (key == ' ')
		exit(0);
	switch(key) {
	case ' ':
		exit(0);
		break;
	case 'w':
		camera.wireframe = !camera.wireframe;
		if (camera.wireframe)
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		else
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case 's':
		camera.smooth = !camera.smooth;
		if (camera.smooth) {
			glShadeModel(GL_SMOOTH);
			glPolygonMode(GL_FRONT_AND_BACK, GL_SMOOTH);
		} else {
			glShadeModel(GL_FLAT);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FLAT);
		}
		break;
	case 'c':
		camera.control = !camera.control;
		break;
	case '=':
	case '+':
		model.z += .1f;
		break;
	case '-':
		model.z -= .1f;
		break;
	case 'z':
	case 'Z':
		model.phi = -90;
		model.theta = 0;
		break;
	case 'y':
	case 'Y':
		model.phi = 0;
		model.theta = 0;
		break;
	}
}

void mySpecialInput(int key, int x, int y) {
	bool shift = glutGetModifiers() == GLUT_ACTIVE_SHIFT;

	if (!shift) {
		switch(key)	{
		case GLUT_KEY_UP:
			model.phi -= 4;
			break;	
		case GLUT_KEY_DOWN:
			model.phi += 4;
			break;
		case GLUT_KEY_LEFT:
			model.theta -= 4;
			break;
		case GLUT_KEY_RIGHT:
			model.theta += 4;
			break;
		}
		printf("phi %f\n", model.phi);
	} else {
		switch(key)	{
		case GLUT_KEY_UP:
			model.y += .1f;
			break;	
		case GLUT_KEY_DOWN:
			model.y -= .1f;
			break;
		case GLUT_KEY_LEFT:
			model.x -= .1f;
			break;
		case GLUT_KEY_RIGHT:
			model.x += .1f;
			break;
		}
	}
}

void myMouseInput(int button, int state, int x, int y) {
	if (state == 1)
		mouse = -1;
	else {
		mouse = button;
		mx = x;
		my = y;
	}
}

void myMotionInput(int x, int y) {
	float dy = (y - my) / 10.;
	float dx = (x - mx) / 10.;
	my = y;
	mx = x;

	switch (mouse) {
	case 0:
		camera.theta += dx / 10.;
		camera.phi += dy / 10.;
		if (camera.phi > PI / 2 - .15)
			camera.phi = PI / 2 - .15;
		if (camera.phi < -PI / 2 + .15)
			camera.phi = -PI / 2 + .15;
		camera.thetac = cos(camera.theta);
		camera.thetas = sin(camera.theta);
		camera.phic = cos(camera.phi);
		camera.phis = sin(camera.phi);
		camera.flatr = camera.r * camera.phic;
		camera.ox = camera.tx + camera.flatr * camera.thetac;
		camera.oy = camera.ty + camera.r * camera.phis;
		camera.oz = camera.tz + camera.flatr * camera.thetas;
		break;
	case 1:
		camera.r += dy;
		camera.flatr = camera.r * camera.phic;
		camera.ox = camera.tx + camera.flatr * camera.thetac;
		camera.oy = camera.ty + camera.r * camera.phis;
		camera.oz = camera.tz + camera.flatr * camera.thetas;
		break;
	case 2:
		float moveh = dy * camera.phis;
		float movev = dy * camera.phic;
		camera.ty += movev;
		camera.oy += movev;
		float movex = dx * camera.thetas + moveh * camera.thetac;
		camera.tx -= movex;
		camera.ox -= movex;
		float movez = dx * camera.thetac - moveh * camera.thetas;
		camera.tz += movez;
		camera.oz += movez;
		break;
	}
}

int main(int argc, char *argv[]) {
	initCommands(argc, argv);
	readInputFile();

	glutInit(&argc, argv);

	initScene(); // quick function to set up scene

	glutDisplayFunc(myDisplay); // function to run when its time to draw something
	glutReshapeFunc(myReshape); // function to run when the window gets resized
	glutIdleFunc(myFrameMove); // function to run when not handling any other task
	glutKeyboardFunc(myInput);
	glutSpecialFunc(mySpecialInput);
	glutMouseFunc(myMouseInput);
	glutMotionFunc(myMotionInput);
	glutMainLoop(); // infinite loop that will keep drawing and resizing and whatever else

	return 0;
}