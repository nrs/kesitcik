
#ifndef __PRGL_H__ 
#define __PRGL_H__

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>         /* glut.h includes gl.h and glu.h*/
#endif

#include "pr1.h"
#include "glutMaster.h"
#include <iostream>

void resize(int width, int height);

void SetColor(int i);



class MeshWindow : public GlutWindow{
public:

  int height, width;
  int initPositionX, initPositionY;

  MeshWindow(GlutMaster * glutMaster,
             int setWidth, int setHeight,
             int setInitPositionX, int setInitPositionY,
             char * title);
  ~MeshWindow();

  void CallBackDisplayFunc(void);
  void CallBackMotionFunc(int x, int y);
  void CallBackMouseFunc(int button, int state, int x, int y);
  void CallBackKeyboardFunc(unsigned char key, int x, int y);
  void init();

//  mesh *glmesh;
//  std::vector<mesh*> glmesh;
  supermesh *glsupermesh;

//  void add_mesh(std::list<mesh> &m);
  void draw_supermesh();
  void draw_mesh(mesh &m);
  void draw_node(node &n, Real rad);
  void draw_node2(node &n, Real rad);
  
  void draw_line(line &l, int c);
  void draw_crosshair(node &n, Real len, int c);

  void draw_triangle(triangle &t, int c);
//   void CallBackReshapeFunc(int w, int h);   
//   void CallBackIdleFunc(void);
};




#endif
// Local Variables:
// mode: c++
// End:
