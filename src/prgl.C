#include "prgl.h"
//#include "zpr.h"


namespace prgl{
  mesh *glmesh;
};

#define VIEWING_DISTANCE_MIN  5.0

//static bool g_bLightingEnabled = true;
//static bool g_bFillPolygons = true;
//static bool g_bTexture = false;
static bool g_bButton1Down = false;
//static GLfloat g_fTeapotAngle = 0.0;
//static GLfloat g_fTeapotAngle2 = 0.0;
static GLfloat g_scaling_factor = 0.1;
static GLfloat g_fViewDistance = 3 * VIEWING_DISTANCE_MIN;
static GLfloat g_nearPlane = 10;
static GLfloat g_farPlane = 10000;
static int g_Width = 600;                          // Initial window width
static int g_Height = 600;                         // Initial window height
static int g_yClick = 0;
//static float g_lightPos[4] = { 10, 10, -100, 1 };  // Position of light

static GLfloat basic_colors[]=
{
    0.0, 0.0, 0.0 ,    1.0, 0.0, 0.0 ,     0.0, 1.0, 0.0 ,
    0.0, 0.0, 1.0 ,    1.0, 1.0, 0.0 ,     0.0, 1.0, 1.0 ,
    1.0, 0.0, 1.0 ,    1.0, 1.0, 1.0
};

void SetColor(int i)
{
    glColor3fv(basic_colors+((i%8)*3));
}


// void resize(int width, int height)
// {
//     const double ar = (double) width / (double) height;

//     glViewport(0, 0, width, height);
//     glMatrixMode(GL_PROJECTION);
//     glLoadIdentity();
//     glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);

//     glMatrixMode(GL_MODELVIEW);
//     glLoadIdentity() ;
// }

// void  glMeshRenderer::reshape(GLint width, GLint height)
// {
//    g_Width = width;
//    g_Height = height;

//    glViewport(0, 0, g_Width, g_Height);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(65.0, (float)g_Width / g_Height, g_nearPlane, g_farPlane);
//    glMatrixMode(GL_MODELVIEW);
// }




void  MeshWindow::CallBackMouseFunc(int button, int state, int x, int y)
{
  // Respond to mouse button presses.
  // If button1 pressed, mark this state so we know in motion function.

  if (button == GLUT_LEFT_BUTTON)
    {
      g_bButton1Down = (state == GLUT_DOWN) ? true : false;
      g_yClick = y - 3 * g_fViewDistance;
    }
}
void MeshWindow::CallBackMotionFunc(int x, int y)
{
  // If button1 pressed, zoom in/out if mouse is moved up/down.

  if (g_bButton1Down)
    {
      g_fViewDistance = (y - g_yClick) / 3.0;
      if (g_fViewDistance < VIEWING_DISTANCE_MIN)
         g_fViewDistance = VIEWING_DISTANCE_MIN;
      glutPostRedisplay();
    }
}

void MeshWindow::draw_node2(node &n, Real rad, int c){
    glPushMatrix();
//    SetColor(c);
    SetColor(n.color);

    glTranslatef(n(0),n(1),0);

    glBegin(GL_QUADS);

    glVertex2f(-1*rad, rad);    // A
    glVertex2f(rad,rad);    // B
    glVertex2f( rad, -1* rad);    // C
    glVertex2f( -1*rad,  -1*rad);

    glEnd();
    glPopMatrix();
}


void MeshWindow::draw_node(node &n, Real rad, int c){

    Real th,x,y;
    glPushMatrix();
//    SetColor(c);
    SetColor(n.color);

    glTranslatef(n(0),n(1),0);
    glBegin(GL_TRIANGLE_FAN);

    for(int j=0; j<=360; j++)
    {
        th=M_PI * j / 180.0;
        x = rad * cos(th);
        y= rad * sin(th);
        glVertex2f(x, y);
    }

    glEnd();
    glPopMatrix();

}

void MeshWindow::draw_line(line &l, int c)
{
    glPushMatrix();
//    glTranslatef((*l.nodes[1])(0),(*l.nodes[1])(1),0);
//    SetColor(c);
    SetColor(l.color);

    glBegin(GL_LINES);
//    glVertex2f(0.0, 0.0);
    glVertex2f((*l.nodes[0])(0),(*l.nodes[0])(1));
    glVertex2f((*l.nodes[1])(0),(*l.nodes[1])(1));
    glEnd();
    glPopMatrix();

}


void MeshWindow::draw_triangle(triangle &t, int c)
{
    glPushMatrix();
    SetColor(t.color);

//    SetColor(c);
    glBegin(GL_TRIANGLES);
    node *n0, *n2, *n1;
    n0=t.nodes[0];
    n1=t.nodes[1];
    n2=t.nodes[2];
    glVertex2f((*n0)(0),(*n0)(1));
    glVertex2f((*n1)(0),(*n1)(1));
    glVertex2f((*n2)(0),(*n2)(1));



    glEnd();
    glEnd();
    glPopMatrix();

}

// void MeshWindow::add_mesh(std::list<mesh> &m)
// {
//   std::list<mesh>::iterator mesh_it;
//   for (mesh_it = m.begin(); mesh_it != m.end(); mesh_it++){
//     glmesh.push_back(&(*mesh_it));
//   }
// }
void MeshWindow::draw_mesh(mesh &m)
{
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;

  std::list<triangle>::iterator tri_it;



//  m.print_info();

    for (tri_it=m.triangles.begin(); 
         tri_it!=m.triangles.end(); tri_it++){
      draw_triangle(*tri_it,4);

    }

    // for (node_it=m.nodes.begin(); node_it!=m.nodes.end(); node_it++){
    //   draw_node2(*node_it,0.5,0);
    // }

    // for (node_it=(*mesh_it)->centroids.begin(); 
    //      node_it!=(*mesh_it)->centroids.end(); node_it++){
    //   draw_node2(*node_it,0.08,0);
    // }


    for (line_it= m.lines.begin(); line_it!= m.lines.end(); line_it++){
      draw_line(*line_it,0);
    }

    // for (bine_it=(*mesh_it)->hull.begin(); bine_it!=(*mesh_it)->hull.end(); bine_it++){
    //   draw_line(*(*bine_it),1);
    // }

  

}
void MeshWindow::draw_supermesh()
{
//  float a = ((*glsupermesh->meshes.front().triangles.front().nodes[2])(0));
//  cout << &((*glsupermesh->meshes.front().triangles.front().nodes[2])(0)) <<endl;
//  cout << a << endl;
//  printf("%f\n",(*glsupermesh->meshes.front().triangles.front().nodes[2])(0) );

  std::list<mesh>::iterator mesh_it;
  std::list<line>::iterator line_it;

  for (mesh_it = glsupermesh->meshes.begin(); 
       mesh_it!= glsupermesh->meshes.end(); mesh_it++){
    if (mesh_it->nodes.size()==0) continue;
    draw_mesh(*mesh_it);
  }

  for (line_it = glsupermesh->misclines.begin(); 
       line_it!= glsupermesh->misclines.end(); line_it++){
    draw_line(*line_it,0);
  }

}



void MeshWindow::CallBackDisplayFunc()
{
/* clear window */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();

  SetColor(0);
  
//  glOrtho(0, 400,  600, 0,  -1.0, -1.0);
  glScaled(g_scaling_factor/g_fViewDistance,
           g_scaling_factor/g_fViewDistance,
           g_scaling_factor/g_fViewDistance);
//  std::cout << g_fViewDistance << std::endl;
//  glLoadIdentity();

//  gluLookAt(0, 0, -g_fViewDistance/10000000, 0, 0, -1, 0, 1, 0);
//  gluLookAt(0, 0, -9, 0, 0, -1, 0, 1, 0);

  draw_supermesh();


/* draw unit square polygon */
/*
  glBegin(GL_POLYGON);
  glVertex2f(-0.5, -0.5);
  glVertex2f(-0.5, 0.5);
  glVertex2f(0.5, 0.5);
  glVertex2f(0.5, -0.5);
  glEnd();
*/

/* flush GL buffers */
   glutSwapBuffers();

//  glFlush();
}
// void  glMeshRenderer::idle(void)
// {
 
//     //timex = glutGet(GLUT_ELAPSED_TIME) ;

  
// //  frame++;
// //  time2=glutGet(GLUT_ELAPSED_TIME);
// //  std::cout << "ASDASDASD";
// //  if (time2 - timebase > 15) {
//     glutPostRedisplay();
//     //fps = frame*1000.0/(time2-timebase);
// //    timebase = time2;
// //    frame = 0;
// //  }
  
// //    Sleep(10);
// }

void  MeshWindow::init()
{


/* set clear color to black */
  glClearColor (1, 1, 1, 1);

/* set fill color to white */
  glColor3f(0,0,0);

/* set up standard orthogonal view with clipping */
/* box as cube of side 2 centered at origin */
/* This is default view and these statement could be removed */
   glMatrixMode (GL_PROJECTION);
   glLoadIdentity ();


//  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
}




void MeshWindow::CallBackKeyboardFunc(unsigned char key, int x, int y)
{

    //usercontrol(key, &bodyarray[0]);

    switch(key)
    {

    case 'g':
//      glmesh->generate1();
      glutPostRedisplay();
      break;

    case 's':
//      glmesh->out_vtk1();
      break;

    case '+':
      g_scaling_factor*=2;
      glutPostRedisplay();
      break;
    case '-':
      g_scaling_factor*=0.5;
      glutPostRedisplay();
      break;
    case 'n':
      std::cout << std::endl << "Enter the number of nodes to  be generated:"<<
        std::endl<< ">>> ";
//      std::cin >> glmesh->nnodes;
//      std::cout << std::endl;
//      glmesh->generate1();
      glutPostRedisplay();
      break;

    default:

      break;
    }
}

MeshWindow::MeshWindow(GlutMaster * glutMaster,
                       int setWidth, int setHeight,
                       int setInitPositionX, int setInitPositionY,
                       char * title){

   width  = setWidth;               
   height = setHeight;

   initPositionX = setInitPositionX;
   initPositionY = setInitPositionY;

   glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(width, height);
   glutInitWindowPosition(initPositionX, initPositionY);
   glViewport(0, 0, width, height);   // This may have to be moved to after the next line on some platforms

   glutMaster->CallGlutCreateWindow(title, this);

//   glEnable(GL_DEPTH_TEST);

//   glMatrixMode(GL_PROJECTION);
//   glOrtho(-80.0, 80.0, -80.0, 80.0, -500.0, 500.0);

//   glMatrixMode(GL_MODELVIEW);
//   glLoadIdentity();

//   glRotatef(60, 1, 1, 1);
//   glColor4f(1.0, 0.0, 0.0, 1.0);
}

MeshWindow::~MeshWindow(){

   glutDestroyWindow(windowID);
}

