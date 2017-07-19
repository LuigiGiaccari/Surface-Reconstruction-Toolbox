//Classe che contiene tutte le funzioni per disegnare un viewr di una superficie triangolarizzata con le opengl

namespace VIEWER
{
using namespace VIEWER;

#include "stdio.h"

#include "windows.h"
#include "gl/gl.h"
#include "gl/glu.h"



// freeglut library
#define FREEGLUT_STATIC
#include "gl/glut.h"
#include "math.h"
#include "util/ArraysLib.h"//for memory allocation and arrys operation


#ifndef M_PI
#define M_PI  3.14159265358979323846f
#endif
#define DEGREES_PER_PIXEL 0.6f
#define RADIANS_PER_PIXEL 0.002f
#define UNITS_PER_PIXEL 0.0025f
#define UNITS_PER_WHEELTICK 0.35f
#define ZOOM_FACTOR .04f



//For the mouse state
struct g_mouseState
{
    bool leftButton;
    bool rightButton;
    bool middleButton;
    int x;
    int y;
} MouseState;


//functions


//Global variables
float eyex, eyey, eyez;  // eye x,y,z values for gluLookAt (location of eye)
float focusx, focusy, focusz; // the point the eye is looking at
GLdouble  CRotX, CRotY, CRotZ;
GLdouble  BarX,  BarY,  BarZ;
double FactScale;
double Diam;
float g_xRotation, g_yRotation;
float lmodel_ambient[4]= {0.1f,0.1f,0.1f,0.1f};
GLfloat * buffer;
int Nfacets;


//FUNCTIONS DECLARATION
void InitRenderer();
void DrawCube();
void InitializeLight();
void HandleMouseState(int button, int state, int x, int y);
void HandleMouseMove(int x, int y);
void keyboard (unsigned char key, int x, int y);
void ChangeSize(int w, int h);
void glutSetupRC();
void DrawClickedPoint();
void GetClickedPoint(int x,int y,GLdouble *CRotX, GLdouble *CRotY, GLdouble *CRotZ);
void sistema_di_riferimento(double X,double Y,double Z);
void DrawText(string* Text,void * font,int x,int y, float r,float g, float b);
void PrintMouseState();
void DrawTriangles();
void RenderScene();
void DrawHelpLines();
void GetModelData(double* p,int N);
void Tnorm(double* p1,double* p2,double* p3,float* tnorm);
void FillBuffer(double* p,int N,int *t,int nt);
//"Public"

void StartViewer(double* p,int N,int *t,int nt);






void StartViewer(double* p,int N,int *t,int nt)
{


    Nfacets=nt;
    GetModelData(p,N);//gets model dimensions

    FillBuffer(p,N,t,nt);


    int argc=1;
    char *argv[] = {"mydemo",NULL};


    // creo la finestra
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800,600);
    glutCreateWindow("Viewer");
    glutReshapeFunc(ChangeSize);
    glutMouseFunc(HandleMouseState);
    glutMotionFunc(HandleMouseMove);
    glutKeyboardFunc (keyboard);
    glutDisplayFunc(RenderScene);
    //glutSetupRC();
    InitRenderer();
    glutMainLoop();


    Deallocate(&buffer);
}




void FillBuffer(double* p,int N,int *t,int nt)
{
    int i,j,idt,idb,P1,P2,P3;



    Allocate(&buffer,nt*3*4);//allocate buffer memory

    //loop trough tirangles and fill buffer

    for(i=0; i<nt; i=i++)
    {

        idt=i*3;
        idb=i*12;
        P1=t[idt]*3;
        P2=t[idt+1]*3;
        P3=t[idt+2]*3;

        Tnorm(&p[P1],&p[P2], &p[P3],&buffer[idb]);

        for(j=0; j<3; j++)buffer[idb+j+3]=p[P1+j];

        for(j=0; j<3; j++)buffer[idb+j+6]=p[P2+j];

        for(j=0; j<3; j++)buffer[idb+j+9]=p[P3+j];

    }

}


/*
 our display function, draws a box with proper rotation and
 translation according to our mouse interaction.
*/
void Tnorm(double* p1,double* p2,double* p3,GLfloat* tnorm)
{

    int i;
    double v21[3];
    double v31[3];
    double leng;

    for(i=0; i<3; i++)
    {
        v21[i]=p2[i]-p1[i];
    }
    for(i=0; i<3; i++)
    {
        v31[i]=p3[i]-p1[i];
    }


    tnorm[0]=v21[1]*v31[2]-v21[2]*v31[1];
    tnorm[1]=v21[2]*v31[0]-v21[0]*v31[2];
    tnorm[2]=v21[0]*v31[1]-v21[1]*v31[0];

    leng=sqrt((p1[0]* p1[0])+(p1[1] * p1[1])+(p1[2] * p1[2])) ;

    for(i=0; i<3; i++)
    {
        tnorm[i]=tnorm[i]/leng;
    }
}
//funzione di esempio che disegna un cubo
void DrawCube()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

    // move our eye to the most recent place
    gluLookAt(eyex, eyey, eyez, focusx, focusy, focusz, 0,1,0);

    // rotate our box in the x and y directions
    glRotatef(g_xRotation, 0,1,0);
    glRotatef(g_yRotation, 1,0,0);

    // set our color to something
    glColor3f(.5f, .5f, .5f);
    // draw a cube
    glutSolidCube(1.0);

    glPopMatrix();
    glutSwapBuffers();
}





/*
 add a light to make the box we're showing a bit more interesting
*/

//Questa funzione non la chiamo più mi creava problem con il rendering era tuto buio
void InitializeLight()
{
    /*glMatrixMode(GL_MODELVIEW);

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_CULL_FACE);*/


    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glShadeModel(GL_FLAT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);







}

/*
  some basic initialization of our GL code
*/
void InitRenderer()
{


    FactScale=10/Diam;


    glClearColor(0.0f, 0.0f, 1.0f, 0.0f);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // set up our viewing frustum
    gluPerspective(50.0f, 1.0f, 1, 30);//zfar non mi funziona

    // starting rotation.
    g_yRotation = -45.0f;
    g_xRotation = 30.0f;

    //initialize the mouse state
    MouseState.leftButton = MouseState.rightButton = MouseState.middleButton = false;
    MouseState.x = MouseState.y = 0;

    // init our eye location
    eyex = BarX;
    eyey = BarY;
    eyez = BarZ+Diam*3;
    focusx =BarX;
    focusy =BarY;
    focusz = BarZ;//focus
    CRotX=BarX;
    CRotY=BarY;
    CRotZ=BarZ;//center of rotation

    // create some lighting
    InitializeLight();
}

void PrintMouseState()
{
    printf("mouse x: %d, y:%d left: %s, right:%s, middle:%s\n", MouseState.x,
           MouseState.y, MouseState.leftButton ? "down":"up",
           MouseState.rightButton ? "down":"up", MouseState.middleButton ? "down":"up");
}

/*
  this function is called when any mouse buttons are pressed
*/
void HandleMouseState(int button, int state, int x, int y)
{
    // update our button state
    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
            MouseState.leftButton = true;
        else
            MouseState.leftButton = false;
    }
    if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
            MouseState.rightButton = true;
        else
            MouseState.rightButton = false;
    }
    if(button == GLUT_MIDDLE_BUTTON)
    {
        if(state == GLUT_DOWN)
            MouseState.middleButton = true;
        else
            MouseState.middleButton = false;
    }



    // update our position so we know a delta when the mouse is moved
    MouseState.x = x;
    MouseState.y = y;

    //PrintMouseState(); per evitare di plottare lo stato del mouse
}

void HandleMouseMove(int x, int y)
{
//   double dist;
    //double ClickedX,ClickedY,ClickedZ;
    // calculate a delta in movement
    int yDelta = MouseState.y - y;
    int xDelta = MouseState.x - x;
    // commit the mouse position
    MouseState.x = x;
    MouseState.y = y;

    // when we need to rotate (only the left button is down)
    if(MouseState.leftButton && !MouseState.rightButton && !MouseState.middleButton)
        //ROTATION
    {




        // rotate by the delta
        g_xRotation -= xDelta * DEGREES_PER_PIXEL;
        g_yRotation -= yDelta * DEGREES_PER_PIXEL;

// translate the eye with a value equals to to focus vector (we put the focus on (0,0,0))


    }
    //PAN
    // if we need to move translate (left and right buttons are down
    else if(MouseState.leftButton && MouseState.rightButton && !MouseState.middleButton)
    {

        //  Il Pan avviene nella funzione glulookat che non sembra risentire dell'effetto della scala
        //Dobbiamo considerare la dimensione dell'oggetto (Diam)

        // move our eye
        eyex += xDelta * UNITS_PER_PIXEL*Diam;
        eyey -= yDelta * UNITS_PER_PIXEL*Diam;

        // move our focus point
        focusx += xDelta * UNITS_PER_PIXEL*Diam;
        focusy -= yDelta * UNITS_PER_PIXEL*Diam;

        //Trova il punto cliccato

    }

    //ROTATION CENTER
    else if(!MouseState.leftButton && !MouseState.rightButton && MouseState.middleButton)
    {
        //Trova il punto cliccato

    }
    else if(!MouseState.leftButton && MouseState.rightButton && !MouseState.middleButton)

        //ZOOM
    {
        // zooming.  we move our eye along the vector between the eye
        //  and the focus point.

        //zoom in
        if(yDelta>0)
        {
            eyex = (1 - ZOOM_FACTOR) * eyex + focusx * ZOOM_FACTOR;
            eyey = (1 - ZOOM_FACTOR) * eyey + focusy * ZOOM_FACTOR;
            eyez = (1 - ZOOM_FACTOR) * eyez + focusz * ZOOM_FACTOR;
            //glutPostRedisplay();
        }
        else
        {
            //zoom out;
            eyex = (1 + ZOOM_FACTOR) * eyex - focusx * ZOOM_FACTOR;
            eyey = (1 + ZOOM_FACTOR) * eyey - focusy * ZOOM_FACTOR;
            eyez = (1 + ZOOM_FACTOR) * eyez - focusz * ZOOM_FACTOR;
            //glutPostRedisplay();
        }
    }

    glutPostRedisplay();
    //PrintMouseState(); //PER EVITARE CHE STAMPI DI CONTINUO LA POSIZIONE DEL MOUSE
}

// This function does any needed initialization on the rendering
// context.
void glutSetupRC()
{
    // Bluish background
    glClearColor(0.0f, 0.0f, .50f, 1.0f );



    // Draw everything as wire frame
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}


//come chiamare la finestra nel main
/*
int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(640,480);
	glutCreateWindow("Mouse Rotation");
	glutDisplayFunc(Draw);
	glutMouseFunc(HandleMouseState);
	glutMotionFunc(HandleMouseMove);

	InitRenderer();

	glutMainLoop();

	return 0;
}
*/


void keyboard (unsigned char key, int x, int y)
{


    double dist;
    double ClickedX,ClickedY,ClickedZ;



    //  Print what key the user is hitting
    // printf ("User is hitting the '%c' key.\n", key);
    // printf ("ASCII code is %d.\n", key);

    switch (key)
    {
    //  User hits A key
    case 'a':


        // zooming.  we move our eye along the vector between the eye
        //  and the focus point.
        eyex = (1 - ZOOM_FACTOR) * eyex + focusx * ZOOM_FACTOR;
        eyey = (1 - ZOOM_FACTOR) * eyey + focusy * ZOOM_FACTOR;
        eyez = (1 - ZOOM_FACTOR) * eyez + focusz * ZOOM_FACTOR;
        //glutPostRedisplay();

        break;

    //  User hits Shift + A key
    case 'z':


        eyex = (1 + ZOOM_FACTOR) * eyex - focusx * ZOOM_FACTOR;
        eyey = (1 + ZOOM_FACTOR) * eyey - focusy * ZOOM_FACTOR;
        eyez = (1 + ZOOM_FACTOR) * eyez - focusz * ZOOM_FACTOR;
        //glutPostRedisplay();

        break;
    case 'f':
        //Trova il punto cliccato
        GetClickedPoint(x,y,&ClickedX,&ClickedY,&ClickedZ);



        // metti a fuoco il centro di rotazione
        dist=(CRotX-ClickedX)*(CRotX-ClickedX)+(CRotY-ClickedY)*(CRotY-ClickedY)+(CRotZ-ClickedZ)*(CRotZ-       ClickedZ);
        if (dist<(Diam*Diam*2))//resetta il fuoco solo se non troppo lontano da quelllo corente
        {
            CRotX=ClickedX;
            CRotY=ClickedY;
            CRotZ=ClickedZ;
            focusx=CRotX;
            focusy=CRotY;
            focusz=CRotZ;
        }
        break;
    case 'r': //Reset view to centroid
    {
        CRotX=BarX;
        CRotY=BarY;
        CRotZ=BarZ;
        focusx=BarX;
        focusy=BarY;
        focusz=BarZ;
    }
    break;

    }

    glutPostRedisplay ();
}


//reimposta la window quando viene modificata in dimensioni
void ChangeSize(int w, int h)
{
    GLfloat fAspect;

    // Prevent a divide by zero, when window is too short
    // (you cant make a window of zero width).


    if(h == 0)
        h = 1;
    if (w==0) w=1;

    glViewport(0, 0, w, h);

    fAspect = (GLfloat)w / (GLfloat)h;

    // Reset the coordinate system before modifying
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Set the clipping volume
    gluPerspective(35.0f, fAspect, 1.0f, 50.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}



void DrawClickedPoint()
{

    glPointSize(5);


    glBegin(GL_POINTS);
    glColor3f(1,0,0);
    glVertex3f( CRotX, CRotY, CRotZ);
    glEnd();//end drawing of points

}


//Trova la coordinata del punto cliccato (note le coordiante trovate sono le originali non scalate)
void GetClickedPoint(int x,int y,GLdouble *CRotX, GLdouble *CRotY, GLdouble *CRotZ)
{
//Variabili per calcolare la proiezione
    int viewport[4];
    double modelview[16];
    double projection[16];
    float winX, winY, winZ;


    //trovo la matrice del modello, di proiezione e divisulaizzazion
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);
    winX = (float)x;
    winY = (float)(viewport[3] - y);

    //Trova la "profondità" del pixel in question ed effettua la proiezione
    glReadPixels((int)winX, (int)winY, 1, 1,
                 GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
    gluUnProject(winX, winY, winZ, modelview, projection, viewport,
                 CRotX, CRotY, CRotZ);
}


void sistema_di_riferimento(double X,double Y,double Z)
{
    char* p;


    glEnable(GL_DEPTH_TEST);

    Coord3D origin;


    origin.x=X;
    origin.y=Y;
    origin.z=Z;



    //quadric to draw cylinder
    GLUquadricObj *pObj;    // Temporary, used for quadrics
    pObj = gluNewQuadric();
    gluQuadricDrawStyle(pObj, GLU_FILL);
    gluQuadricNormals(pObj, GLU_SMOOTH);
    gluQuadricOrientation(pObj, GLU_OUTSIDE);
    gluQuadricTexture(pObj, GLU_FALSE);



    glColor3f(0.0f, 0.0f, 0.0f);
    glPushMatrix();
    glTranslatef(origin.x, origin.y, origin.z);
    glutSolidSphere(0.02,25, 25);
    glPopMatrix();



    glColor3f(1.0f, 1.0f, 0.0f);

    //Axes' label
    p="X";
    glRasterPos3d(origin.x+ZOOM_FACTOR*30, 0, origin.z);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *p);

    p="Y";
    glRasterPos3d(origin.x, origin.y+30*ZOOM_FACTOR, 0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *p);

    p="Z";
    glRasterPos3d(origin.x, 0, origin.z+30*ZOOM_FACTOR);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *p);



    //X Axis
    glColor3f(1.0f, 0.00f, 0.00f);
    glPushMatrix();
    glTranslatef(origin.x, 0.0f, 0.0f);
    glRotatef(90.0f,0.0f,1.0f,0.0f);
    gluCylinder(pObj, 0.01, 0.01, ZOOM_FACTOR*20, 10, 1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(origin.x+ZOOM_FACTOR*20, 0.0f, 0.0f);
    glRotatef(90.0f,0.0f,1.0f,0.0f);
    glutSolidCone(2.7*0.025,9*0.025,10,10);
    glPopMatrix();

    //Y Axis
    glColor3f(0.0f, 1.00f, 0.00f);
    glPushMatrix();
    glTranslatef(origin.x, origin.y, 0.0f);
    glRotatef(-90.0f,1.0f,0.0f,0.0f);
    gluCylinder(pObj, 0.01, 0.01, ZOOM_FACTOR*20, 10, 1);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(origin.x, ZOOM_FACTOR*20, 0.0f);
    glRotatef(-90.0f,1.0f,0.0f,0.0f);
    glutSolidCone(2.7*ZOOM_FACTOR,9*ZOOM_FACTOR,10,10);
    glPopMatrix();

    //Z Axis
    glColor3f(0.0f, 0.00f, 1.00f);
    glPushMatrix();
    glTranslatef(origin.x, 0.0f, origin.z);
    gluCylinder(pObj, 0.01, 0.01, ZOOM_FACTOR*20, 10, 1);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(origin.x, 0.0f, ZOOM_FACTOR*20);
    glutSolidCone(2.7*ZOOM_FACTOR,9*ZOOM_FACTOR,10,10);
    glPopMatrix();

    glDisable(GL_DEPTH_TEST);
    gluDeleteQuadric(pObj);
}


void DrawText(string* Text,void*font,int x,int y, float r,float g, float b)
{

    int viewport[4];

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

//dimensioni window
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluOrtho2D(0.0,viewport[2], 0.0, viewport[3]);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glColor3f(r, g, b); // Green
    glRasterPos2i(x, y);

    for (string::iterator i = Text->begin(); i != Text->end(); ++i)
    {
        char c = *i;
        glutBitmapCharacter(font, c);
    }
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}


void DrawTriangles()
{


    int i;
    glClearDepth(1);
    glEnable(GL_DEPTH_TEST);

    //glPushMatrix();



    // Create light components

    /*
    GLfloat  ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat  diffuseLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat  specularLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat  position[] = { 0, 0, 4,1 };




    // Assign created components to GL_LIGHT0
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);*/

//GLfloat global_ambient[]= { 0.5f, 0.5f, 0.5f, 1.0f };
//glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

//GLfloat specular[] = {1.0, 1.0, 1.0, 1.0};
//glLightfv(GL_LIGHT0, GL_SPECULAR, specular);


// glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);//annulla l'effeto della trasformazione di scala
    //glEnable(GL_RESCALE_NORMAL);//from opengl 1.2

    /*
    GLfloat mat_green_diffuse[] = { 0.0, 0.7, 0.1, 1.0 };
    glEnable(GL_COLOR_MATERIAL);
     glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_green_diffuse);
    */

    glColor3f(0.5f, 1.0f, 0.0f);


    for(i=0; i<Nfacets*12; i=i+12)
    {
        glBegin(GL_TRIANGLES);
        glNormal3f(buffer[i],buffer[i+1],buffer[i+2]);
        glVertex3f(buffer[i+3],buffer[i+4],buffer[i+5]);
        glVertex3f(buffer[i+6],buffer[i+7],buffer[i+8]);
        glVertex3f(buffer[i+9],buffer[i+10],buffer[i+11]);
        glEnd();



    }





    //glPopMatrix();
    glDisable(GL_DEPTH_TEST);


    //cout<<"Draw"<<endl;

}

// Called to draw scene
void RenderScene()
{

    GLfloat lightSourcePosition[] = {0, 0, 3, 1};

    // Clear the window with current clearing color
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Save the matrix state and do the rotations

    //glMatrixMode(GL_MODELVIEW); glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();






    glLightfv(GL_LIGHT0, GL_POSITION, lightSourcePosition);

    // move our eye to the most recent place
//	gluLookAt(0, 0, -5, 0,0,0, 0,1,0);
//gluLookAt(eyex, eyey, eyez, 0, 0, 0, 0,1,0);
    glScalef(FactScale,FactScale,FactScale);
    gluLookAt(eyex, eyey, eyez, focusx, focusy, focusz, 0,1,0);

    // Rotate the camera
    // rotate our box in the x and y directions

//Rotation around rotation center

    glTranslatef(CRotX, CRotY, CRotZ);
    glRotatef(g_xRotation, 0,1,0);
    glRotatef(g_yRotation, 1,0,0);
    glTranslatef(-CRotX, -CRotY, -CRotZ);


    // Draw something

    DrawTriangles();


    DrawClickedPoint();





    DrawHelpLines();







    // Restore the matrix state
    //glPopMatrix();

    // Buffer swap
    glutSwapBuffers();

}


void GetModelData(double* p,int N)
{
    //For graphical purpose

    double Mx,My,Mz,mx,my,mz;



    Mx=-HUGE_VAL;
    My=-HUGE_VAL;
    Mz=-HUGE_VAL;
    mx=HUGE_VAL;
    my=HUGE_VAL;
    mz=HUGE_VAL;


    double dist;



    BarX=Mean(p,N,3);
    BarY=Mean(&p[1],N,3);
    BarZ=Mean(&p[2],N,3);
    MinMax(p,N,&Mx,&mx,3);
    MinMax(&p[1],N,&My,&my,3);
    MinMax(&p[2],N,&Mz,&mz,3);

    //Calcolo la diagonale del modello

    Diam=sqrt((Mx-mx)*(Mx-mx)+(My-my)*(My-my)+(Mz-mz)*(Mz-mz));
}


void DrawHelpLines()
{

    /* Tipi di testo supportati da glut
    #  GLUT_BITMAP_8_BY_13
    # GLUT_BITMAP_9_BY_15
    # GLUT_BITMAP_TIMES_ROMAN_10
    # GLUT_BITMAP_TIMES_ROMAN_24
    # GLUT_BITMAP_HELVETICA_10
    # GLUT_BITMAP_HELVETICA_12
    # GLUT_BITMAP_HELVETICA_18
    */

    void * font = GLUT_BITMAP_HELVETICA_12;



    int y=70;//inizio tabella in y
    int h=15;//spazio fra un riga e l'altra
    int x=10;//sapzio con il margine sinistro

    float r=0;//colorazione testo
    float g=1;
    float b=0;
    string Text;

    Text = "ZOOM: mouse R";
    DrawText(&Text,font,x,y, r,g,b);
    Text = "PAN: mouse R&L";
    y=y-h;
    DrawText(&Text,font,x,y, r,g,b);
    Text = "ROTATE: mouse L";
    y=y-h;
    DrawText(&Text,font,x,y, r,g,b);
    Text = "FOCUS: f";
    y=y-h;
    DrawText(&Text,font,x,y, r,g,b);
    Text = "RESET: r";
    y=y-h;
    DrawText(&Text,font,x,y, r,g,b);

}

}