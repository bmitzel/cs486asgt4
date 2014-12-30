/* Programmer: Brian Mitzel
 * Email: bmitzel@csu.fullerton.edu
 * Course: CPSC 486
 * Assignment: 4 - Separating Axes
 * Due Date: 12/17/2014
 *
 * Filename: moving_boxes.cpp
 *
 * This program demonstrates the separating axes algorithm
 * for collision detection between 2D polygons. A number of
 * squares are rendered on the screen with variable starting
 * locations and velocities. When a square collides with
 * another square, its velocity vector is reversed. Also,
 * when a square collides with an edge of the window, it is
 * reflected.
 *
 * As the squares move around the screen, they will slowly
 * fade out to zero intensity (black). However, when one
 * square collides with another square, it will change color
 * and appear at full intensity again.
 *
 * This source code was adapted from a demo by Professor
 * Michael Shafae.
 */

#ifdef _WIN32
#include <Windows.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#define _USE_MATH_DEFINES
#include <cmath>
#define sprintf sprintf_s
#define strdup _strdup
#else
#include <cstdio>
#include <cstdlib>
#include <ctime>
#endif

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glut.h>
#else
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#ifdef FREEGLUT
#include <GL/freeglut_ext.h>
#endif

#include <list>
#include <sys/time.h>

#ifndef NDEBUG
#define msglError( ) _msglError( stderr, __FILE__, __LINE__ )
#else
#define msglError( )
#endif

typedef struct rectangle{
  float center[3];
  float half_width;
  float half_height;
  int   hue;
  float intensity;
  float translation[3];
  float velocity[3];
}rectangle_t;

/* Identify the four edges of the screen */
enum screen_edges {LEFT, RIGHT, BOTTOM, TOP};

/* Identify the different types of collisions */
enum collision_types {NONE, HORIZONTAL_SCREEN_EDGE, VERTICAL_SCREEN_EDGE, RECTANGLE};

#define NUMRECTS 22
rectangle_t* eRects[4];        /* screen edge rectangles */
rectangle_t* gRects[NUMRECTS]; /* graphical rectangles */
int gTimer;
bool gPrintRectangleFlag;

/* The colors of the graphical squares */
static const float hues[2][3] = {{0.40f, 1.00f, 0.30f} /* light green */
                                ,{0.20f, 0.00f, 1.00f} /* azure */
                                };

/**
 * Returns the number of seconds that have elapsed since the previous call to the function
 * @return - The number of seconds that have elapsed since the previous call to the function
 */
double getDeltaTime()
{
    static bool firstTime = true;
    static timeval previous;
    struct timeval now;
    struct timeval elapsed;

    if (firstTime)
    {
        firstTime = false;
        gettimeofday(&previous, NULL);
    }

    gettimeofday(&now, NULL);
    elapsed.tv_sec = now.tv_sec - previous.tv_sec;
    elapsed.tv_usec = now.tv_usec - previous.tv_usec;

    if (elapsed.tv_usec < 0)
    {
        --elapsed.tv_sec;
        elapsed.tv_usec += 1000000;
    }

    previous = now;
    double n = (static_cast<uint64_t>(elapsed.tv_sec) * static_cast<uint64_t>(1000000) +
            static_cast<uint64_t>(elapsed.tv_usec)) / 1000000.0;

    return n;
} /* getDeltaTime() */

void fcopy(float* src, float* dst, size_t n){
  for(size_t i = 0; i < n; i++){
    dst[i] = src[i];
  }
}

int matPrint4d(FILE *f, double *m){
  int i;
  int accum = 0;
  for( i = 0; i < 4; i++ ){
    accum += fprintf( f, "%.16f %.16f %.16f %.16f\n", m[i + 0], m[i + 4], m[i + 8], m[i + 12] );
  }
  return(accum);
}

void vecPerp2f(float* p, float* v)
{
    float temp[] = {v[0], v[1]};

    p[0] =  temp[1];
    p[1] = -temp[0];
}

void vecDifference3f(float *c, float *a, float *b){
  int i;
  for(i = 0; i < 3; i++){
    c[i] = a[i] - b[i];
  }
}

float vecDot3f(float *a, float *b){
  float accumulate = 0.0;
  int i;
  for(i = 0; i < 3; i++){
    accumulate += a[i] * b[i];
  }
  return( accumulate );
}

void vecScalarMult3f(float *b, const float *a, float s){
  int i;
  for(i = 0; i < 3; i++){
    b[i] = a[i] * s;
  }
}

bool _msglError( FILE *out, const char *filename, int line ){
  bool ret = false;
  GLenum err = glGetError( );
  while( err != GL_NO_ERROR ) {
    ret = true;
    fprintf( out, "GL ERROR:%s:%d: %s\n",
      filename, line, (char*)gluErrorString( err ) );
    err = glGetError( );
  }
  return( ret );
}

void printTopOfBothStacks( char *msg ){
  double m[16];
  msglError( );
  if( msg != NULL ){
    fprintf( stderr, "%s\n", msg );
  }
  fprintf( stderr, "Top of Projection Stack:\n" );
  glGetDoublev( GL_PROJECTION_MATRIX, m );
  puts("proj");
  matPrint4d( stdout, m );
  fprintf( stderr, "\n\n" );
  glGetDoublev( GL_MODELVIEW_MATRIX, m );
  puts("mv");
  matPrint4d( stdout, m );
  msglError( );
}

float frand( ){
  float r = float(rand( )) / float(RAND_MAX);
  if( rand( ) % 2 ){
    r *= -1;
  }
  return r;
}

void printRectangle(int i, rectangle_t* r){
  printf("rectangle %d{\n", i);
  printf("\tcenter = {%g, %g, %g}\n", r->center[0], r->center[1], r->center[2]);
  printf("\thalf width = %g\n", r->half_width);
  printf("\thalf height = %g\n", r->half_height);
  printf("\thue = %d\n", r->hue);
  printf("\tintensity = %g\n", r->intensity);
  printf("\ttranslation = {%g, %g, %g}\n", r->translation[0], r->translation[1], r->translation[2]);
  printf("\tvelocity = {%g, %g, %g}\n", r->velocity[0], r->velocity[1], r->velocity[2]);
  printf("}\n");
}

void rectangleVertices(float v[12], rectangle_t* r){
  v[0] = r->center[0] + r->half_width + r->translation[0];
  v[1] = r->center[1] + r->half_height + r->translation[1];
  v[2] = r->center[2] + r->translation[2];
  
  v[3] = r->center[0] - r->half_width + r->translation[0];
  v[4] = r->center[1] + r->half_height + r->translation[1];
  v[5] = r->center[2] + r->translation[2];
  
  v[6] = r->center[0] - r->half_width + r->translation[0];
  v[7] = r->center[1] - r->half_height + r->translation[1];
  v[8] = r->center[2] + r->translation[2];
  
  v[9] = r->center[0] + r->half_width + r->translation[0];
  v[10] = r->center[1] - r->half_height + r->translation[1];
  v[11] = r->center[2] + r->translation[2];
}

/**
 * Returns an array of vertices for a given rectangle
 * @param rectangle - The rectangle used for calculating the vertices
 * @param vertices[] - The array of vertices to return
 */
void getVertices(rectangle_t* rectangle, float vertices[][3])
{
    vertices[0][0] = rectangle->center[0] + rectangle->translation[0] - rectangle->half_width;
    vertices[0][1] = rectangle->center[1] + rectangle->translation[1] - rectangle->half_height;
    vertices[0][2] = 0.0f;
    vertices[1][0] = rectangle->center[0] + rectangle->translation[0] + rectangle->half_width;
    vertices[1][1] = rectangle->center[1] + rectangle->translation[1] - rectangle->half_height;
    vertices[1][2] = 0.0f;
    vertices[2][0] = rectangle->center[0] + rectangle->translation[0] + rectangle->half_width;
    vertices[2][1] = rectangle->center[1] + rectangle->translation[1] + rectangle->half_height;
    vertices[2][2] = 0.0f;
    vertices[3][0] = rectangle->center[0] + rectangle->translation[0] - rectangle->half_width;
    vertices[3][1] = rectangle->center[1] + rectangle->translation[1] + rectangle->half_height;
    vertices[3][2] = 0.0f;
} /* getVertices() */

/**
 * Calculates if a vertex from one rectangle overlaps any of the vertices from a second rectangle
 * @param vertices - The set of vertices for the second rectangle
 * @param direction - The axis used to test for separation
 * @param point - The vertex from the first rectangle
 * @return - 1 if all of the vertices are separated along the given axis in the positive direction;
 * -1 if all of the vertices are separated along the given axis in the negative direction;
 * otherwise, 0
 */
int whichSide(float vertices[][3], float direction[], float point[])
{
    int positive = 0;
    int negative = 0;

    /* Check for separating axes with each vertex of the second rectangle */
    for (int i = 0; i < 4; i++)
    {
        float vector[3]; /* the vector given by vertex - point */
        float t; /* used to determine the sign of the angle between the axis and the vector */

        vecDifference3f(vector, vertices[i], point);
        t = vecDot3f(direction, vector);

        if (t > 0)
        {
            ++positive;
        }
        else if (t < 0)
        {
            ++negative;
        }

        if (positive && negative)
        {
            return 0;
        }
    }

    return positive ? +1 : -1;
} /* whichSide() */

/**
 * Checks for 2D intersection between two rectangles using separating axes
 * @param rOne - The first rectangle
 * @param rTwo - The second rectangle
 * @return - True if the rectangles overlap;
 * otherwise, false
 */
bool intersection2D(rectangle_t* rOne, rectangle_t* rTwo)
{
    /* Check if the rectangles are the same */
    if (rOne == rTwo)
    {
        return false;
    }
    /* Check if the rectangles overlap */
    else
    {
        float vOne[4][3]; /* rectangle rOne's vertices */
        float vTwo[4][3]; /* rectangle rTwo's vertices */

        /* Get the vertices of each rectangle */
        getVertices(rOne, vOne);
        getVertices(rTwo, vTwo);

        /* Separate axes for each vertex in the first rectangle
         * versus all the vertices in the second rectangle */
        for (int iOne = 0, iTwo = 3; iOne < 4; iTwo = iOne, iOne++)
        {
            float direction[3];

            vecDifference3f(direction, vOne[iOne], vOne[iTwo]);

            if (whichSide(vTwo, direction, vOne[iOne]) > 0)
            {
                return false;
            }
        }

        /* Separate axes for each vertex in the second rectangle
         * versus all the vertices in the first rectangle */
        for (int iOne = 0, iTwo = 3; iOne < 4; iTwo = iOne, iOne++)
        {
            float d[3];

            vecDifference3f(d, vOne[iOne], vOne[iTwo]);
            vecPerp2f(d, d);

            if (whichSide(vTwo, d, vOne[iOne]) > 0)
            {
                return false;
            }
        }
    }

    return true;
} /* intersection2D() */

/**
 * Initializes the four squares representing the edges of the screen
 * @param translations - An array containing the translations for each square
 */
void initEdgeSquares(float translations[][3])
{
    for (int i = 0; i < 4; i++)
    {
        fcopy(translations[i], eRects[i]->translation, 3);
        eRects[i]->half_width = 1.0f;
        eRects[i]->half_height = 1.0f;

        eRects[i]->hue = 0;
        eRects[i]->intensity = 0.0f;

        for (int j = 0; j < 3; j++)
        {
            eRects[i]->center[j] = 0.0f;
            eRects[i]->velocity[j] = 0.0f;
        }
    }
} /* initEdgeSquares() */

#define DIM 0.1 /* the dimension of each of the graphical squares */
/**
 * Initializes a graphical square
 * @param r - The square (rectangle) to initialize
 * @param color - The RGB color values for the square
 */
void initUnitSquare(rectangle_t* r, int colorHue, float colorIntensity)
{
    /* Initialize the half-dimensions */
    r->half_width = DIM * 0.5;
    r->half_height = DIM * 0.5;

    /* Initialize the hue and intensity */
    r->hue = colorHue;
    r->intensity = colorIntensity;

    /* Initialize the center, velocity, and translation;
     * ensure no rectangles are overlapping when spawned
     */
    bool overlapping;

    do
    {
        overlapping = false;

        for (int i = 0; i < 2; i++)
        {
            /* Initialize the center */
            r->center[i] = 0.0;

            /* Initialize the velocity */
            float temp = frand();

            if (temp < 0.0f)
            {
                r->velocity[i] = (temp * 0.01f - 0.01f) / 3.333f; /* range [-0.006, -0.003] */
            }
            else /* temp >= 0.0f */
            {
                r->velocity[i] = (temp * 0.01f + 0.01f) / 3.333f; /* range [0.003, 0.006] */
            }

            /* Initialize the translation */
            r->translation[i] = frand();

            /* Avoid overlapping another rectangle */
            for (int j = 0; r != gRects[j]; j++)
            {
                if (intersection2D(r, gRects[j]))
                {
                    overlapping = true;
                }
            }

            /* Avoid overlapping a screen edge */
            for (int j = 0; j < 4; j++)
            {
                if (intersection2D(r, eRects[j]))
                {
                    overlapping = true;
                }
            }
        }
    } while (overlapping);

    r->center[2] = 0.0;
    r->velocity[2] = 0.0;
    r->translation[2] = 0.0;
} /* initUnitSquare() */

/**
 * Resets a graphical square
 * @param r - The square (rectangle) to reset
 */
void resetUnitSquare(rectangle_t* r)
{
    int hue = 0;

    for (int i = 0; i < NUMRECTS; i++)
    {
        /* The intensity of the color as a float in the range [0.3, 1.0], incrementing by 0.05 */
        float intensity = static_cast<float>(rand() % 15 + 6) / 20.0f;

        initUnitSquare(gRects[i], hue, intensity);

        /* Alternate hues when spawning the squares */
        hue = (hue + 1) % 2;
    }
//  for(int i = 0; i < 2; i++){
//    r->center[i] = 0.0;
//    r->velocity[i] = frand( ) * 0.01;
//    r->translation[i] = frand( );
//  }
//  r->center[2] = 0.0;
//  r->velocity[2] = 0.0;
//  r->translation[2] = 0.0;
} /* resetUnitSquare() */

/**
 * Updates a graphical rectangle
 * @param r - The rectangle to update
 * @param deltaTime - The number of seconds since the last update
 */
void updateRectangle(rectangle_t* r, double deltaTime)
{
    /* Update the rectangle's intensity */
    r->intensity = r->intensity - (0.000166667f / deltaTime); /* divide to scale with gTime */
    fprintf(stderr, "deltaTime = %f\n", deltaTime);

    /* Update the rectangle's translation */
    for (int i = 0; i < 3; i++)
    {
        r->translation[i] += r->velocity[i];

        if (gPrintRectangleFlag)
        {
            printRectangle(i, r);
        }
    }
} /* updateRectangle() */

void drawRectangle(rectangle_t* r){
  float *c = r->center;
  float hw = r->half_width;
  float hh = r->half_height;
  float color[3];
  vecScalarMult3f(color, ::hues[r->hue], r->intensity);
  glColor3fv(color);
  glPushMatrix( );
  glTranslatef(r->translation[0], r->translation[1], r->translation[2]);

  glBegin(GL_TRIANGLES);
  glVertex3f(c[0] + hw, c[1] + hh, 0.0);
  glVertex3f(c[0] - hw, c[1] - hh, 0.0);
  glVertex3f(c[0] + hw, c[1] - hh, 0.0);
  
  glVertex3f(c[0] + hw, c[1] + hh, 0.0);
  glVertex3f(c[0] - hw, c[1] + hh, 0.0);
  glVertex3f(c[0] - hw, c[1] - hh, 0.0);
  glEnd( );

  glPopMatrix( );
}

/**
 * Deallocates memory before exiting
 */
void cleanUp()
{
    /* Deallocate the screen edge squares */
    for (int i = 0; i < 4; i++)
    {
      delete eRects[i];
    }

    /* Deallocate the graphical squares */
    for (int i = 0; i < NUMRECTS; i++)
    {
      delete gRects[i];
    }
} /* cleanUp() */

void keyboardUp( unsigned char key, int x, int y ){
  switch( key ){
    //default:
    //fprintf( stderr, "The key '%c' (%d) just went up.\n", key, key );
  }
}

void keyboard( unsigned char key, int x, int y ){
  //fprintf( stderr, "You pushed '%c' (%d).\n", key, key );
  switch( key ){
  case 27: // The 'esc' key
  case 'q':
    cleanUp();
#ifdef FREEGLUT
    glutLeaveMainLoop( );
#else
    exit( 0 );
#endif
    break;
  case '1':
    gTimer = 500;
    break;
  case '2':
    gTimer = 16;
    break;
  case 'd':
    gPrintRectangleFlag = !gPrintRectangleFlag;
    break;
  case 'r':
    for(int i = 0; i < NUMRECTS; i++){
      resetUnitSquare(gRects[i]);
    }
    break;
  default:
    fprintf( stderr, "Unregistered: You pushed '%c' (%d).\n", key, key );
    break;
  }
}

void display( ){
  msglError( );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  for(int i = 0; i < NUMRECTS; i++){
    drawRectangle(gRects[i]);
  }
  glutSwapBuffers( );
  msglError( );
}

void reshape( int width, int height ){
  fprintf( stderr, "reshape\n" );
  if (height == 0){
    height = 1;
  }
  glViewport( 0, 0, (GLsizei)width, (GLsizei)height );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );
}

void initOpenGL( ){
  glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
  glEnable( GL_DEPTH_TEST );
  glDepthRange(0.0, 1.0);
  glDepthFunc(GL_LEQUAL);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );
  msglError( );
}

/**
 * Initializes all of the program data, including the graphical and screen edge squares
 */
void initData()
{
    /* Initialize the program data */
    srand(time(NULL));
    gTimer = 16;
    gPrintRectangleFlag = false;

    /* Translations of the squares representing the four edges of the screen */
    float eTranslations[4][3] = {{-2.0f,  0.0f, 0.0f}
                                ,{ 2.0f,  0.0f, 0.0f}
                                ,{ 0.0f, -2.0f, 0.0f}
                                ,{ 0.0f,  2.0f, 0.0f}
                                };

    /* Allocate memory for the screen edge squares */
    for (int i = 0; i < 4; i++)
    {
        eRects[i] = new rectangle_t;
    }

    /* Initialize the screen edge squares */
    initEdgeSquares(eTranslations);

    /* Initialize the graphical squares */
    int hue = 0;

    for (int i = 0; i < NUMRECTS; i++)
    {
        /* The intensity of the color as a float in the range [0.3, 1.0], incrementing by 0.05 */
        float intensity = static_cast<float>(rand() % 15 + 6) / 20.0f;

        gRects[i] = new rectangle_t;
        initUnitSquare(gRects[i], hue, intensity);

        /* Alternate hues when spawning the squares */
        hue = (hue + 1) % 2;
    }
} /* initData() */

/**
 * Resolves collisions between pairs of overlapping rectangles
 * @param resolvedCollisions - An array storing previously resolved collisions
 * @return - True if there was a new collision detected;
 * otherwise, false
 */
bool resolveRectangleCollisions(char resolvedCollisions[])
{
    bool isColliding = false;
    char collisions[NUMRECTS] = {0}; /* counts the number of collisions for each rectangle */

    /* Initialize the array of colliding rectangles */
    for (int i = 0; i < NUMRECTS - 1; i++)
    {
        for (int j = i + 1; j < NUMRECTS; j++)
        {
            if (intersection2D(gRects[i], gRects[j]))
            {
                /* Set the colliding flag to true */
                isColliding = true;

                /* Add the collision to the collisions array */
                ++collisions[i];
                ++collisions[j];
            }
        }
    }

    /* Resolve the collisions */
    for (int i = 0; i < NUMRECTS; i++)
    {
        if (collisions[i])
        {
            /* Check if this rectangle already has a resolved collision */
            if (RECTANGLE != resolvedCollisions[i])
            {
                /* Do nothing if the rectangle's collision has already been resolved;
                 * otherwise, resolve the new collision
                 */
                if (NONE == resolvedCollisions[i])
                {
                    /* Undo the rectangle's last translation to avoid overlapping rectangles */
                    for (int j = 0; j < 3; j++)
                    {
                        gRects[i]->translation[j] -= gRects[i]->velocity[j];
                    }

                    /* Reverse the rectangle's velocity */
                    vecScalarMult3f(gRects[i]->velocity, gRects[i]->velocity, -1.0f);
                }
                else if (HORIZONTAL_SCREEN_EDGE == resolvedCollisions[i])
                {
                    /* Reverse the y component of the rectangle's velocity */
                    gRects[i]->velocity[1] *= -1.0f;
                }
                else if (VERTICAL_SCREEN_EDGE == resolvedCollisions[i])
                {
                    /* Reverse the x component of the rectangle's velocity */
                    gRects[i]->velocity[0] *= -1.0f;
                }

                /* Adjust the rectangle's hue and intensity */
                gRects[i]->hue = (gRects[i]->hue + 1) % 2; /* alternate the hue */
                gRects[i]->intensity = 1.0f;               /* maximize the intensity */

                /* Update the rectangle in the resolved collisions array */
                resolvedCollisions[i] = RECTANGLE;
            }
        }
    }

    return isColliding;
} /* resolveRectangleCollisions() */

/**
 * Resolves collisions of rectangles overlapping screen edges
 * @param resolvedCollisions - An array storing previously resolved collisions
 * @return - True if there was a new collision detected;
 * otherwise, false
 */
bool resolveEdgeCollisions(char resolvedCollisions[])
{
    bool isColliding = false;

    /* Check if a rectangle is colliding with a screen edge */
    for (int i = 0; i < NUMRECTS; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (intersection2D(gRects[i], eRects[j]))
            {
                /* Set the colliding flag to true */
                isColliding = true;

                /* Undo the last translation to avoid overlapping the screen edge */
                for(int k = 0; k < 3; k++)
                {
                    gRects[i]->translation[k] -= gRects[i]->velocity[k];
                }

                /* Reflect the rectangle off the screen edge */
                if (LEFT == j || RIGHT == j)
                {
                    gRects[i]->velocity[0] *= -1.0f;

                    /* Update the rectangle in the resolved collisions array */
                    resolvedCollisions[i] = HORIZONTAL_SCREEN_EDGE;
                }
                else /* BOTTOM or TOP */
                {
                    gRects[i]->velocity[1] *= -1.0f;

                    /* Update the rectangle in the resolved collisions array */
                    resolvedCollisions[i] = VERTICAL_SCREEN_EDGE;
                }
            }
        }
    }

    return isColliding;
} /* resolveEdgeCollisions() */

/**
 * Updates the rectangles, resolves any collisions, and redraws the display
 * @param x - Unused value parameter
 */
void update(int x)
{
    double deltaTime = getDeltaTime();

    /* Update the rectangles */
    for(int i = 0; i < NUMRECTS; i++)
    {
        updateRectangle(gRects[i], deltaTime);
    }

    /* Initialize an empty array of resolved collisions */
    char resolvedCollisions[NUMRECTS] = {NONE};

    /* Test for collisions with other rectangles or screen edges */
    bool rectanglesCollided = resolveRectangleCollisions(resolvedCollisions);
    bool screenEdgesCollided = resolveEdgeCollisions(resolvedCollisions);

    /* Continue resolving collisions until no new collisions occur */
    while (rectanglesCollided || screenEdgesCollided)
    {
        /* The test for collisions with screen edges only needs to be called once */
        screenEdgesCollided = false;

        /* Repeat the test for new collisions with other rectangles */
        rectanglesCollided = resolveRectangleCollisions(resolvedCollisions);
    }

    glutTimerFunc(gTimer, update, 0);
    glutPostRedisplay();

    msglError();
} /* update() */

int main(int argc, const char* argv[]){
  glutInit( &argc, (char**)argv );
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowSize( 800, 800);
  glutCreateWindow( "Moving Boxes" );
#ifdef FREEGLUT
  //glutFullScreen( );
#endif
  initOpenGL( );
  initData( );
  glutKeyboardFunc( keyboard );
  glutKeyboardUpFunc( keyboardUp );
  glutDisplayFunc( display );
  glutReshapeFunc( reshape );
  glutTimerFunc( gTimer, update, 0 );
  glutMainLoop( );
  return(0);
}
