/* vim: set sw=4 sts=4 et nofoldenable : */
#include <QtGui>
#include <QtOpenGL>
#include <QMouseEvent>

#include <math.h>

#include "glwidget.hh"


// switch on the following if you want the scene to be drawn in wireframe
#undef WIREFRAME
#define ANTIALIAS

static const float SphereColor[3] = { 0.2f, 0.1f, 0.9f };
const float PIf = 3.14159265358979323846f;
const unsigned int NumSphereSlices = 30;
const float RotSpeed = 0.1f;
const unsigned int RotTimer = 0;	// 0 = workproc, i.e., when there are no more UI events in the queue


    GLWidget::GLWidget(QWidget * father )
: QGLWidget( QGLFormat(QGL::SampleBuffers), father),			// enables multi-sampling
      m_object(0), m_xRot(0), m_yRot(0), m_zRot(0),
      m_xTrans(0.0), m_yTrans(0.0), m_zTrans(-5.0),
      m_curr_rot(0.0), m_numFlakeRec(2)
      //m_timer()
{
    if ( ! format().sampleBuffers() )
        fprintf(stderr, "Could not get sample buffer; no polygon anti-aliasing!");

    m_backgroundcolor[0] = 0.2f;  m_backgroundcolor[1] = 1.0f;  m_backgroundcolor[2] = 0.4f;
}


GLWidget::~GLWidget()
{
    makeCurrent();
    glDeleteLists(m_object, 1);
}


QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(400, 400);
}


/*void GLWidget::setXRotation( int angle )
{
    normalizeAngle( &angle );
    if (angle != m_xRot)
    {
        m_xRot = angle;
        mit xRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::setYRotation( int angle )
{
    normalizeAngle( &angle );
    if (angle != m_yRot)
    {
        m_yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::normalizeAngle( int *angle ) const
{
    while (*angle < 0)
        *angle += 360;
    while (*angle > 360 )
        *angle -= 360;
}


void GLWidget::setFlakeRec(int rec)
{
    printf("rec=%d \n", rec );
    m_numFlakeRec = rec;
    updateGL();
}


void GLWidget::setBackgroundColor( int v )
{
    m_backgroundcolor[2] = static_cast<float>( v ) / 100.0f;

    updateGL();
}
*/

void GLWidget::initializeGL()
{
    // query OpenGL
    const GLubyte * strVersion = glGetString( GL_VERSION );
    printf("\nGL version = %s\n", strVersion );
    const GLubyte * strExt = glGetString( GL_EXTENSIONS );
    fputs("GL extensions: ", stdout);
    puts( reinterpret_cast<const char *>(strExt) );
    GLint multisamplebufs;
    glGetIntegerv( GL_SAMPLE_BUFFERS, &multisamplebufs );
    GLint multisamples;
    glGetIntegerv( GL_SAMPLES, &multisamples );
    printf("multisampling = %d\nnum samples = %d\n",
            static_cast<int>(multisamplebufs), static_cast<int>(multisamples) );

    // set opengl states
    glEnable( GL_DEPTH_TEST );	// is already default in Qt

    glCullFace( GL_BACK );		// back-face culling
    glEnable(GL_CULL_FACE);		// how much performance does this yield?
    glFrontFace( GL_CCW );

    // specify lighting
    GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat mat_shininess[] = { 50.0 };
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
    glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
    GLfloat ambient_light[] = { 0.1f, 0.1f, 0.5f, 1.0f};
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, ambient_light );

    GLfloat light0_pos[] = { -1.0f, 1.0f, -1.0f, 0.0 };
    GLfloat light0_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };		// yes, there are defaults for light 0 --
    GLfloat light0_diffuse[] = { 0.5f, 0.5f, 0.5f, 1.0f };		// but better make it explicit
    GLfloat light0_specular[] = { 0.9f, 0.9f, 0.9f, 1.0f };
    glLightfv( GL_LIGHT0, GL_AMBIENT, light0_ambient );
    glLightfv( GL_LIGHT0, GL_DIFFUSE, light0_diffuse );
    glLightfv( GL_LIGHT0, GL_SPECULAR, light0_specular );
    glLightfv( GL_LIGHT0, GL_POSITION, light0_pos );
    glEnable( GL_LIGHT0 );

    GLfloat light1_pos[] = { 1.0f, -1.0f, 1.0f, 0.0f };			// remember: we look along -z
    GLfloat light1_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat light1_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat light1_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv( GL_LIGHT1, GL_AMBIENT, light1_ambient );
    glLightfv( GL_LIGHT1, GL_DIFFUSE, light1_diffuse );
    glLightfv( GL_LIGHT1, GL_SPECULAR, light1_specular );
    glLightfv( GL_LIGHT1, GL_POSITION, light1_pos );
    glEnable( GL_LIGHT1 );

    glEnable( GL_LIGHTING );

    // we ask GL to renormalize normal vectors,
    // because glScale() does not preserve normal vectors ;-(
    //glEnable( GL_NORMALIZE );
    glEnable( GL_RESCALE_NORMAL );		// no difference on my Mac ...

#ifdef ANTIALIAS
#ifdef WIREFRAME
    // anti-aliasing of lines
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
#else /* no WIREFRAME */
    // anti-aliasing of polygons and lines
    glEnable( GL_MULTISAMPLE );
#endif
#endif

    // create object
    /*m_object = glGenLists(1);
    glNewList(m_object, GL_COMPILE);
    makeSphere( NumSphereSlices, NumSphereSlices, 1.0 );
    glEndList();
*/
    // now we can start the timer for the animation function timerEvent()
    startTimer(RotTimer);
}


// Warning: don't make paintGL() and resizeGL() const!
// otherwise, the rendering will produce complete garbage.
// (is this a bug in Qt? g++? opengl?)

void GLWidget::paintGL()
{
    //m_timer.stop();		// make sure that glFinish() or glSwapBuffers() is counted, too;
    //m_timer.start();	// thus, we stop the timer *after* Qt has given control back to us.

    float r, g, b;
    HsvToRgb( m_backgroundcolor[0], m_backgroundcolor[1], m_backgroundcolor[2], &r, &g, &b );
    glClearColor( r, g, b, 0.0 );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef( m_xTrans, m_yTrans, m_zTrans );
    glRotatef( m_xRot, 1.0, 0.0, 0.0 );				// angle is in degrees
    glRotatef( m_yRot, 0.0, 1.0, 0.0 );
    glRotatef( m_zRot, 0.0, 0.0, 1.0 );

    // animate spheres
    //renderSphereFlake( 0.5, m_numFlakeRec );

    // display framerate, not counting current frame (which is not being displayed yet)
    //display_info( 1.0 / (m_timer.average() / 1E9) );
}


// Display some text as 2D overlay on the 3D scene
/*void GLWidget::display_info( float framerate )
{
    const char * hud_text = "frame rate = %3d";
    char txt[100];
    sprintf( txt, hud_text, static_cast<int>( framerate ) );

    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode( GL_MODELVIEW );

    glDisable( GL_DEPTH_TEST );
    glColor3f( 1.0, 0.0, 0.0 );

    glPushMatrix();
    glLoadIdentity();

    renderText( 10, 20, QString(txt), QFont("System", 12) );		// Qt function

    glPopMatrix();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );

    glEnable(GL_DEPTH_TEST);
}*/



void GLWidget::resizeGL( int w, int h )		// = width & height
{
    int side = qMin(w, h);
    glViewport((w - side) / 2, (h - side) / 2, side, side);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(-1.5, +1.5, +1.5, -1.5, 1.0, 150.0);
    glFrustum(-1.0, +1.0, -1.0, +1.0, 1.0, 150.0);
}


/*void GLWidget::timerEvent( QTimerEvent * )
{
    m_curr_rot += RotSpeed;
    if ( m_curr_rot > 360.0 )
        m_curr_rot = 0.0;

    updateGL();
}*/


/*void GLWidget::mousePressEvent( QMouseEvent * e )
{
    m_lastPos = e->pos();
}


void GLWidget::mouseMoveEvent( QMouseEvent * e )
{
    int dx = e->x() - m_lastPos.x();
    int dy = e->y() - m_lastPos.y();

    bool ctrl_key = e->modifiers() & Qt::MetaModifier;		// only needed for Mac OS X, but doesn't hurt on other OSes

    if ( (e->buttons() & Qt::RightButton) ||
            ctrl_key )
    {
        m_zTrans += 0.5 * dy;
        m_xTrans += 0.5 * dx;
    }
    else if ( e->buttons() & Qt::LeftButton )
    {
        setXRotation(m_xRot + 8 * dy);
        setYRotation(m_yRot + 8 * dx);
    }
    m_lastPos = e->pos();
    e->accept();
    updateGL();
}
*/



/*
void GLWidget::renderSphereFlake( const float scaling, unsigned int n_recursions ) const
{
    glCallList(m_object);

    if ( n_recursions <= 1 )
        return;

    glPushMatrix();
    glRotatef( m_curr_rot, 1.0, 1.0, 1.0 );
    glTranslatef( 1.5, 0.0, 0.0 );
    glScalef( scaling, scaling, scaling );
    renderSphereFlake( scaling, n_recursions-1 );
    glPopMatrix();

    glPushMatrix();
    glRotatef( m_curr_rot, 1.0, 1.0, 1.0 );
    glTranslatef( -1.5, 0.0, 0.0 );
    glScalef( scaling, scaling, scaling );
    renderSphereFlake( scaling, n_recursions-1 );
    glPopMatrix();

    glPushMatrix();
    glRotatef( m_curr_rot, 1.0, 1.0, 1.0 );
    glTranslatef( 0.0, -1.5, 0.0 );
    glScalef( scaling, scaling, scaling );
    renderSphereFlake( scaling, n_recursions-1 );
    glPopMatrix();

    glPushMatrix();
    glRotatef( m_curr_rot, 1.0, 1.0, 1.0 );
    glTranslatef( 0.0, 1.5, 0.0 );
    glScalef( scaling, scaling, scaling );
    renderSphereFlake( scaling, n_recursions-1 );
    glPopMatrix();

    glPushMatrix();
    glRotatef( m_curr_rot, 1.0, 1.0, 1.0 );
    glTranslatef( 0.0, 0.0, 1.5 );
    glScalef( scaling, scaling, scaling );
    renderSphereFlake( scaling, n_recursions-1 );
    glPopMatrix();

    glPushMatrix();
    glRotatef( m_curr_rot, 1.0, 1.0, 1.0 );
    glTranslatef( 0.0, 0.0, -1.5 );
    glScalef( scaling, scaling, scaling );
    renderSphereFlake( scaling, n_recursions-1 );
    glPopMatrix();
}*/



/*typedef float Point[3];						// should be done using proper class
typedef double Pointd[3];					// and double precision
*/

/**  Compute a normal
 *
 * @param n			the normal (out)
 * @param p1,..,p3	the three points
 *
 * Computes a normal from the three points, i.e.,
 * n = (p2-p1) x (p3-p1).
 *
 * @implementation
 *   Of course, in real life, this fct would be in a utility library ...
 **/

/*void calcNormal( Point n, const Point p1, const Point p2, const Point p3 )
{
    Pointd a, b, nt;
    for ( unsigned int i = 0; i < 3; i ++ )
        a[i] = p2[i] - p1[i];
    for ( unsigned int i = 0; i < 3; i ++ )
        b[i] = p3[i] - p1[i];

    // cross product
    nt[0] = a[1] * b[2] - a[2] * b[1];
    nt[1] = a[2] * b[0] - a[0] * b[2];
    nt[2] = a[0] * b[1] - a[1] * b[0];
    // normalize
    double l = nt[0]*nt[0] + nt[1]*nt[1] + nt[2]*nt[2];
    if ( l > 1E-12 )
    {
        l = sqrtf( l );
        n[0] = nt[0] / l;
        n[1] = nt[1] / l;
        n[2] = nt[2] / l;
    }
    else
    {
        n[0] = nt[0];
        n[1] = nt[1];
        n[2] = nt[2];
    }
}

*/

/**  Send a normal to OpenGL
 *
 * @param vertex		the vertex array
 * @param p1,..,p3	the indices of the first three points into the vertex array
 *
 * Computes a normal from the three points specified by p1-p3 and vertex.
 * It is assumed that these three points are the first three vertices of the
 * primitive to be sent after this to OpenGL.
 * It is also assumed that the three vertices are in counter-clockwise order.
 *
 * @warning
 *   No bounds check on the indices p1 .. p3 is done!
 **/

/*void sendNormal( const Point * vertex,
        unsigned int p1, unsigned int p2, unsigned int p3 )
{
    Point n;
    calcNormal( n, vertex[p1], vertex[p2], vertex[p3] );
    glNormal3f( n[0], n[1], n[2] );
}
*/


/**  Render a quad
 *
 * @param vertex		the vertex array
 * @param p1,..,p4	the indices of the four points into the vertex array
 *
 * Emits 4 glVertex3f commands, where the coords come from the vertex array.
 *
 * @warning
 *   No bounds check on the indices p1 .. p4 is done!
 **/
/*
void renderQuad( const Point * vertex, unsigned int p1, unsigned int p2,
        unsigned int p3, unsigned int p4 )
{
    sendNormal( vertex, p1, p2, p3 );
    glVertex3fv( vertex[p1] );
    glVertex3fv( vertex[p2] );
    glVertex3fv( vertex[p3] );
    glVertex3fv( vertex[p4] );
}
*/
/**  Render a triangle
 *
 * @param vertex		the vertex array
 * @param p1,..,p3	the indices of the three points into the vertex array
 *
 * Emits 3 glVertex3f commands, where the coords come from the vertex array.
 *
 * @warning
 *   No bounds check on the indices p1 .. p3 is done!
 **/
/*
void renderTri( const Point * vertex,
        unsigned int p1, unsigned int p2, unsigned int p3 )
{
    sendNormal( vertex, p1, p2, p3 );
    glVertex3fv( vertex[p1] );
    glVertex3fv( vertex[p2] );
    glVertex3fv( vertex[p3] );
}
*/


/**  Render a sphere
 *
 * @param lati,longi	number of latitudes and longitudes
 * @param radius		radius of the sphere
 *
 * @implementation
 *   copied from Y ... (and found a bug -- after 15 years ...)
 *
 **/
/*
void GLWidget::makeSphere( unsigned int lati, unsigned int longi,
        const float radius ) const
{
    if ( lati < 2 )
        lati = 2;
    if ( longi < 3 )
        longi = 3;
    unsigned int np = 2 + longi*(lati - 1);		// num vertices

    // create points between equator and poles
    Point * vertex = new Point[np] ;
    float cd = cosf( 2*PIf / longi );
    float sd = sinf( 2*PIf / longi );
    for ( unsigned int j = 0; j < lati-1; j ++ )
    {
        float h = cosf( static_cast<float>(j+1) / lati * PIf );
        float ci = sinf( static_cast<float>(j+1) / lati * PIf );
        float si = 0.0;

        for ( unsigned int i = 0; i < longi; i ++ )
        {
            vertex[ j*longi+i ][0] = ci * radius;
            vertex[ j*longi+i ][1] = h  * radius;
            vertex[ j*longi+i ][2] = si * radius;
            float ci1 = ci*cd - si*sd;
            float si1 = ci*sd + si*cd;
            ci = ci1;
            si = si1;
        }
    }

    // create pole points
    vertex[np-2][0] = 0.0;			// north
    vertex[np-2][1] = radius;
    vertex[np-2][2] = 0.0;
    vertex[np-1][0] = 0.0;			// south
    vertex[np-1][1] = -radius;
    vertex[np-1][2] = 0.0;

    // for ( unsigned int i = 0; i < np; i ++ )
    // printf("v %d: % 6f % 6f % 6f\n", i, vertex[i][0], vertex[i][1], vertex[i][2] );

    glColor3fv( SphereColor );

    // create quadrangles
#ifdef WIREFRAME
    glBegin( GL_LINE_LOOP );
#else
    glBegin( GL_QUADS );
#endif
    for ( unsigned int j = 1; j < lati-1; j ++ )
    {
        for ( unsigned int i = 0; i < longi-1; i ++ )
        {
            unsigned int k = (j-1) * longi + i;
            renderQuad( vertex, k, k+1, k+longi+1, k+longi );
        }
        unsigned int k = (j-1) * longi;
        renderQuad( vertex, k+longi-1, k, k+longi, k+2*longi-1 );
    }
    glEnd();

    // create pole triangles
#ifdef WIREFRAME
    glBegin( GL_LINE_LOOP );
#else
    glBegin( GL_TRIANGLES );
#endif
    for ( unsigned int i = 0; i < longi-1; i ++ )
        renderTri( vertex, np-2, i+1, i );				// north cap
    renderTri( vertex, np-2, 0, longi-1 );

    for ( unsigned int i = 0; i < longi-1; i ++ )
        renderTri( vertex, np-2-longi + i, np-2-longi + i+1, np-1 );	// south cap
    renderTri( vertex, np-3, np-2-longi, np-1 );
    glEnd();

    delete [] vertex;
}
*/

float clamp( float low, float val, float high)
{
    if ( val < low )
        return low;
    else if ( val > high )
        return high;
    else
        return val;
}



/**  Convert from HSV color space to RGB color space
 *
 * @param h,s,v		hue, saturation, value (in)
 * @param r,g,b		red, green, blue       (out)
 *
 * Reminder: saturation=1.0 is border of color cone,
 *     saturation=0.0 is grey, value=0.0 is apex of cone which is black
 * The values of h,s,v will be clamped to [0,1].
 *
 * @implementation
 *   copied from Y ...
 *
 **/

void GLWidget::HsvToRgb( float h, float s, float v,
        float *r, float *g, float *b ) const
{
    clamp( 0.0f, v, 1.0f );
    clamp( 0.0f, s, 1.0f );

    if ( s < 1E-5 )
        *r = *g = *b = v;
    else
    {
        clamp( 0.0f, h, 1.0f );

        h = h * 360 / 60.0f;
        int i = static_cast<int>( floor(h) );
        float f = h - i;
        float p = v * (1 - s);
        float q = v * (1 - (s*f) );
        float t = v * (1 - (s*(1-f)) );
        switch( i )
        {
            case 0: *r = v; *g = t; *b = p; break;
            case 1: *r = q; *g = v; *b = p; break;
            case 2: *r = p; *g = v; *b = t; break;
            case 3: *r = p; *g = q; *b = v; break;
            case 4: *r = t; *g = p; *b = v; break;
            case 5: *r = v; *g = p; *b = q; break;
            default: fprintf( stderr, "GLWidget::HsvToRgb: bug!"
                             " i=%d shouldn't occur!\n" , i );
                     fprintf( stderr, "( h=%f s=%f v=%f )\n", h, s, v );
        }
    }
}

