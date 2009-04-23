/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of HONEI. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <QtGui>
#include <QtOpenGL>

#include <cmath>

#include <clients/qt/gl_widget.hh>
#include <iostream>
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
      m_object(0), m_xRot(-45.), m_yRot(0), m_zRot(45),
      m_xTrans(-4.), m_yTrans(-10.0), m_zTrans(-45.0),
      m_curr_rot(0), m_numFlakeRec(2), _solver_precision_flag(true), _solver_start_stop_flag(false)
      //m_timer()
{
    if ( ! format().sampleBuffers() )
        fprintf(stderr, "Could not get sample buffer; no polygon anti-aliasing!");

    m_backgroundcolor[0] = 0.2f;  m_backgroundcolor[1] = 1.0f;  m_backgroundcolor[2] = 0.4f;

    if(_solver_precision_flag)
    {
        _sim_control_float = new SimulationController<float>();
        _sim_control_double = 0;
    }
    else
    {
        _sim_control_float = 0;
        _sim_control_double = new SimulationController<double>();
    }
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
    return QSize(640, 480);
}

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

    //Insert HONEI gl initialization here
    glEnable( GL_DEPTH_TEST ); //Later needed for lighting

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

    //start the solver timer
    _solver_timer = new QTimer(this);
    connect(_solver_timer, SIGNAL(timeout()), SLOT(solver_event()));
    _solver_timer->start();
}

void GLWidget::paintGL()
{
    glClearColor(0., 0., 0., 0.);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //render HONEI solver output
    glTranslatef(m_xTrans, m_yTrans, m_zTrans);
    glRotatef(m_xRot,1.0, 0.0, 0.0);
    glRotatef(m_yRot,0.0, 1.0, 0.0);
    glRotatef(m_zRot,0.0, 0.0, 1.0);

    glScalef(1.0f, 1.0f, 100.0f);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if (_solver_precision_flag)
    {
        _render_matrix(_sim_control_float->get_b(), 0.0, 0.8, 0.0, 1.);
        _render_matrix(_sim_control_float->get_hb(), 0.0, 0.0, 0.8, 0.5);
    }
    else
    {
        _render_matrix(_sim_control_double->get_b(), 0.0, 0.8, 0.0, 1.);
        _render_matrix(_sim_control_double->get_hb(), 0.0, 0.0, 0.8, 0.5);
    }
}

void GLWidget::resizeGL( int w, int h )		// = width & height
{
    int side = qMin(w, h);
    glViewport((w - side) / 2, (h - side) / 2, side, side);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(-1.5, +1.5, +1.5, -1.5, 1.0, 150.0);
    glFrustum(-1.1, +1.0, -1.0, +1.0, 1.0, 210.0);
}

template <typename Prec_>
void GLWidget::_render_matrix(DenseMatrix<Prec_> & matrix, float r, float g, float b, float a)
{
    //glPushMatrix();
    //glTranslatef(-matrix.columns()/2, -matrix.rows()/2, 0.);

    glBegin(GL_QUADS);
    for(unsigned int i = 0 ; i < matrix.rows() - 1 ; ++i)
    {
        for(unsigned int j = 0 ; j < matrix.columns() - 1 ; ++j)
        {
            glColor4f(r, g, b, a);
            glVertex3d(i,j, matrix[i][j]);
            glColor4f(r, g + 0.2, b + 0.2, a);
            glVertex3d(i,j+1, matrix[i][j+1]);
            glVertex3d(i+1,j+1, matrix[i+1][j+1]);
            glVertex3d(i+1,j, matrix[i+1][j]);
        }
    }
    //glPopMatrix();
    glEnd();
}

void GLWidget::solver_event()
{
    if(_solver_start_stop_flag)
    {
        if(_solver_precision_flag)
            _sim_control_float->do_timestep();
        else
            _sim_control_double->do_timestep();

        updateGL();
    }
}

void GLWidget::solver_start_stop()
{
    _solver_start_stop_flag = !_solver_start_stop_flag;
}


void GLWidget::mousePressEvent( QMouseEvent * e )
{
    m_lastPos = e->pos();
}

void GLWidget::mouseMoveEvent( QMouseEvent * e )
{
    int dx = e->x() - m_lastPos.x();
    int dy = e->y() - m_lastPos.y();

    bool ctrl_key = e->modifiers() & Qt::MetaModifier;

    if ( (e->buttons() & Qt::RightButton) ||
            ctrl_key )
    {
        m_zTrans += 0.5 * dy;
        m_xTrans += 0.5 * dx;
    }
    else if ( e->buttons() & Qt::LeftButton )
    {
        _set_x_rotation(m_xRot + 8 * dy);
        _set_y_rotation(m_yRot + 8 * dx);
    }
    m_lastPos = e->pos();
    e->accept();
    updateGL();
}

void GLWidget::_set_x_rotation( int angle )
{
    _normalize_angle( &angle );
    if (angle != m_xRot)
    {
        m_xRot = angle;
        updateGL();
    }
}

void GLWidget::_set_y_rotation( int angle )
{
    _normalize_angle( &angle );
    if (angle != m_yRot)
    {
        m_yRot = angle;
        updateGL();
    }
}
void GLWidget::_normalize_angle( int *angle ) const
{
    while (*angle < 0)
        *angle += 360;
    while (*angle > 360 )
        *angle -= 360;
}
