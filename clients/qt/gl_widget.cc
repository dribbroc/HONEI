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

const unsigned int RotTimer = 0;	// 0 = workproc, i.e., when there are no more UI events in the queue


    GLWidget::GLWidget(QWidget * father )
: QGLWidget( QGLFormat(QGL::SampleBuffers), father),			// enables multi-sampling
      m_object(0), m_xRot(0), m_yRot(0), m_zRot(0),
      m_xTrans(0), m_yTrans(0.), m_zTrans(1.),
      m_curr_rot(0), _solver_precision_flag(true), _solver_start_stop_flag(false), _render_idle_flag(false), _hud_on_flag(true),
      _sim_w(50), _sim_h(50), _sim_id(0)
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

#ifdef HONEI_SSE
    _current_solver = sse_full_dry;
#else
    _current_solver = cpu_full_dry;
#endif

    _idle_hb = new DenseMatrix<float>(50, 50, float(0));
    _idle_b = new DenseMatrix<float>(50, 50, float(0));
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
    glEnable( GL_DEPTH_TEST );

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
    glTranslatef(m_xTrans, m_yTrans, 1.);
    glRotatef(m_xRot,1.0, 0.0, 0.0);
    glRotatef(m_yRot,0.0, 1.0, 0.0);
    glRotatef(m_zRot,0.0, 0.0, 1.0);

    //mz_Trans is now used for scaling the scene aka zoom-in/zoom-out
    glScalef(1.f * m_zTrans, 1.f * m_zTrans, 100.0f * m_zTrans);
    glEnable (GL_BLEND);

    if(_solver_precision_flag)
    {
        if(!_render_idle_flag)
            glTranslatef(float(-(_sim_control_float->get_b()).columns())/2., float(-(_sim_control_float->get_b()).rows())/2., 0.);
        else
            glTranslatef(float(-(_idle_b->columns()))/2., float(-(_idle_b->rows()))/2., 0.);
    }
    else
    {
        if(!_render_idle_flag)
            glTranslatef(float(-(_sim_control_double->get_b()).columns())/2., float(-(_sim_control_double->get_b()).rows())/2., 0.);
        else
            glTranslatef(float(-(_idle_b->columns()))/2., float(-(_idle_b->rows()))/2., 0.);
    }

    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    if(!_render_idle_flag)
    {
        if (_solver_precision_flag)
        {
            DenseMatrix<float> b(_sim_control_float->get_b());
            _render_matrix(b, 0.0, 0.8, 0.0, 1.);
            DenseMatrix<float> hb(_sim_control_float->get_hb());
            _render_matrix(hb, 0.0, 0.0, 0.8, 0.5);
        }
        else
        {
            DenseMatrix<double> b(_sim_control_double->get_b());
            _render_matrix(b, 0.0, 0.8, 0.0, 1.);
            DenseMatrix<double> hb(_sim_control_double->get_hb());
            _render_matrix(hb, 0.0, 0.0, 0.8, 0.5);
        }
    }
    else
    {
            _render_matrix(*_idle_hb, 0.0, 0.8, 0.0, 1.);
            _render_matrix(*_idle_b, 0.0, 0.0, 0.8, 0.5);
    }

    _render_hud();
}

void GLWidget::resizeGL( int w, int h )		// = width & height
{
    int side = qMin(w, h);
    glViewport((w - side) / 2, (h - side) / 2, side, side);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-100, +100, +100, -100, -500., 100.0);
    //glFrustum(-100., +100.0, -100.0, +100.0, -1000.0, 1000.0);
}

    template <typename Prec_>
void GLWidget::_render_matrix(DenseMatrix<Prec_> & matrix, float r, float g, float b, float a)
{
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

void GLWidget::simulation_reload()
{
    if(_solver_start_stop_flag)
        _solver_start_stop_flag = false;

    if(!_render_idle_flag)
        _render_idle_flag = true;

    if(_solver_precision_flag)
        _sim_control_float->reload_simulation();
    else
        _sim_control_double->reload_simulation();

    _render_idle_flag = false;

    updateGL();
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
        m_xTrans -= 0.5 * dy;
        m_xTrans += 0.5 * dx;
    }
    else if ( e->buttons() & Qt::LeftButton )
    {
        _set_x_rotation(m_xRot + 4 * dy);
        _set_z_rotation(m_zRot + 4 * dx);
    }
    else if ( e->buttons() & Qt::MidButton )
    {
        m_yTrans += 0.5 * dy;
        m_yTrans -= 0.5 * dx;
    }

    m_lastPos = e->pos();
    e->accept();
    updateGL();
}

void GLWidget::wheelEvent( QWheelEvent * e )
{
    m_zTrans += float((e->delta() / 120) * 0.05);

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

void GLWidget::_set_z_rotation( int angle )
{
    _normalize_angle( &angle );
    if (angle != m_yRot)
    {
        m_zRot = angle;
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

void GLWidget::_render_hud()
{
    if(_hud_on_flag)
    {
        std::string hud_text_1("backend:     ");
        if(_solver_precision_flag)
            hud_text_1 += _sim_control_float->get_backend_info();
        else
            hud_text_1 += _sim_control_double->get_backend_info();

        std::string hud_text_2("precision:    ");
        if(_solver_precision_flag)
            hud_text_2 += "single";
        else
            hud_text_2 += "double";

        std::string hud_text_3("simulation:  ");
        if(_solver_precision_flag)
            hud_text_3 += _sim_control_float->get_simulation_info();
        else
            hud_text_2 += _sim_control_float->get_simulation_info();

        const char * txt_1 = hud_text_1.c_str();
        const char * txt_2 = hud_text_2.c_str();
        const char * txt_3 = hud_text_3.c_str();

        glMatrixMode( GL_PROJECTION );
        glPushMatrix();
        glLoadIdentity();
        glMatrixMode( GL_MODELVIEW );

        glDisable( GL_DEPTH_TEST );
        glColor3f( 1.0, 1.0, 1.0 );

        glPushMatrix();
        glLoadIdentity();

        renderText( 10, 20, QString(txt_1), QFont("System", 8) );
        renderText( 10, 35, QString(txt_2), QFont("System", 8) );
        renderText( 10, 50, QString(txt_3), QFont("System", 8) );

        glPopMatrix();

        glMatrixMode( GL_PROJECTION );
        glPopMatrix();
        glMatrixMode( GL_MODELVIEW );
        glEnable(GL_DEPTH_TEST);
    }
}

void GLWidget::hud_on_off()
{
    _hud_on_flag = !_hud_on_flag;

    updateGL();
}

unsigned long GLWidget::get_sim_id()
{
    return _sim_id;
}

void GLWidget::simulation_load(unsigned long id)
{
    _sim_id = id;
}
