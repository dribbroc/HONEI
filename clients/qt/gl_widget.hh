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


#ifndef HONEI_GUARD_GL_WIDGET_HH
#define HONEI_GUARD_GL_WIDGET_HH

#include <QGLWidget>
#include <QTimer>
#include <simulation_controller.hh>
class GLWidget : public QGLWidget
{
    Q_OBJECT

    public:
        GLWidget(QWidget * father = 0);
        ~GLWidget();

        QSize minimumSizeHint() const;
        QSize sizeHint() const;

        unsigned long get_sim_id();
        void solver_start_stop();
        void simulation_reload();
        void simulation_load(unsigned long id);
        void hud_on_off();

    protected:
        virtual void initializeGL();
        virtual void paintGL();
        virtual void resizeGL( int width, int height );
        virtual void mousePressEvent( QMouseEvent *event );
        virtual void mouseMoveEvent( QMouseEvent *event );
        virtual void wheelEvent( QWheelEvent *event );

    private:
        GLuint m_object;
        int m_xRot, m_yRot, m_zRot;
        float m_xTrans, m_yTrans, m_zTrans;
        QPoint m_lastPos;
        float m_backgroundcolor[3];
        float m_curr_rot;

        DenseMatrix<float> * _idle_hb;
        DenseMatrix<float> * _idle_b;

        QTimer * _solver_timer;

        bool _solver_precision_flag;
        bool _solver_start_stop_flag;
        bool _render_idle_flag;
        bool _hud_on_flag;
        unsigned long _sim_w, _sim_h, _sim_id;
        float _sim_dx, _sim_dy, _sim_dt, _sim_tau;

        solver_type _current_solver;

        SimulationController<float>* _sim_control_float;
        SimulationController<double>* _sim_control_double;

        void _set_x_rotation(int value);
        void _set_z_rotation(int value);
        void _normalize_angle(int * angle) const;

        template <typename Prec_>
        void _render_matrix( DenseMatrix<Prec_> & matrix, float r, float g, float b, float a);

        void _render_hud();

    private slots:
        void solver_event();
};

#endif
