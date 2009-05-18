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

#ifndef HONEI_GUARD_WINDOW_HH
#define HONEI_GUARD_WINDOW_HH

#include <QMainWindow>
#include <QScrollArea>
#include <QWidget>
#include <QMenu>
#include <QAction>

#include <clients/qt/gl_widget.hh>

class QSlider;
class GLWidget;

class Window : public QMainWindow
{
    Q_OBJECT

    public:
        Window();

    private:
        QWidget * centralWidget;
        QScrollArea * glArea;
        GLWidget *glWidget;

        ///Simulation menu related members:
        QMenu * _simulation_menu;
        QAction * _solver_start_stop_action;
        QAction * _simulation_reload_action;
        QMenu * _simulation_load_submenu;
        std::vector<QAction *> _simulation_load_actions;
        QAction * _exit_action;

        ///Solver menu related members:
        QMenu * _solver_menu;
        QMenu * _solver_backend_submenu;
        QAction * _solver_backend_cpu_action;
#ifdef HONEI_SSE
        QAction * _solver_backend_sse_action;
#endif
#ifdef HONEI_CUDA
        QAction * _solver_backend_cuda_action;
#endif

        ///HUD menu related members:
        QMenu * _hud_menu;
        QAction * _hud_on_off_action;

        void create_menu();

    private slots:
        void _solver_start_stop();
        void _simulation_reload();
        void _hud_on_off();
        void _simulation_load();
        void _solver_backend_cpu();
        void _solver_backend_sse();
        void _solver_backend_cuda();
};

#endif
