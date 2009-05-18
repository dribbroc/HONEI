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

#include <clients/qt/gl_widget.hh>
#include <clients/qt/window.hh>
#include <honei/lbm/scenario_collection.hh>

Window::Window()
{
    centralWidget = new QWidget;
    setCentralWidget(centralWidget);

    glArea = new QScrollArea;
    glWidget = new GLWidget;

    glArea->setWidget(glWidget);
    glArea->setWidgetResizable(true);
    glArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    glArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    glArea->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    glArea->setMinimumSize(50, 50);

    QGridLayout *mainLayout = new QGridLayout;

    mainLayout->addWidget(glArea);
    centralWidget->setLayout(mainLayout);
    setWindowTitle( "HONEI QT OpenGL client" );

    create_menu();

    resize(640,480);
}

void Window::create_menu()
{

    ///Simulation menu:
    _simulation_menu = menuBar() -> addMenu(tr("&Simulation"));
    _solver_start_stop_action = new QAction(tr("&Start/Stop"), this);
    _solver_start_stop_action->setShortcut(tr("Ctrl+S"));
    connect(_solver_start_stop_action, SIGNAL(triggered()),
            this, SLOT(_solver_start_stop()));

    _simulation_menu->addAction(_solver_start_stop_action);


    _simulation_reload_action = new QAction(tr("&Reload"), this);
    _simulation_reload_action->setShortcut(tr("Ctrl+R"));
    connect(_simulation_reload_action, SIGNAL(triggered()),
            this, SLOT(_simulation_reload()));

    _simulation_menu->addAction(_simulation_reload_action);

    _simulation_load_submenu = _simulation_menu -> addMenu(tr("&Load predefined"));
    for (unsigned long i(0) ; i <= ScenarioCollection::get_stable_scenario_count() ; ++i)
    {
        QAction * action = new QAction(tr(ScenarioCollection::get_scenario_descr(i).c_str()), this);
        action->setCheckable(true);
        _simulation_load_actions.push_back(action);
        _simulation_load_submenu->addAction(action);
        connect(action, SIGNAL(triggered()),
                this, SLOT(_simulation_load()));
        if(glWidget->get_sim_id() == i)
            action->setChecked(true);
    }

    ///HUD menu:
    _hud_menu = menuBar() -> addMenu(tr("&HUD"));
    _hud_on_off_action = new QAction(tr("&On/Off"), this);
    _hud_on_off_action->setShortcut(tr("Ctrl+I"));
    connect(_hud_on_off_action, SIGNAL(triggered()),
            this, SLOT(_hud_on_off()));

    _hud_menu->addAction(_hud_on_off_action);
}

void Window::_solver_start_stop()
{
    glWidget->solver_start_stop();
}

void Window::_simulation_reload()
{
    glWidget->simulation_reload();
}

void Window::_hud_on_off()
{
    glWidget->hud_on_off();
}

void Window::_simulation_load()
{
    unsigned long current_sim_id(glWidget->get_sim_id());

    _simulation_load_actions[current_sim_id]->setChecked(false);

    for(unsigned long i(0) ; i <= ScenarioCollection::get_stable_scenario_count() ; ++i)
    {
        if (_simulation_load_actions[i]->isChecked() == true)
        {
            glWidget->simulation_load(i);
            break;
        }
    }
}
