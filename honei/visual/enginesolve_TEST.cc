/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LibVisual C++ library. LibVisual is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibVisual is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/visual/enginesolve.hh>
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>


using namespace honei;
using namespace tests;
using namespace std;
using namespace gl_globals;

template <typename Tag_, typename DataType_>
class EngineTest :
    public BaseTest
{
    public:
        EngineTest(const std::string & type) :
            BaseTest("Engine test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            int i =1;
            int * pi = &i;

            char c[] = "SWE";
            char * cc = c;
            char ** cp = &cc;
            glutInit(pi,cp);
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
            //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);
            glutInitWindowSize(screen_width, screen_height);
            glutInitWindowPosition(0,0);
            glutCreateWindow("LibSWE Render- engine 1.0 (c) 2007 Markus Geveler");
            glutDisplayFunc(Engine::display);
            glutIdleFunc(Engine::display);
            glutReshapeFunc(Engine::resize);
            glutKeyboardFunc(Engine::keyboard);
            glutSpecialFunc(Engine::keyboard_s);
            glutMouseFunc(Engine::mouse);
            menu_id_rendering = glutCreateMenu(Engine::menu_rendering);
            glutAddMenuEntry("Toggle fill mode", 2);
            glutAddMenuEntry("Toggle ground", 3);
            glutAddMenuEntry("Toggle water", 4);
            glutAddMenuEntry("Toggle shading", 5);
            glutAddMenuEntry("Toggle primitive type", 6);
            glutAddMenuEntry("Toggle alpha blending", 7);
            menu_id_scenario = glutCreateMenu(Engine::menu_scenario);
            glutAddMenuEntry("todo", 1);
            menu_id_main = glutCreateMenu(Engine::menu_main);
            glutAddMenuEntry("Restart scenario", 0);
            glutAddSubMenu("Rendering", menu_id_rendering);
            glutAddSubMenu("Scenarios", menu_id_scenario);
            glutAddMenuEntry("Exit programm", 10);
            glutAttachMenu(GLUT_RIGHT_BUTTON);
            glutMotionFunc(NULL);
            glutPassiveMotionFunc(NULL);
            Engine::init();
            glutMainLoop();
            TEST_CHECK(true);
        }
};
EngineTest<tags::CPU, double> engine_test_double("double");
