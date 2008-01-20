/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/libvisual/engine.hh>
#include <honei/libswe/solver.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
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

            char * c = "Test: Engine";
            char ** cp = &c;
            glutInit(pi,cp);
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
            //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);
            glutInitWindowSize(screen_width, screen_height);
            glutInitWindowPosition(0,0);
            glutCreateWindow("LibSWE Render- engine 1.0 (c) 2007 Markus Geveler");
            glutDisplayFunc(Engine::display);
            glutIdleFunc(Engine::display);
            glutReshapeFunc(Engine::resize);
            glutKeyboardFunc (Engine::keyboard);
            glutSpecialFunc (Engine::keyboard_s);
            Engine::init();
            glutMainLoop();
            TEST_CHECK(true);
        }
};
EngineTest<tags::CPU, double> engine_test_double("double");
