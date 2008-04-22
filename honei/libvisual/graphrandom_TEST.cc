/* vim: set sw=4 sts=4 et foldmethod=syntax nu : */

#include <honei/libvisual/enginegraph.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>
#include <honei/libgraph/graph.hh>
#include <honei/libgraph/position.hh>
#include <cstdlib>

//
using namespace honei;
using namespace tests;
using namespace std;
using namespace gl_globals;

template <typename Tag_, typename DataType_, typename GraphTag_>
class EngineGraphTest :
    public BaseTest
{ 
    private:
        typedef EngineGraph<Tag_, DataType_, GraphTag_> Engine;
        int _nodeCount;
        double _p;
        static inline DataType_ randomWeight()
        {
            return DataType_((std::rand() % 4) + 1);
        }
    public:
        EngineGraphTest(const std::string & type, int nodeCount, double p) :
            BaseTest("EngineGraph test<" + type + ">"),
            _nodeCount(nodeCount),
            _p(p)
        {
            register_tag(Tag_::name);
        }

        
        virtual void run() const
        {
            int i =1;
            int * pi = &i;
            Graph<DataType_>  g(_nodeCount, 2);
            for (int i = 0; i < _nodeCount; i++)
                g.add_node(i, randomWeight());
                
            //for (int j = 1; j < _nodeCount; j++)
            //        g.addEdge(0, j, randomWeight());

        for (int i = 0; i < _nodeCount; ++i)
            for (int j = i+1; j < _nodeCount; ++j)
                        if (std::rand() < RAND_MAX * _p)
                            g.add_edge(i, j, randomWeight());

            std::cout << "\nCalculate Position\n";
            std::cout << *g.edges();
            
            Engine::setTestCase(g, new Positions<Tag_, DataType_, GraphTag_>(g, (DataType_)1), int(_nodeCount * 0.1));
            
            
            char * c = "Test: Engine";
            char ** cp = &c;
            glutInit(pi,cp);
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
            glutInitWindowSize(screen_width, screen_height);
            glutInitWindowPosition(0,0);
            glutCreateWindow("Libgraph Render engine 2007 Thorsten Deinert based on: Render- engine 1.0 (c) 2007 Markus Geveler");
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
//EngineGraphTest<tags::CPU::SSE, float, methods::WeightedKamadaKawai> engine_test_double("wkk float", 11);
//EngineGraphTest<tags::CPU::SSE, double, methods::WeightedFruchtermanReingold> engine_test_double2("wkk double big", 200, .01345);
EngineGraphTest<tags::CPU::SSE, float, methods::WeightedKamadaKawai> engine_test_double2("wkk float big", 200, 0.33);
