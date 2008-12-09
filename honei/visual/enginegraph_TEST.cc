/* vim: set sw=4 sts=4 et foldmethod=syntax nu : */

#include <honei/visual/enginegraph.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>
#include <honei/graph/graph.hh>
#include <honei/graph/position.hh>
#include <honei/graph/test_scenario.hh>


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
        Graph<DataType_> * _graph;
    public:
        EngineGraphTest(const std::string & type, Graph<DataType_> * graph) :
            BaseTest("EngineGraph test<" + type + ">"),
            _graph(graph)
        {
            register_tag(Tag_::name);
        }


        virtual void run() const
         {
            int i =1;
            int * pi = &i;
   /*         Graph<DataType_> g(_nodeCount, 2);
            for (int i = 0; i < _nodeCount; i++)
                g.add_node(i, 1 + std::rand() % 4);

            for (int j = 1; j < _nodeCount; j++)
                    g.add_edge(0, j, 1);
   //         std::cout << "coordinates g: " << *g->coordinates();
   //         std::cout << "edge matrix\n" << *g->edges();

            std::cout << "\nCalculate Position";
     */



            Positions<Tag_, DataType_, GraphTag_> * fine(new Positions<Tag_, DataType_, GraphTag_>(*_graph, (DataType_)1));
        //    fine->update(0, 256);
            Engine::setTestCase(*_graph, fine, 1, -140);


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
#ifdef HONEI_SSE
//EngineGraphTest<tags::CPU::SSE, float, methods::WeightedKamadaKawai> engine_test_double("wkk float", 11);
EngineGraphTest<tags::CPU::SSE, float, methods::WeightedFruchtermanReingold> engine_test_float("WFR Test Grid", TestScenario<float>::Grid(26,26));
//EngineGraphTest<tags::CPU::SSE, float, methods::WeightedKamadaKawai> engine_test_double2("wkk float big", 200);
#endif
