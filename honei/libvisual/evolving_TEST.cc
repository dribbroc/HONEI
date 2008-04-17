/* vim: set sw=4 sts=4 et foldmethod=syntax nu : */

#include <honei/libvisual/enginegraph.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
#include <string>
#include <honei/libgraph/graph.hh>
#include <honei/libgraph/evolving_graph.hh>
#include <honei/libgraph/position.hh>
#include <honei/libgraph/test_scenario.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace gl_globals;

template <typename Tag_, typename DataType_, typename GraphTag_>
class EngineEvolvingGraphTest :
    public BaseTest
{
    private:
        typedef EngineGraph<Tag_, DataType_, GraphTag_> Engine;
        typedef Node<DataType_> ND;
        int _nodes, _slices;
    public:
        EngineEvolvingGraphTest(const std::string & type, int nodes, int slices) :
            BaseTest("EngineGraph test<" + type + ">"),
            _nodes(nodes),
            _slices(slices)
            
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            int i =1;
            int * pi = &i;
            int nps = _nodes * _slices;
            int sizes[] =  {2, 3, 4, 5, 6};
            //EvolvingGraph<DataType_> * eg = TestScenario<DataType_>::Evolving(5, sizes);
            
            EvolvingGraph<DataType_> eg(2, DataType_(5));
            for (int t(0); t < _slices; ++t)
            {
                Graph<DataType_> & g(eg.add_timeslice((t+1) * _nodes - t*0));
                for (int n(t*0); n < (t+1)*_nodes; ++n)
                    g.add_node(t*_nodes  + n);
                for (int n(t*0); n < (t+1)*_nodes; ++n)
                   // for (int m(n+1); m < (t+1)*_nodes; ++m)
                        g.add_edge(t*_nodes  + n, t*_nodes  +  (n+1) % ((t+1)*_nodes), 1);
            }
            /*

            Graph<DataType_>  t1(eg.add_timeslice(4));
            t1.add_node(1,1);
            t1.add_node(2,1);
            t1.add_node(5,2);
            t1.add_node(6,2);
            
            t1.add_edge(1, 2, 1);
            t1.add_edge(2, 5, 1);
            t1.add_edge(5, 6, 1);
            t1.add_edge(6, 1, 1);
            //eg.addTimeslice(t1);

            Graph<DataType_>  t2(eg.add_timeslice(4));
            t2.add_node(1,1);
            t2.add_node(2,1);
            t2.add_node(3,1);
            t2.add_node(6,3);
            
            t2.add_edge(1,2,1);
            t2.add_edge(1,3,1);
            t2.add_edge(1,6,2);
            //eg.addTimeslice(t2);   

            Graph<DataType_>  t3(eg.add_timeslice(3));   
            t3.add_node(1,1);
            t3.add_node(3,1);
            t3.add_node(4,1);
            t3.add_edge(1,4,1);
            t3.add_edge(3,4,1);
            t3.add_edge(1,3,1);
            //eg.addTimeslice(t3);     
            */       
            
            std::cout << "Scenario created!\n";
            //eg.reassemble_graph();
            std::cout << "Assembled\n";
            std::cout << "Evolving: " << eg.coordinates()->rows() << " Nodes, " << eg.edges()->rows() << "² Edges\n";
         //   std::cout << "coordinates eg: " << *eg.coordinates();
         //   std::cout << "edge matrix\n" << *eg.edges();
            
            
            std::cout << "\nCalculate Position\n";
            //Positions<Tag_, DataType_, methods::WeightedKamadaKawai> positions(&gl_globals::graph, (DataType_)1);
            
            Engine::setTestCase(eg, new Positions<Tag_, DataType_, GraphTag_>(eg, (DataType_)1), 1,-40);
            

            char * c = "Test: Engine";
            char ** cp = &c;
            glutInit(pi,cp);
            glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
            //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);
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
EngineEvolvingGraphTest<tags::CPU::SSE, float, methods::WeightedFruchtermanReingold> engine_test_double("wkk float", 7, 6);
//EngineEvolvingGraphTest<tags::CPU::SSE, float, methods::WeightedFruchtermanReingold> engine_test_double("wkk double");
