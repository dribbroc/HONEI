/* vim: set sw=4 sts=4 et foldmethod=syntax nu : */

#include <honei/libvisual/enginegraph.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
#include <string>
#include <honei/libgraph/graph.hh>
#include <honei/libgraph/evolving_graph.hh>
#include <honei/libgraph/position.hh>

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
    public:
        EngineEvolvingGraphTest(const std::string & type) :
            BaseTest("EngineGraph test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            int i =1;
            int * pi = &i;
            
            EvolvingGraph<DataType_> eg(2, 1000);

            eg.addNode(new Node<DataType_>(1, 1));
            eg.addNode(new Node<DataType_>(2, 1));
            eg.addNode(new Node<DataType_>(3, 1));
            eg.addNode(new Node<DataType_>(4, 1));
            eg.addNode(new Node<DataType_>(5, 2));
            eg.addNode(new Node<DataType_>(6, 2));

            Graph<DataType_> * t1 = new Graph<DataType_>(4, 2);
            t1->addNode(eg.getNode(1));
            t1->addNode(eg.getNode(2));
            t1->addNode(eg.getNode(5));
            t1->addNode(eg.getNode(6));
            
            t1->addEdge(1, 2, 1);
            t1->addEdge(2, 5, 1);
            t1->addEdge(5, 6, 1);
            t1->addEdge(6, 1, 1);
            eg.addTimeslice(t1);

            Graph<DataType_> * t2 = new Graph<DataType_>(4, 2);
            t2->addNode(eg.getNode(1));
            t2->addNode(eg.getNode(2));
            t2->addNode(eg.getNode(3));
            eg.getNode(6)->setWeight(3);
            t2->addNode(eg.getNode(6));
            
            t2->addEdge(1,2,1);
            t2->addEdge(1,3,1);
            t2->addEdge(1,6,2);
            eg.addTimeslice(t2);   

            Graph<DataType_> * t3 = new Graph<DataType_>(3, 2);   
            t3->addNode(eg.getNode(1));
            t3->addNode(eg.getNode(3));
            t3->addNode(eg.getNode(4));
            t3->addEdge(1,4,1);
            t3->addEdge(3,4,1);
            t3->addEdge(1,3,1);
            eg.addTimeslice(t3);            
            
            std::cout << "coordinates eg: " << *eg.coordinates();
            std::cout << "edge matrix\n" << *eg.edges();

            std::cout << "\nCalculate Position";
            //Positions<Tag_, DataType_, methods::WeightedKamadaKawai> positions(&gl_globals::graph, (DataType_)1);
            
            Engine::setTestCase(eg, new Positions<Tag_, DataType_, GraphTag_>(eg, (DataType_)1), 5);
            

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
//EngineEvolvingGraphTest<tags::CPU::SSE, float, methods::WeightedKamadaKawai> engine_test_double("wkk double");
EngineEvolvingGraphTest<tags::CPU::SSE, float, methods::WeightedFruchtermanReingold> engine_test_double("wkk double");
