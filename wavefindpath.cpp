/*!
 * @file main cpp file for finding intersection of two graphs
 * @author Malsim Orlov
 * @date 06.01.2020
 */
#include <iostream>
#include <vector>
#include "cpp/sc_addr.hpp"
#include "cpp/sc_memory.hpp"
#include "cpp/sc_iterator.hpp"
#include "utils.h"


ScAddr graph1, graph2, rrel_arcs, rrel_nodes;

/*!
 * @brief Define if vertex in set
 * @param context Sc memory functions
 * @param vertex Vertex we want to find in set
 * @param set Set where to search vertex
 * @return True if vertex is in set and False if not
 */
bool find_vertex_in_set (const std::unique_ptr<ScMemoryContext>& context, ScAddr vertex, ScAddr set)
{
    ScIterator3Ptr location = context->Iterator3(set, ScType::EdgeAccessConstPosPerm,ScType(0)); // Create iterator to find elements belongs to the set

    while (location->Next()) // Getting every element in set
    {
        ScAddr element = location->Get(2);

        if (vertex == element)
            return true;
    }

    return false;
}

/*!
 * @brief Writes the begin of the edge to v1, the end of the edge to v2
 * @param context Sc memory functions
 * @param edge The edge from we want to get vertices
 * @param v1 The first (Source) vertex the begin of the edge will be written
 * @param v2 The second (Target) vertex the end of the edge will be written
 */
void get_edge_vertexes (const std::unique_ptr<ScMemoryContext>& context, ScAddr edge, ScAddr &v1, ScAddr &v2)
{
    v1 = context->GetEdgeSource(edge);
    v2 = context->GetEdgeTarget(edge);
}

/*!
 * @brief Function to print Graph to console
 * @param context Sc memory functions
 * @param graph The Graph we want to print
 */
void print_graph (const std::unique_ptr<ScMemoryContext>& context, ScAddr graph)
{
    ScAddr arcs, nodes, v1, v2, printed_vertex;
    printed_vertex = context->CreateNode(ScType::Const); // Creating node (set) to store printed vertex

    printEl(context, graph); // printing name (identifier) of Graph
    std::cout << std::endl << "----------------------" << std::endl;

    ScIterator5Ptr it = context->Iterator5(graph, ScType::EdgeAccessConstPosPerm,ScType(0), ScType::EdgeAccessConstPosPerm, rrel_arcs); // Iterator to find arcs set

    if (it->Next())
    {
        arcs = it->Get(2);

        ScIterator3Ptr arcs_it = context->Iterator3(arcs, ScType::EdgeAccessConstPosPerm,ScType(0)); // Iterator to arcs

        while (arcs_it->Next()) // Getting every arc from arcs set
        {
            ScAddr t_arc = arcs_it->Get(2);

            get_edge_vertexes(context, t_arc, v1, v2); // Getting vertex of the arc
            printEl(context, t_arc); // Printing edge like (v1 --> v2)
            std::cout << std::endl;

            if (!find_vertex_in_set(context,v1, printed_vertex)) // Adding v1 and v2 to printed vertices to avoid double printing same vertices
                context->CreateEdge(ScType::EdgeAccessConstPosPerm, printed_vertex, v1);
            if (!find_vertex_in_set(context,v2, printed_vertex))
                context->CreateEdge(ScType::EdgeAccessConstPosPerm, printed_vertex, v2);
        }
    }

    it = context->Iterator5(graph, ScType::EdgeAccessConstPosPerm,ScType(0), ScType::EdgeAccessConstPosPerm, rrel_nodes); // Getting set of vertices

    if (it->Next())
    {
        nodes = it->Get(2);
        ScIterator3Ptr nodes_it = context->Iterator3(nodes, ScType::EdgeAccessConstPosPerm,ScType(0));

        while (nodes_it->Next()) // Getting all vertices from set
        {
            ScAddr t_node = nodes_it->Get(2);

            if (!find_vertex_in_set(context, t_node, printed_vertex))
            {
                printEl(context, t_node); // Printing vertices that haven't been printed earlier
                std::cout << std::endl;
            }
        }
    }
    std::cout << "----------------------" << std::endl;
}

/*!
 * @brief Function that defines intersection of graphs
 * @param context Sc memory functions
 * @param graph1 The first Graph of intersection
 * @param graph2 The second Graph of intersection
 * @param result_graph_name String identifier for result Graph
 * @return result Sc address of graph that is an intersection of graph1 and graph2
 */
ScAddr graph_intersection (const std::unique_ptr<ScMemoryContext>& context, ScAddr graph1, ScAddr graph2, std::string result_graph_name)
{
    /*! Creating graph structure (node for graph, for nodes set, arcs set and relations between them) */
    ScAddr result = context->CreateNode(ScType::Const);
    context->HelperSetSystemIdtf(result_graph_name, result);
    ScAddr nodes = context->CreateNode(ScType::Const);
    context->HelperSetSystemIdtf("Result_nodes", nodes);
    ScAddr arcs = context->CreateNode(ScType::Const);
    context->HelperSetSystemIdtf("Result_arcs", arcs);

    ScAddr arc_to_nodes = context->CreateEdge(ScType::EdgeAccessConstPosPerm, result, nodes);
    ScAddr arc_to_arcs = context->CreateEdge(ScType::EdgeAccessConstPosPerm, result, arcs);

    context->CreateEdge(ScType::EdgeAccessConstPosPerm, rrel_arcs, arc_to_arcs);
    context->CreateEdge(ScType::EdgeAccessConstPosPerm, rrel_nodes, arc_to_nodes);

    /*! creating Iterators to get sets for nodes and arcs of graph1 and graph2 */
    ScIterator5Ptr arcs_five1 = context->Iterator5(graph1, ScType::EdgeAccessConstPosPerm,ScType(0), ScType::EdgeAccessConstPosPerm, rrel_arcs);
    ScIterator5Ptr nodes_five1 = context->Iterator5(graph1, ScType::EdgeAccessConstPosPerm,ScType(0), ScType::EdgeAccessConstPosPerm, rrel_nodes);
    ScIterator5Ptr arcs_five2 = context->Iterator5(graph2, ScType::EdgeAccessConstPosPerm,ScType(0), ScType::EdgeAccessConstPosPerm, rrel_arcs);
    ScIterator5Ptr nodes_five2 = context->Iterator5(graph2, ScType::EdgeAccessConstPosPerm,ScType(0), ScType::EdgeAccessConstPosPerm, rrel_nodes);

    if (nodes_five1->Next() and nodes_five2->Next())
    {
        ScAddr nodes2 = nodes_five2->Get(2);
        ScAddr nodes1 = nodes_five1->Get(2);

        ScIterator3Ptr nodes_it1 = context->Iterator3(nodes1, ScType::EdgeAccessConstPosPerm,ScType(0));
        ScIterator3Ptr nodes_it2;

        while (nodes_it1->Next()) // Getting every vertex of the first graph
        {
            nodes_it2 = context->Iterator3(nodes2, ScType::EdgeAccessConstPosPerm,ScType(0));

            while (nodes_it2->Next()) // Getting every vertex of the second graph
            {
                if (nodes_it1->Get(2) == nodes_it2->Get(2)) // identifier of the first graph vertex and of the second graph vertex are the same
                {
                    ScAddr new_vertex = nodes_it1->Get(2);
                    context->CreateEdge(ScType::EdgeAccessConstPosPerm, nodes, new_vertex); // Adding this vertex to the result graph nodes set
                }
            }
        }

        if (arcs_five1->Next() and arcs_five2->Next())
        {
            ScAddr arcs1 = arcs_five1->Get(2);
            ScAddr arcs2 = arcs_five2->Get(2);

            ScIterator3Ptr arcs_it1 = context->Iterator3(arcs1, ScType::EdgeAccessConstPosPerm, ScType(0));
            ScIterator3Ptr arcs_it2;

            while (arcs_it1->Next()) // Getting every arc of the first Graph
            {
                arcs_it2 = context->Iterator3(arcs2, ScType::EdgeAccessConstPosPerm, ScType(0));
                ScAddr temp_arc1 = arcs_it1->Get(2);

                while (arcs_it2->Next()) // Getting every arc of the second Graph
                {
                    ScAddr temp_arc2 = arcs_it2->Get(2);
                    ScAddr v1, v2;
                    get_edge_vertexes(context, temp_arc1, v1, v2); // Getting begin and end of the first graph edge
                    ScAddr v3, v4;
                    get_edge_vertexes(context, temp_arc2, v3, v4); // Getting begin and end of the second graph edge

                    if (v1 == v3 and v2 == v4) // edges are connected with the same vertices
                    {
                        context->CreateEdge(ScType::EdgeAccessConstPosPerm, arcs, temp_arc1); // Adding the edge to the result graph arcs set
                    }
                }
            }
        }

    }

    return result;
}

/*! @brief function call the intersection function and printing result
 * @param context Sc memory functions
 * @param number_test String to count and differ tests
 * @param graph_name1 String to held identifier of the first graph
 * @param graph_name2 String to held identifier of the second graph
 */
void run_test (const std::unique_ptr<ScMemoryContext>& context, std::string number_test, std::string graph_name1, std::string graph_name2)
{
    std::string result_graph_name = "Graph (" + graph_name1 + " \\cap " + graph_name2 + ")"; // Result graph name consists of source graph names and the intersection tag (\cup)

    std::cout << "Test" << number_test << ":" << std::endl;

    /*! Initializing graphs from sc-memory */
    graph1 = context->HelperResolveSystemIdtf(graph_name1);
    graph2 = context->HelperResolveSystemIdtf(graph_name2);

    rrel_arcs = context->HelperResolveSystemIdtf("rrel_arcs");
    rrel_nodes = context->HelperResolveSystemIdtf("rrel_nodes");

    /*! Printing  source graphs and result graph*/
    std::cout << "Graph1: ";
    print_graph(context, graph1);

    std::cout << std::endl;

    std::cout << "Graph2: ";
    print_graph(context, graph2);

    ScAddr result_graph = graph_intersection(context, graph1, graph2, result_graph_name);

    std::cout <<  std::endl << "Result Graph: ";
    print_graph(context, result_graph);
    std::cout << std::endl << std::endl;
}

/*! Configure params and running tests */
int main()
{
    sc_memory_params params;

    sc_memory_params_clear(&params);
    params.repo_path = "/home/mqsm/cursuch/project/ostis-web-platform/kb.bin";
    params.config_file = "/home/mqsm/cursuch/project/ostis-web-platform/config/sc-web.ini";
    params.ext_path = "/home/mqsm/cursuch/project/ostis-web-platform/sc-machine/bin/extensions";
    params.clear = SC_FALSE;

    ScMemory mem;
    mem.Initialize(params);

    const std::unique_ptr<ScMemoryContext> context(new ScMemoryContext(sc_access_lvl_make_max,"example"));

    run_test(context,"0", "G0", "G1");
    run_test(context,"1", "G1", "G2");
    run_test(context,"2", "G2", "G2");
    run_test(context,"3", "G3", "G4");
    run_test(context,"4", "G4", "G5");


    std::cout << "The end" << std::endl;

    mem.Shutdown(true);

    return 0;
}

