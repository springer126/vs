#ifndef _VESSELGRAPH_H_
#define _VESSELGRAPH_H_

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <map>

#include "mitkPointSet.h"
#include "mitkColorSequenceRainbow.h"
#include <mitkImage.h>
#include <itkImage.h>


struct vertex_Xcoordinate_t {
	typedef boost::vertex_property_tag kind;
};
struct vertex_Ycoordinate_t {
	typedef boost::vertex_property_tag kind;
};
struct vertex_Zcoordinate_t {
	typedef boost::vertex_property_tag kind;
};
struct vertex_radius_t {
	typedef boost::vertex_property_tag kind;
};
struct edge_meanradius_t {
	typedef boost::edge_property_tag kind;
};
struct edge_distance_t {
	typedef boost::edge_property_tag kind;
};
struct edge_length_t {
	typedef boost::edge_property_tag kind;
};	
struct edge_treeno_t {
	typedef boost::edge_property_tag kind;
};
struct edge_angle_t {
	typedef boost::edge_property_tag kind;
};
struct edge_visited_t {
	typedef boost::edge_property_tag kind;
};

typedef boost::property<vertex_Xcoordinate_t, double> XCoordProperty;
typedef boost::property<vertex_Ycoordinate_t,double, XCoordProperty> XYCoordProperty;
typedef boost::property<vertex_Zcoordinate_t,double, XYCoordProperty> XYZCoordProperty;
typedef boost::property<vertex_radius_t,double, XYZCoordProperty> VertexProperty;
typedef boost::property<edge_visited_t,int> VisProperty;
typedef boost::property<edge_meanradius_t, double,VisProperty> MeanrProperty;
typedef boost::property<edge_distance_t,double, MeanrProperty> DisProperty;
typedef boost::property<edge_length_t,double, DisProperty> LengthProperty;
typedef boost::property<edge_treeno_t,int, LengthProperty> TreenoProperty;
typedef boost::property<edge_angle_t,double, TreenoProperty> EdgeProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperty, EdgeProperty> Graph;

typedef boost::property_map<Graph, vertex_Xcoordinate_t>::type TXCoord;
typedef boost::property_map<Graph, vertex_Ycoordinate_t>::type TYCoord;
typedef boost::property_map<Graph, vertex_Zcoordinate_t>::type TZCoord;
typedef boost::property_map<Graph, vertex_radius_t>::type TRadius;
typedef boost::property_map<Graph, edge_visited_t>::type TEdgeVisit;
typedef boost::property_map<Graph, edge_meanradius_t>::type TMeanRadius;
typedef boost::property_map<Graph, edge_distance_t>::type TDistance;
typedef boost::property_map<Graph, edge_length_t>::type TLength;
typedef boost::property_map<Graph, edge_treeno_t>::type TTreeno;
typedef boost::property_map<Graph, edge_angle_t>::type TAngle;
typedef boost::graph_traits<Graph> GraphTraits;
//

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
typedef boost::graph_traits<Graph>::in_edge_iterator InEdgeIterator;



class VesselGraph
{
	typedef std::map<int , std::list<Vertex> > VertexMap;
public:
	
	void ReadFile(const char * Path);
	unsigned int GetVertexIndex(Vertex v) const ;
	unsigned int GetEdgeSourceIndex(Edge v) const ;
	unsigned int GetEdgeTargetIndex(Edge v) const ;
	void OutputGraph();
	Vertex GetMaxRadiusEdgeNode();
	void CreateSubGraph(Vertex& v);
	bool GetSubGraph();
	void DivideVesselTree();
	mitk::Image::Pointer VoxelDivision(mitk::Image *image);
	void RadiusFilter(double threshold);
	void KMeans(unsigned int clsNumber,double theta);
	void SetSeparateArg(double arg)
	{
		this->separateArgument = arg;
	}
	void SetStartVertex(int start)
	{
		this->start = start;
	}
	void SetSubGraphCount(int arg)
	{
		this->subGraphCount = arg;
	}

	int GetSubGraphCount()
	{
		return this->subGraphCount;
	}
	
	//距离归类
	template <typename TPixel, unsigned int VImageDimension>
	void NNSA(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Image::Pointer *pointer);

	//膨胀
	template <typename TPixel, unsigned int VImageDimension>
	void AnotherNNSA(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Image::Pointer *pointer);

	void Reset();
	VesselGraph(const char * Path);
	~VesselGraph();
public:
	Graph g;
	std::vector<Vertex> vs;
	std::vector<Edge> es;

	int subGraphCount;
	int start;
	
	//QDiVertex[2-subGraphCount]存放各个分支节点集合,QDivVertex[1]存放一级分支
	//交互式选定分支时，起始的二级分支为QDivVertex[2]
	VertexMap QDiVertex;
	std::queue<Vertex> vq;
	double separateArgument;
	
	TMeanRadius edgeMrMap;
	TEdgeVisit edgeVisitMap;
	TXCoord txMap;
	TYCoord tyMap;
	TZCoord tzMap;
	mitk::Point3D point3D;
	mitk::ColorSequenceRainbow m_RainbowColor;
};

#endif