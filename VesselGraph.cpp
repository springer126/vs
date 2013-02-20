#pragma comment(lib,"boost_graph-vc80-1_38.lib")

#include "VesselGraph.h"
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map.hpp>
#include <mitkPointSet.h>
#include "mitkImageAccessByItk.h"
#include "mitkITKImageImport.h"
#include "mitkPointSet.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "mitkITKImageImport.h"
#include "mitkDataStorage.h"
#include "mitkDataTreeHelper.h"
#include "QmitkStdMultiWidget.h"
#include "mitkProperties.h"
#include <qmessagebox.h> 
#include <qdatetime.h> 
#include <limits>
#include <math.h>
#include <stack>
#include <fstream>
#include <time.h>
#include <algorithm>

using namespace boost;
typedef boost::property_map<Graph, vertex_index_t>::type TVertexIndex;
typedef boost::property_map<Graph, edge_index_t>::type TEdgeIndex;

int id;
bool isCls(int cls)
{
	if(cls==id)
		return true;
	else
		return false;
};

template <class Graph> 
struct exercise_vertex
{
	exercise_vertex(Graph& g_) : g(g_) {}
	Graph& g;
	void operator()(const Vertex& v) const
	{
		TVertexIndex index = get(vertex_index, g);
		TRadius radius = get(vertex_radius_t(), g);
		std::cout <<"Index:"<< index[v]<<" Radius:"<<radius[v]<<" Out-Edges: ";
		OutEdgeIterator out_i, out_end;
		Edge e;
		for (boost::tie(out_i, out_end) = out_edges(v, g); 
			out_i != out_end; ++out_i) {
				e = *out_i;
				Vertex src = source(e, g), targ = target(e, g);
				std::cout << "(" << index[src] << "," 
					<< index[targ] << ") ";
		}
		std::cout << "  In-edges: ";
		InEdgeIterator in_i, in_end;
		for (tie(in_i, in_end) = in_edges(v,g); 
			in_i != in_end; ++in_i) {
				e = *in_i;
				Vertex src = source(e, g), targ = target(e, g);
				std::cout << "(" << index[src] << "," << index[targ] << ") ";
		}
		std::cout << "  Adjacent Vertexs:";
		AdjacencyIterator ad_i,ad_end;
		for(tie(ad_i,ad_end) = adjacent_vertices(v,g);ad_i!=ad_end;ad_i++)
		{
			std::cout <<" "<< index[*ad_i];	
		}
		std::cout << std::endl;
	}
};

VesselGraph::VesselGraph(const char * Path):g(0),subGraphCount(0),separateArgument(0.6),start(-1)
{
	ReadFile(Path);
	edgeMrMap = get(edge_meanradius_t(),g);
	edgeVisitMap = get(edge_visited_t(),g);
	txMap = get(vertex_Xcoordinate_t(),g);
	tyMap = get(vertex_Ycoordinate_t(),g);
	tzMap = get(vertex_Zcoordinate_t(),g);

	EdgeIterator e_i,e_end;
	for (boost::tie(e_i,e_end)=edges(g);e_i!=e_end;++e_i)
	{
		edgeVisitMap[*e_i] = 0;
	}
}

VesselGraph::~VesselGraph()
{
	
}

void VesselGraph::ReadFile(const char * Path)
{
	std::ifstream infile;
	infile.open( Path, std::ios::in);
	/** Read graphml */
	boost::dynamic_properties dp;
	dp.property("XCoord", get(vertex_Xcoordinate_t(), g));
	dp.property("YCoord", get(vertex_Ycoordinate_t(), g));
	dp.property("ZCoord", get(vertex_Zcoordinate_t(), g));
	dp.property("Radius", get(vertex_radius_t(), g));
	dp.property("MeanRadius", get(edge_meanradius_t(), g));
	dp.property("Distance", get(edge_distance_t(), g));
	dp.property("Length", get(edge_length_t(), g));	
	dp.property("Treeno", get(edge_treeno_t(), g));
	dp.property("Angle", get(edge_angle_t(), g));
	read_graphml(infile, g, dp); 
	

	VertexIterator v_i,v_end;
	EdgeIterator e_i,e_end;
	vs.clear();
	es.clear();
	for (boost::tie(v_i,v_end) = vertices(g);v_i!=v_end;v_i++)
	{
		vs.push_back(*v_i);
	}
	for (boost::tie(e_i,e_end) = edges(g);e_i!=e_end;++e_i)
	{
		es.push_back(*e_i);
	}
	return;
}

Vertex VesselGraph::GetMaxRadiusEdgeNode()
{
	assert(num_vertices(g)!=0);
	assert(num_edges(g)!=0);

	TVertexIndex vertexIndexMap = get(vertex_index,g);
	EdgeIterator e_i,e_end;
	Edge e;
	e = *(edges(g).first);
	for (boost::tie(e_i,e_end) = edges(g);e_i!=e_end;++e_i)
	{
		if(edgeMrMap[*e_i]>edgeMrMap[e] && edgeVisitMap[*e_i]==0)
			e = *e_i;
	}
	std::cout << "Max Edge Radius:"<<edgeMrMap[e] <<" Source:"<<vertexIndexMap[source(e,g)]<<" Target:"<<vertexIndexMap[target(e,g)]<< std::endl;
	Vertex v = source(e,g);
	return v;
	//return source(e,g);
}

//???????????????,?????????????????separateArgument?,???????push?vq?????
void VesselGraph::CreateSubGraph(Vertex& v)
{
	assert(num_vertices(g)!=0);
	assert(num_edges(g)!=0);
    TVertexIndex vertexIndexMap = get(vertex_index,g);

	//Reset();
	++subGraphCount;
	QDiVertex[subGraphCount].push_back(v);
	Edge e;
	bool finish = false;
	std::list<Vertex>::iterator v_i;
	while (true && !finish)
	{
		finish = true;
		for (v_i=QDiVertex[subGraphCount].begin();v_i!=QDiVertex[subGraphCount].end();++v_i)
		{
			OutEdgeIterator e_i,e_end;
			for(boost::tie(e_i,e_end) = out_edges(*v_i,g);e_i!=e_end;++e_i)
			{
				if(edgeVisitMap[*e_i]==0)
				{
					finish = false;	
					e = *e_i;	
					break;
				}
			}
		}
		if(finish)
		{
			break;	
		}

		for (v_i=QDiVertex[subGraphCount].begin();v_i!=QDiVertex[subGraphCount].end();++v_i)
		{
			OutEdgeIterator e_i,e_end;
			for(boost::tie(e_i,e_end) = out_edges(*v_i,g);e_i!=e_end;++e_i)
			{
				if(edgeVisitMap[*e_i]==0)
				{
					if(edgeMrMap[*e_i]>edgeMrMap[e])
						e = *e_i;
				}
			}
		}
		if(source(e,g)!=v)
		{
			Edge parentEdge = *(in_edges(source(e,g),g).first);
			if(edgeMrMap[parentEdge]*separateArgument < edgeMrMap[e])
			{
				    //std::cout << "Edge Vertex Index: (" << vertexIndexMap[source(e,g)] <<","<<vertexIndexMap[target(e,g)]<<")"<< std::endl;
					QDiVertex[subGraphCount].push_back(target(e,g));
					edgeVisitMap[e] = subGraphCount;
			}
			else
			{
				edgeVisitMap[e] = -1;
				if(out_degree(target(e,g),g)!=0)
					vq.push(target(e,g));
				else
				{
					edgeVisitMap[e] = subGraphCount;
					//std::cout << "Edge Vertex Index: (" << vertexIndexMap[source(e,g)] <<","<<vertexIndexMap[target(e,g)]<<")"<< std::endl;
					QDiVertex[subGraphCount].push_back(target(e,g));
				}
			}
		}
		else
		{
			QDiVertex[subGraphCount].push_back(target(e,g));
			edgeVisitMap[e] = subGraphCount;
		}
	}

	/*mitk::PointSet::Pointer points = mitk::PointSet::New();
	int index = 0;
	std::cout <<"No."<<subGraphCount<<" SubGraph ,Vertex Number: " <<QDiVertex[subGraphCount].size()<< std::endl;
	for (v_i=QDiVertex[subGraphCount].begin();v_i!=QDiVertex[subGraphCount].end();++v_i)
	{
		std::cout <<vertexIndexMap[*v_i]<<" ";
		point3D[0] = txMap[*v_i];
		point3D[1] = tyMap[*v_i];
		point3D[2] = tzMap[*v_i];
		points->GetPointSet()->GetPoints()->InsertElement( index++, point3D );
	}

	std::cout << std::endl;
	mitk::DataTreeNode::Pointer pNode = mitk::DataTreeNode::New();
	pNode->SetData( points );
	char name[20];
	sprintf(name,"No.%d SubGraph",subGraphCount);
	pNode->SetProperty("name", mitk::StringProperty::New(name));
	pNode->SetProperty("layer", mitk::IntProperty::New(1));
	if (subGraphCount%8==0)
	{
		m_RainbowColor.GoToBegin();
	}
	mitk::Color color = m_RainbowColor.GetNextColor();
	pNode->SetProperty("color", mitk::ColorProperty::New(color.GetRed(),color.GetGreen(),color.GetBlue()));
	pNode->SetProperty("pointsize", mitk::FloatProperty::New(1));
	mitk::DataStorage::GetInstance()->Add( pNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();*/
}


bool VesselGraph::GetSubGraph()
{
	//Reset();
	if(start==-1)
		return false;

	Vertex v = vs.at(start);
	std::cout << txMap[v]<<","<<tyMap[v]<<","<<tzMap[v]<<std::endl;
	++subGraphCount;
	
	std::queue<Vertex> vq;
	QDiVertex[subGraphCount].clear();
	QDiVertex[subGraphCount].push_back(v);
	vq.push(v);
	
	OutEdgeIterator e_beg,e_end;
	while (!vq.empty())
	{
		v = vq.front();
		vq.pop();
		for (tie(e_beg,e_end) = out_edges(v,g);e_beg!=e_end;e_beg++)
		{
			v = target(*e_beg,g);
			QDiVertex[subGraphCount].push_back(v);
			vq.push(v);
		}
	}
	start = -1;
	if(QDiVertex[subGraphCount].size()<=1)
	{
		QDiVertex[subGraphCount].clear();
		--subGraphCount;
		std::cout << "[Deprecated:1 point.]" <<std::endl;
		return false;
	}
	else
	{
		std::cout << "SubGraph: "<<subGraphCount;
		std::cout <<" Vertex Numbers: "<< QDiVertex[subGraphCount].size()<<std::endl;
		std::list<Vertex>::iterator begin;
		TVertexIndex index = get(vertex_index, g);
		for (begin = QDiVertex[subGraphCount].begin();begin!=QDiVertex[subGraphCount].end();begin++)
		{
			std::cout << index[*begin]<<" ";
		}
		std::cout << std::endl;

		return true;
	}
	
}


// ?separateArgument????????
void VesselGraph::DivideVesselTree()
{
	vq.push(this->GetMaxRadiusEdgeNode());
	while(!vq.empty())
	{
		Vertex v = vq.front();
		vq.pop();
		this->CreateSubGraph(v);
	}
}

mitk::Image::Pointer VesselGraph::VoxelDivision(mitk::Image *image)
{
	mitk::Image::Pointer resultImage;
	//std::cout << image->GetPixelType().GetItkTypeAsString() << std::endl;
	if(subGraphCount<2 || image==NULL)
		return NULL;
	else
	{
		AccessByItk_1(image,NNSA,&resultImage);
		return resultImage;
	}
}

void VesselGraph::Reset()
{
	/*if(subGraphCount!=0)
	{
		for (int i=1;i<=subGraphCount;i++)
		{
			QDiVertex[i].clear();		
		}
	}*/
	if(subGraphCount>=2)
	{
		for (int i=1;i<=subGraphCount;i++)
		{
			QDiVertex[i].clear();
			QDiEdge[i].clear();
		}
	}
	QDiVertex.clear();
	QDiEdge.clear();
	subGraphCount = 0;
}

template <typename TPixel, unsigned int VImageDimension>
void VesselGraph::NNSA(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Image::Pointer *pointer)
{
	std::cout << "[NNSA Begin]" << std::endl;
	QTime time;
	time.start();

	int count = 0;
	typedef itk::Image<TPixel, VImageDimension> TImageType;
	typedef typename TImageType::RegionType TRegionType;
	typedef itk::Point<double,VImageDimension> PointType;
	typedef typename TImageType::SpacingType Space;
	typedef TImageType::SizeType Size;
	typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(itkImage);
	duplicator->Update();
	typename TImageType::Pointer clonedImage = duplicator->GetOutput();
	
	typedef itk::ImageRegionIterator<TImageType> RegionIteratorType;
	RegionIteratorType iterator1(itkImage, itkImage->GetLargestPossibleRegion());
	RegionIteratorType iterator2(clonedImage, clonedImage->GetLargestPossibleRegion());

	TRegionType region = itkImage->GetLargestPossibleRegion();
	Size size = region.GetSize();
	std::cout <<"Size:" << size<< std::endl;
	Space space;
	space = itkImage->GetSpacing();
	std::cout <<"Space :" << space << std::endl;

	for (iterator1.GoToBegin(),iterator2.Begin(); !iterator1.IsAtEnd(); ++iterator1,++iterator2)
	{
		count++;
		if(count%100000==0)
		{
			std::cout <<"Solving "<< count/100000 <<"00000 : " <<count/6360300.0<<std::endl;
		}

		if (iterator1.Get()==-10000)
		{
			iterator2.Set(-10000);
			continue;
		}

		PointType point;
		float* dis = new float[subGraphCount+1];
		itkImage->TransformIndexToPhysicalPoint( iterator1.GetIndex(), point );
		std::list<Vertex>::iterator v_i;
		for (int i=2;i<=subGraphCount;i++)
		{
			float distance = FLT_MAX,tmp;
			for (v_i=QDiVertex[i].begin();v_i!=QDiVertex[i].end();++v_i)
			{
				tmp = sqrt(pow(txMap[*v_i]-point[0],2)+pow(tyMap[*v_i]-point[1],2)+pow(tzMap[*v_i]-point[2],2));
				if(tmp<distance)
					distance = tmp;
			}
			dis[i] = distance;
		}

		float minDist = dis[2];
		int clsId = 2;
		for(int i=2;i<=subGraphCount;i++)
		{
			if(dis[i]<minDist)
			{
				minDist = dis[i];	
				clsId = i;
			}
		}
		//??8????:20-90
		iterator2.Set((TPixel)clsId*10);
		//iterator2.Set(100);
	}

	
	*pointer = mitk::ImportItkImage( clonedImage );
	mitk::Image::Pointer result = mitk::ImportItkImage( clonedImage );
	mitk::DataTreeNode::Pointer pNode = mitk::DataTreeNode::New();
	pNode->SetData( result );
	pNode->SetProperty("name", mitk::StringProperty::New("Liver Image"));
	pNode->SetProperty("layer", mitk::IntProperty::New(1));
	pNode->SetProperty("color", mitk::ColorProperty::New(1.0,0.0,0.0));
	mitk::DataStorage::GetInstance()->Add( pNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	
	char buffer[64];
	sprintf(buffer, "Time elapsed: %.03fs", time.elapsed()/1000.0);
	printf("%s\n", buffer);
	QMessageBox::information(NULL, "Voxel Classify", QString(buffer), QMessageBox::Ok);

	std::cout << "[NNSA End]" << std::endl;
}



template <typename TPixel, unsigned int VImageDimension>
void VesselGraph::AnotherNNSA(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Image::Pointer *pointer)
{
	std::cout << "[NNSA1 Begin]" << std::endl;
	QTime time;
	time.start();
	typedef itk::Image<TPixel, VImageDimension> TImageType;
	typedef typename TImageType::RegionType TRegionType;
	typedef itk::Point<double,VImageDimension> PointType;
	typedef TImageType::SizeType Size;
	typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	typedef itk::ImageRegionIterator<TImageType> RegionIteratorType;
	
	typedef itk::ConstantBoundaryCondition< TImageType > ConstBoundaryConditionType;
	typedef itk::NeighborhoodIterator<TImageType, ConstBoundaryConditionType> NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType radius;
	ConstBoundaryConditionType boundaryCondition;
	boundaryCondition.SetConstant( 0 );
	radius.Fill(1);
	

	NeighborhoodIteratorType it;
	bool isFinished = false;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(itkImage);
	duplicator->Update();
	typename TImageType::Pointer clonedImage = duplicator->GetOutput();

	RegionIteratorType iterator1(itkImage, itkImage->GetLargestPossibleRegion());
	RegionIteratorType iterator2(clonedImage, clonedImage->GetLargestPossibleRegion());
	

	for (iterator1.GoToBegin(),iterator2.Begin(); !iterator1.IsAtEnd(); ++iterator1,++iterator2)
	{
		if (iterator1.Get()==-10000)
		{
			iterator2.Set(-10000);
		}
		else
		{
			iterator2.Set(0);	
		}
	}
	PointType point,ipoint;
	typename TImageType::IndexType index;
	std::list<Vertex>::iterator v_i;
	for (int i=2;i<=subGraphCount;i++)
	{
		for (v_i=QDiVertex[i].begin();v_i!=QDiVertex[i].end();++v_i)
		{
			point[0]  = txMap[*v_i];
			point[1]  = tyMap[*v_i];
			point[2]  = tzMap[*v_i];
			clonedImage->TransformPhysicalPointToIndex(point,index);
			clonedImage->SetPixel(index,i*10);
		}
	}

	it.Initialize( radius, clonedImage, clonedImage->GetLargestPossibleRegion() );
	it.SetBoundaryCondition( boundaryCondition );  
	std::cout << "Size: "<<it.Size()<< std::endl;
	long int step = 10;
	long int n_step = 0;
	while(!isFinished)
	{
		n_step++;
		std::cout << "Step : "<< n_step << std::endl;
		isFinished = true;
		for (iterator2.GoToBegin(); !iterator2.IsAtEnd(); ++iterator2)
		{
			if (iterator2.Get()==0)
			{
				it.SetLocation( iterator2.GetIndex());
				isFinished = false;
				float distance  = FLT_MAX,temp;
				//clonedImage->TransformIndexToPhysicalPoint( iterator2.GetIndex(), point );
				typename NeighborhoodIteratorType::OffsetType offset;
				for (unsigned int i = 0; i < it.Size() && i!=it.GetCenterNeighborhoodIndex(); ++i)
				{    
					if(it.GetPixel(i)>0)
					{
						//itkImage->TransformIndexToPhysicalPoint( 
							//it.GetNeighborhoodIndex(it.ComputeInternalIndex(i)), ipoint );
						//temp = sqrt(pow(ipoint[0]-point[0],2)+pow(ipoint[1]-point[1],2)+pow(ipoint[2]-point[2],2));
						offset = it.GetOffset(i);
						temp = 0;
						for (int k=0;k<offset.GetOffsetDimension();k++)
						{
							temp += pow((double)offset[k],2);
						}
						if(temp < distance)
						{
							distance = temp;
							iterator2.Set(it.GetPixel(i));
						}
					} 
				}
			}
		}
	}

	*pointer = mitk::ImportItkImage( clonedImage );
	mitk::Image::Pointer result = mitk::ImportItkImage( clonedImage );
	mitk::DataTreeNode::Pointer pNode = mitk::DataTreeNode::New();
	pNode->SetData( result );
	pNode->SetProperty("name", mitk::StringProperty::New("Liver Image"));
	pNode->SetProperty("layer", mitk::IntProperty::New(1));
	pNode->SetProperty("color", mitk::ColorProperty::New(1.0,0.0,0.0));
	mitk::DataStorage::GetInstance()->Add( pNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();

	char buffer[64];
	sprintf(buffer, "Time elapsed: %.03fs", time.elapsed()/1000.0);
	printf("%s\n", buffer);
	QMessageBox::information(NULL, "Voxel Classify", QString(buffer), QMessageBox::Ok);
	std::cout << "[NNSA1 End]" << std::endl;
}


void VesselGraph::OutputGraph()
{
	std::cout << "==============Graph:=============" << std::endl;
	std::for_each(vertices(g).first,vertices(g).second,exercise_vertex<Graph>(g));
	std::cout <<"==============Graph============="<< std::endl;
}

int VesselGraph::GetVertexIndex(Vertex v) const 
{
	TVertexIndex vertexIndexMap = get(vertex_index,g);
	return vertexIndexMap[v];
}



int VesselGraph::GetEdgeSourceIndex(Edge e) const 
{
	 return GetVertexIndex(source(e,g));
}

int VesselGraph::GetEdgeTargetIndex(Edge e) const
{
	return GetVertexIndex(target(e,g));
}

void VesselGraph::RadiusFilter(double lowRadius,int lowNodes)
{
	assert((vs.size()>0));
	assert((es.size()>0));
	TVertexIndex vertexIndexMap = get(vertex_index,g);
	std::vector<Edge>::iterator e_beg;
	std::vector<Vertex>::iterator v_beg,beg;
	QDiVertex[1].clear();
	//erase edge whose radius more than threshold and its source vertex
	//if sink point is leaf then delete it
	for (e_beg = es.begin();e_beg!=es.end();e_beg++)
	{
		if(edgeMrMap[*e_beg]<lowRadius)
		{
			res.push_back(*e_beg);
		}
	}

	for (e_beg = res.begin();e_beg!=res.end();e_beg++)
	{
		for (v_beg = rvs.begin();v_beg!=rvs.end() && vertexIndexMap[*v_beg] !=source(*e_beg,g);v_beg++);
		if(v_beg==rvs.end())
		{
			for (beg = vs.begin();beg!=vs.end();beg++)
			{
				if(vertexIndexMap[*beg]==source(*e_beg,g))
				{
					rvs.push_back(*beg);
					break;
				}
			}
		}

		for (v_beg = rvs.begin();v_beg!=rvs.end() && vertexIndexMap[*v_beg] !=target(*e_beg,g);v_beg++);
		if(v_beg==rvs.end())
		{
			for (beg = vs.begin();beg!=vs.end();beg++)
			{
				if(vertexIndexMap[*beg]==target(*e_beg,g))
				{
					rvs.push_back(*beg);
					break;
				}
			}
		}
	}


	/*for (v_beg = rvs.begin();v_beg!=rvs.end();v_beg++)
	{
		for (e_beg = es.begin();e_beg!=es.end();e_beg++)
		{
			if(target(*e_beg,g)==vertexIndexMap[*v_beg] && find(res.begin(),res.end(),*e_beg)==res.end())
			{
				
			}
		}
	}*/

	std::cout <<"Steps 1: Vertex Number :" << rvs.size()<<" , Edge Number:" << res.size() << std::endl;
	vs.clear();

	mitk::Point3D point3D;
	mitk::PointSet::Pointer roots = mitk::PointSet::New();
	int index = 0;
	for (v_beg = rvs.begin();v_beg!=rvs.end();v_beg++)
	{
		for (e_beg = es.begin();e_beg!=es.end();e_beg++)
		{	
			std::vector<Vertex>::iterator v_iter;
			std::vector<Edge>::iterator e_iter;
			
			//??*v_beg ??????
			if(target(*e_beg,g)==vertexIndexMap[*v_beg] && find(res.begin(),res.end(),*e_beg)==res.end())
			{
				int numOfVertexs = 0;
				AdjacencyIterator ad_i,ad_end;
	
				while (!vq.empty())
				{
					vq.pop();
				}
				vq.push(*v_beg);
				while(!vq.empty())
				{
					++numOfVertexs;
					Vertex v = vq.front();
					vq.pop();
					for(tie(ad_i,ad_end) = adjacent_vertices(v,g);ad_i!=ad_end;ad_i++)
					{
						for (e_iter = res.begin();e_iter!=res.end();e_iter++)
						{
							if(GetEdgeTargetIndex(*e_iter)==vertexIndexMap[*ad_i] && GetEdgeSourceIndex(*e_iter)==vertexIndexMap[v])
								vq.push(*ad_i);
						}
						
					}
				}
				
				std::cout << "Index:"<<vertexIndexMap[*v_beg] << " , Number Of Nodes:" << numOfVertexs<<std::endl;
				//?????????lowNodes,????????(???)
				if(numOfVertexs>lowNodes)
				{
					vs.push_back(*v_beg);
					point3D[0] = txMap[*v_beg];point3D[1] = tyMap[*v_beg];point3D[2] = tzMap[*v_beg];
					roots->GetPointSet()->GetPoints()->InsertElement( index++, point3D );
				}
				else//??,?res?rvs??????
				{
					while (!vq.empty())
					{
						vq.pop();
					}
					vq.push(*v_beg);
					while(!vq.empty())
					{
						
						Vertex v = vq.front();
						for (v_iter = rvs.begin();v_iter!=rvs.end();v_iter++)
						{
							if(vertexIndexMap[*v_iter]==vertexIndexMap[v])
							{
								rvs.erase(v_iter);
								break;
							}
						}

						vq.pop();
						for(tie(ad_i,ad_end) = adjacent_vertices(v,g);ad_i!=ad_end;ad_i++)
						{
							for (e_iter = res.begin();e_iter!=res.end();e_iter++)
							{
								if(GetEdgeTargetIndex(*e_iter)==vertexIndexMap[*ad_i] && GetEdgeSourceIndex(*e_iter)==vertexIndexMap[v])
									vq.push(*ad_i);
							}
						}
						
						for (e_iter = res.begin();e_iter!=res.end();e_iter++)
						{
							if(GetEdgeSourceIndex(*e_iter)==vertexIndexMap[v])
							{
								res.erase(e_iter);
							}
						}

					}
				}
				
			}
		}
	}

	std::cout <<"Steps 2: Vertex Number :" << rvs.size()<<" , Edge Number:" << res.size() << std::endl;

	mitk::DataTreeNode::Pointer pNode = mitk::DataTreeNode::New();
	pNode->SetData( roots );
	char name[20];
	sprintf(name,"%d subtree",index);
	pNode->SetProperty("name", mitk::StringProperty::New(name));
	pNode->SetProperty("layer", mitk::IntProperty::New(1));
	pNode->SetProperty("color", mitk::ColorProperty::New(1.0,0.0,0.0));
	pNode->SetProperty("pointsize", mitk::FloatProperty::New(2));
	mitk::DataStorage::GetInstance()->Add( pNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();




	/*ofstream out_file("D:\\vertex.scale",fstream::app);
	for (i_beg2 = vs.begin();i_beg2!=vs.end();i_beg2++)
	{
		out_file << txMap[*i_beg2]<<","<<tyMap[*i_beg2]<< "," << tzMap[*i_beg2]<<std::endl;
	}

	for (i_beg = es.begin();i_beg!=es.end();i_beg++)
	{
		std::cout << "("<<source(*i_beg,g)<<","<<target(*i_beg,g)<<")"<<std::endl;
	}*/
}

/************************************************************************/
/*   K-means to classify nodes of vascular  
	 Option: subtree roots || subtree all-nodes
/************************************************************************/
double VesselGraph::KMeans(unsigned int clsNumber,double theta)
{
	assert(clsNumber>0 && theta>0);
	assert(vs.size()>0);
	int size = vs.size();
	//find seeds
	Vertex *seeds = new Vertex[clsNumber];
	double sum,*d = new double[size];
	seeds[0] = vs.at(rand()%size);
	
	for (int n_cluster = 1;n_cluster<clsNumber;n_cluster++)
	{
		sum = 0;
		for (int i = 0;i<size;i++)
		{
			float dis = FLT_MAX;
			for (int j=0;j<n_cluster;j++)
			{
				//std::cout << "Cluster index:" << n_cluster << ",Vertex index:" << i << "Current Cluster index:" << j << std::endl;
				float temp = sqrt(pow(txMap[seeds[j]]-txMap[vs.at(i)],2)
					+pow(tyMap[seeds[j]]-tyMap[vs.at(i)],2)+pow(tzMap[seeds[j]]-tzMap[vs.at(i)],2));
				if(temp<dis)
				{
					dis = temp;
				}
			}
			d[i] = dis; 
			sum += dis;
		}

		sum = sum*rand()/(RAND_MAX-1);
		for (int k=0;k<size;k++)
		{
			if((sum-=d[k])>=0)continue;
			seeds[n_cluster] = vs.at(k);
			break;
		}
	}

	int iter = 0;
	double **centerPosition = new double*[clsNumber] ,**curPosition = new double*[clsNumber];
	
	for (int i=0;i<clsNumber;i++)
	{
		centerPosition[i]= new double[3];
		curPosition[i] = new double[3];
	}
	double scope[3][2];
	scope[0][0] = scope[0][1] = txMap[vs.at(0)];
	scope[1][0] = scope[1][1] = tyMap[vs.at(0)];
	scope[2][0] = scope[2][1] = tzMap[vs.at(0)];
	int *cls = new int[vs.size()];
	std::vector<Vertex>::iterator i_beg;
	
	//memset(centerPosition,0,sizeof(centerPosition));
	for (i_beg = vs.begin();i_beg!=vs.end();i_beg++)
	{
		if(txMap[*i_beg]>scope[0][0])
			scope[0][0] = txMap[*i_beg];
		if(txMap[*i_beg]<scope[0][1])
			scope[0][1] = txMap[*i_beg];
		if(tyMap[*i_beg]>scope[1][0])
			scope[1][0] = tyMap[*i_beg];
		if(tyMap[*i_beg]<scope[1][1])
			scope[1][1] = tyMap[*i_beg];
		if(tzMap[*i_beg]>scope[2][0])
			scope[2][0] = tzMap[*i_beg];
		if(tzMap[*i_beg]<scope[2][1])
			scope[2][1] = tzMap[*i_beg];
	}
	/*std::cout << scope[0][0]<<" "<< scope[0][1]<< " "<<scope[1][0]<< 
		" "<<scope[1][1]<<" "<< scope[2][0]<<" "<< scope[2][1]<<std::endl;*/
	/*srand( (unsigned)time( NULL ) );
	double ratio;
	for (int k=0;k<clsNumber;k++)
	{
	    ratio = ((k+0.0)/clsNumber)*((rand()%10)/10.0)+0.1;
		std::cout << "Ratio:" << ratio << std::endl; 
		*(*(centerPosition+k)+0) = ratio*(scope[0][0]-scope[0][1])+scope[0][1];
		*(*(centerPosition+k)+1) = ratio*(scope[1][0]-scope[1][1])+scope[1][1];
		*(*(centerPosition+k)+2) = ratio*(scope[2][0]-scope[2][1])+scope[2][1];
	}*/
	
	//double ratio[8][3] = {1/4.0,1/4.0,1/4.0, 1/4.0,3/4.0,1/4.0, 3/4.0,1/4.0,1/4.0, 3/4.0,3/4.0,1/4.0,
		//		1/4.0,1/4.0,3/4.0, 3/4.0,1/4.0,3/4.0,1/4.0,3/4.0,3/4.0, 3/4.0,3/4.0,3/4.0};
	
	//for (int i=0;i<clsNumber;i++)
	//{
	//	centerPosition[i][0] = scope[0][1]+(scope[0][0]-scope[0][1])*ratio[i][0];
	//	centerPosition[i][1] = scope[1][1]+(scope[1][0]-scope[1][1])*ratio[i][1];
	//	centerPosition[i][2] = scope[2][1]+(scope[2][0]-scope[2][1])*ratio[i][2];
	//}
	std::cout << "Center Positions:";
	for (int i=0;i<clsNumber;i++)
	{
		centerPosition[i][0] = txMap[seeds[i]];
		centerPosition[i][1] = tyMap[seeds[i]];
		centerPosition[i][2] = tzMap[seeds[i]];
		std::cout << txMap[seeds[i]]<<tyMap[seeds[i]]<<tzMap[seeds[i]]<<std::endl;
	}
	
	double thresh = 10;
	int number = 0;
	/************************************************************************/
	/* ??clsNumber?????                                                                     */
	/************************************************************************/
	mitk::Point3D point3D;
	mitk::PointSet::Pointer OriginSeeds = mitk::PointSet::New();
	for (int i = 0;i<clsNumber;i++)
	{
		point3D[0] = centerPosition[i][0];
		point3D[1] = centerPosition[i][1];
		point3D[2] = centerPosition[i][2];
		OriginSeeds->GetPointSet()->GetPoints()->InsertElement( i, point3D );
	}
	mitk::DataTreeNode::Pointer pNode = mitk::DataTreeNode::New();
	pNode->SetData( OriginSeeds );
	char name[20];
	sprintf(name,"%d Seeds",clsNumber);
	pNode->SetProperty("name", mitk::StringProperty::New(name));
	pNode->SetProperty("layer", mitk::IntProperty::New(1));
	pNode->SetProperty("color", mitk::ColorProperty::New(1.0,0.0,0.0));
	pNode->SetProperty("pointsize", mitk::FloatProperty::New(2));
	mitk::DataStorage::GetInstance()->Add( pNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();



	while (thresh>theta)
	{
		iter++;
		/*for (int k=0;k<clsNumber;k++)
		{
			std::cout << centerPosition[k][0]<<","<< centerPosition[k][1]<<","<< centerPosition[k][2]<<std::endl;
		}*/

		for (int index = 0;index<vs.size();index++)
		{
			float distance = FLT_MAX,tmp;
			double xcoord = txMap[vs.at(index)],ycoord = tyMap[vs.at(index)],zorord = tzMap[vs.at(index)];
			for (int i = 0;i<clsNumber;i++)
			{
				tmp = sqrt(pow(xcoord-centerPosition[i][0],2)+
					pow(ycoord-centerPosition[i][1],2)+pow(zorord-centerPosition[i][2],2));
				if(tmp<distance)
				{
					distance = tmp;
					cls[index] = i; 
				}
			}
		}
		
		for (int i = 0;i<clsNumber;i++)
		{
			for (int j = 0;j<3;j++)
			{
				curPosition[i][j] = 0;
			}
		}
		//memset(curPosition,0,sizeof(curPosition));
		for (int i=0;i<vs.size();i++)
		{
			curPosition[cls[i]][0]+=txMap[vs.at(i)];
			curPosition[cls[i]][1]+=tyMap[vs.at(i)];
			curPosition[cls[i]][2]+=tzMap[vs.at(i)];
		}

		for (int i = 0;i<clsNumber;i++)
		{
			id = i;
			number = std::count_if(cls,cls+vs.size(),isCls);
			if(number>0)
			{
				curPosition[i][0] /= number;
				curPosition[i][1] /= number;
				curPosition[i][2] /= number;		
			}
			std::cout << i+1<<" Classes: " << number<<std::endl;
		}
		
		thresh = 0;
		for (int i = 0;i<clsNumber;i++)
		{
			thresh += abs(curPosition[i][0]-centerPosition[i][0]);
			thresh += abs(curPosition[i][1]-centerPosition[i][1]);
			thresh += abs(curPosition[i][2]-centerPosition[i][2]);
		}
		
		for (int i = 0;i<clsNumber;i++)
		{
			for (int j = 0;j<3;j++)
			{
				centerPosition[i][j] = curPosition[i][j];
			}
		}
		//memcpy(centerPosition,curPosition,sizeof(centerPosition));
	}

	std::cout << iter << " iterations." << std::endl;
	if(subGraphCount>1)
	{
		for (int i = 2;i<=subGraphCount;i++)
		{
			QDiVertex[i].clear();
		}
	}


	subGraphCount = 1;
	sum = 0;
	for (int i = 0;i<clsNumber;i++)
	{
		id = i;
		subGraphCount++;
		number = std::count_if(cls,cls+vs.size(),isCls);
		std::cout << i+1<<" Classes: " <<number <<std::endl;
		if(number>0)
		{
			for (int j = 0;j<vs.size();j++)
			{
				if(cls[j]==id)
				{
					QDiVertex[subGraphCount].push_back(vs.at(j));
					sum+= pow(txMap[vs.at(j)]-centerPosition[i][0],2)+pow(tyMap[vs.at(j)]-centerPosition[i][1],2)+pow(tzMap[vs.at(j)]-centerPosition[i][2],2);
				}
			}
			
		}
	}
	

	/************************************************************************/
	/*???????????,?????????QDiVertex?,????????                                                                      */
	/************************************************************************/
	std::list<Vertex>::iterator i_beg2;
	std::vector<Edge>::iterator e_iter;
	AdjacencyIterator ad_i,ad_end;
	TVertexIndex vertexIndexMap = get(vertex_index,g);
	for (int i=2;i<=subGraphCount;i++)
	{
		while (!vq.empty())
		{
			vq.pop();
		}
		std::cout << i-1<<" subtree set root number:" << QDiVertex[i].size();
		for (i_beg2 = QDiVertex[i].begin();i_beg2!=QDiVertex[i].end();i_beg2++)
		{
			vq.push(*i_beg2);
		}
		QDiVertex[i].clear();
		while (!vq.empty())
		{
			Vertex v = vq.front();
			vq.pop();
			QDiVertex[i].push_back(v);
			for(tie(ad_i,ad_end) = adjacent_vertices(v,g);ad_i!=ad_end;ad_i++)
			{
				for (e_iter = res.begin();e_iter!=res.end();e_iter++)
				{
					if(GetEdgeTargetIndex(*e_iter)==vertexIndexMap[*ad_i] && GetEdgeSourceIndex(*e_iter)==vertexIndexMap[v])
					{
						vq.push(*ad_i);
						QDiEdge[i].push_back(*e_iter);
					}
				}
			}
		}
		std::cout <<" , " << QDiVertex[i].size()<< std::endl;
		
	}
	/************************************************************************/
	/*                                                                      */
	/************************************************************************/
	

	
	sum = sqrt(sum);

	delete[] seeds;
	delete[] d;
	delete[] cls;
	for (int i = 0;i<clsNumber;i++)
	{
		delete[] centerPosition[i];
		delete[] curPosition[i];
	}
	delete[] centerPosition;
	delete[] curPosition;
	
	return sum;
}

int VesselGraph::GetOutDegree(const Vertex &v) const
{
	int number = 0;
	GraphTraits::adjacency_iterator beg,end;
	for (boost::tie(beg,end) = adjacent_vertices(v,g);beg!=end;beg++)
	{
		number++;
	}
	return number;
}


//set vertex order using strahler number
void VesselGraph::SetVertexOrder()
{
	assert(g!=NULL);
	
	if(vertexOrder.size()!=0)
		vertexOrder.clear();
	
	TVertexIndex vertexIndexMap = get(vertex_index,g);
	VertexIterator v_beg,v_end;
	EdgeIterator e_beg,e_end;
	bool isFinish = false;
	while(!isFinish)
	{
		isFinish = true;
		for (boost::tie(v_beg,v_end) = vertices(g);v_beg!=v_end;v_beg++)
		{
			if(vertexOrder.find(*v_beg)==vertexOrder.end())
			{
				isFinish = false;
				if(GetOutDegree(*v_beg)==0)
					vertexOrder.insert(std::make_pair(vertexIndexMap[*v_beg],1));
				else
				{
					GraphTraits::adjacency_iterator beg,end;
					int minOder = INT_MAX,maxOrder = INT_MIN;
					boost::tie(beg,end) = adjacent_vertices(*v_beg,g);
					for (;beg!=end;beg++)
					{
						if(vertexOrder.find(*beg)==vertexOrder.end())
							break;
						else
						{
							if(vertexOrder[*beg]<minOder)
								minOder = vertexOrder[*beg];
							if(vertexOrder[*beg]>maxOrder)
								maxOrder = vertexOrder[*beg];
						}
					}

					if(beg == end)
					{
						if(minOder==maxOrder)
						{
							vertexOrder.insert(std::make_pair(vertexIndexMap[*v_beg],maxOrder+1));
						}
						else
						{
							vertexOrder.insert(std::make_pair(vertexIndexMap[*v_beg],maxOrder));
						}
					}
				}
			}
		}

	}

	/*
	for (boost::tie(v_beg,v_end) = vertices(g);v_beg!=v_end;v_beg++)
	{
		std::cout <<"Index: " <<vertexIndexMap[*v_beg]<<" , "<<vertexOrder[*v_beg]<<std::endl;
	}
	*/
	
}