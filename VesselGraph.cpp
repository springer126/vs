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

using namespace boost;
typedef boost::property_map<Graph, vertex_index_t>::type TVertexIndex;
typedef boost::property_map<Graph, edge_index_t>::type TEdgeIndex;

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

VesselGraph::VesselGraph(const char * Path):g(0),subGraphCount(0),separateArgument(0.6)
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

	mitk::PointSet::Pointer points = mitk::PointSet::New();
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
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


bool VesselGraph::GetSubGraph()
{
	//Reset();
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
	
	if(QDiVertex[subGraphCount].size()<=1)
	{
		QDiVertex[subGraphCount].clear();
		--subGraphCount;
		std::cout << "[Deprecated:1 point.]" <<std::endl;
		return false;
	}

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
	AccessByItk_1(image,NNSA,&resultImage);
	return resultImage;
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
	QDiVertex.clear();
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
		//Ã¿¶Î8¶ÎÌåËØÖµ£º20-90
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

unsigned int VesselGraph::GetVertexIndex(Vertex v) const 
{
	TVertexIndex vertexIndexMap = get(vertex_index,g);
	return vertexIndexMap[v];
}

unsigned int VesselGraph::GetEdgeSourceIndex(Edge e) const 
{
	 return GetVertexIndex(source(e,g));
}

unsigned int VesselGraph::GetEdgeTargetIndex(Edge e) const
{
	return GetVertexIndex(target(e,g));
}

void VesselGraph::RadiusFilter(double threshold)
{
	assert((vs.size()>0));
	assert((es.size()>0));
	TVertexIndex vertexIndexMap = get(vertex_index,g);
	std::vector<Edge>::iterator i_beg;
	std::vector<Vertex>::iterator i_beg2;
	for (i_beg = es.begin();i_beg!=es.end();i_beg++)
	{
		if(edgeMrMap[*i_beg]>threshold)
		{
			for (i_beg2 = vs.begin();i_beg2!=vs.end();i_beg2++)
			{
				if(vertexIndexMap[*i_beg2] ==source(*i_beg,g))
				{
					vs.erase(i_beg2);
					break;
				}
			}
			es.erase(i_beg);
		}
	}

	std::cout <<"Vertex Number :" << vs.size()<<" , Edge Number:" << es.size() << std::endl;
	for (i_beg = es.begin();i_beg!=es.end();i_beg++)
	{
		std::cout << "("<<source(*i_beg,g)<<","<<target(*i_beg,g)<<")"<<std::endl;
	}
}