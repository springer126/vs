#ifndef _QmitkAnnotation_H
#define _QmitkAnnotation_H

#include "render.h"
#include <qobject.h>
#include <qlistview.h>
#include <qcheckbox.h>

#include <vtkRenderer.h>
#include <vtkCellPicker.h>
#include "mitkDataTree.h"
#include "QmitkTreeNodeSelector.h"
#include "DataTreeViewItem.h"


//mitk headers
#include "mitkColorSequenceRainbow.h"
#include "mitkCommon.h"
#include "mitkDataTreeNode.h"
#include "mitkDataTree.h"
#include <mitkSurface.h>
#include <mitkStandaloneDataStorage.h>
#include <itkImage.h>
#include <vtkPolyData.h>
#include <vtkPlane.h>
#include <vtkImplicitPlaneWidget.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAxesActor.h>
#include <vtkImageExport.h>
#include <vtkImageImport.h> 
#include <itkVTKImageImport.h>
#include <itkVTKImageExport.h>
#include "VesselGraph.h"
#include <map>
#include <vector>

const int segExp = 8;

class QmitkAnnotation : public Render
{   Q_OBJECT

public:
	static QmitkAnnotation* m_pInstance;

	vtkRenderer*    m_pRenderer;
	vtkRenderWindow* m_pRenWin; 
	//vtkRenderWindowInteractor *m_interactor;
	
	
	static QmitkAnnotation* Instance(QObject* parent, mitk::DataTreeIteratorBase* iterator = NULL, mitk::DataTreeNode* SelectedNode = NULL);
	static void Destroy();
   //Constructor
   QmitkAnnotation(QObject* parent);

   //Destructor
   ~QmitkAnnotation();


   /*!     \brief method for creating the connections of main and control widget     */     
   virtual void CreateConnections();

private:
	mitk::DataTreeNode* SelectedNode;
	bool isInterpolate;
	bool isInitDataTree;
	template <typename TPixel, unsigned int VImageDimension>
	void CastPixelType(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);
	/**2012-09-11**/
	template <typename TPixel, unsigned int VImageDimension>
	void MyCastPixelType(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);

	template <typename TPixel, unsigned int VImageDimension>
	void GetMaxMinPixelIntensity(itk::Image<TPixel, VImageDimension> *itkImage, int &maxIntensity, int &minIntensity);

	template <typename TPixel, unsigned int VImageDimension>
	void Interpolate(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);

	template <typename TPixel, unsigned int VImageDimension>
	void CurvatureFlow(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);

	template <typename TPixel, unsigned int VImageDimension>
	void CopyImage(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);

	template <typename TPixel, unsigned int VImageDimension>
	void SeparateLabelImage(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);

	template <typename ImageType>
	void ConnectVTKToITK(vtkImageExport* in, itk::VTKImageImport<ImageType>* out);

protected:
	mitk::Image* m_MitkImage;
	mitk::Surface* m_MitkSurface;
	vtkPolyData *polydata;
	vtkPlane *plane;
	vtkImplicitPlaneWidget *pWidget;
	vtkActor *actor;
	vtkAxesActor *axes;
	VesselGraph* vesselGraph;
	long voxelNumberOfImage;
	vtkCellPicker *picker;
	vtkUnstructuredGrid *grid;
	std::map<int,vtkActor*> segMap;
	std::map<const QListViewItem*,vtkActor*> actorMap;
	//vtkActor *axes;
	double opacityValue;
	char* clipPolyData;
	char* clipImageData;
	
	int numOfActors;

	//for cacl  information of segment
	mitk::DataTreeIteratorBase *m_DataTreeIteratorBase;
	mitk::DataTreeIteratorClone selectedImage;
	mitk::DataTreeNode::ConstPointer m_Node;
	
	mitk::ColorSequenceRainbow m_RainbowColor;
	mitk::Color color;

	virtual void closeEvent(QCloseEvent*);

	void ClipPolydata(double (&normal)[3],double (&origin)[3],vtkPolyData *source,vtkPolyData *one,vtkPolyData *another);
	void ClipImagedata(double (&normal)[3],double (&origin)[3], mitk::Image *image,mitk::DataTreeNode *parent,char *name1,char *name2);
	void InitDataTree();
	void setNode(vtkPolyData *polydata,mitk::DataTreeNode *parent,char *nodeName);
	void addImageNode(mitk::Image *imagedata,mitk::DataTreeNode *parent,char *nodeName);
	
	template <typename TPixel, unsigned int VImageDimension>
	long VoxelNumberOfImage(itk::Image<TPixel, VImageDimension> *itkImage, long &param1, long &param2);

	void activeActorProperty(bool);

	void DisplayVesselTree(std::vector<Vertex> vs,std::vector<Edge> es);
	void DisplaySubtreeNodes();
protected slots:  

   void RenderButtonClicked();

   void SetActorColor();

   void HighlightButtonClicked();
   void setOpacity();
   void EnableInterpolation();
   void EnableAxes();
   void SegmentByMiddleHepaticVein();
   void SegmentByLeftHepaticVein();
   void SegmentByRightHepaticVein();
   void SegmentByPortalVein();
   void SegmentISeg();
   void Clip();
   void SetDataNodeChecked(DataTreeViewItem*, bool);

   void CaclwithAnnotation();
   
   void DoDisplayVesselTree();
   void DoVesselDivision();
   void AutomaticLiverSegment();
   void ManualLiverSegment();
   void CreateNewSubTree();
   void TreeChanged();
   void itemSelected(QListViewItem*);
   void actorItemSelected(QListViewItem  *);
   void ImageSelected(mitk::DataTreeIteratorClone imageIt);

   void EnableVisible();

   void addActorNode(vtkActor *actor,QString str,mitk::Color color);

   
};
#endif 
