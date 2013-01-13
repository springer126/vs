#include "QmitkAnnotation.h"

#include <qfiledialog.h>
#include <qaction.h>
#include <qmessagebox.h>
#include <qlayout.h>

#include <iostream>
#include <string>
#include <expat.h>
#include <sstream>
#include <fstream>
#include <utility>
#include <Windows.h>
#include <tchar.h>

#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkCutter.h>
#include <vtkDecimatePro.h>
#include <vtkImageResample.h>
#include <vtkImageGaussianSmooth.h> 
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>
#include <vtkSmartPointer.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkTextProperty.h>
#include <vtkCaptionActor2D.h>
#include <vtkPropAssembly.h>
#include <vtkImplicitPlaneWidget.h>
#include "vtkTIPWCallback.h"
#include <vtkImageExport.h>
#include <vtkExtractEdges.h>
#include <vtkTubeFilter.h>
#include <vtkImageCast.h>
#include <vtkImageChangeInformation.h>
#include <vtkDICOMImageReader.h>
#include <vtkSphereSource.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkLineSource.h>
#include <vtkCamera.h>
#include <vtkPointSet.h>
#include <vtkPoints.h> 
#include <vtkCellArray.h> 
#include <vtkVertex.h>
#include <vtkLine.h> 
#include <vtkFloatArray.h>
#include <vtkSphereSource.h>
#include <vtkThresholdPoints.h>
#include <vtkGlyph3D.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkSTLReader.h>
#include <vtkContourFilter.h>
#include <vtkTextMapper.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkTextActor.h>
#include <vtkExtractEdges.h> 
#include <vtkTubeFilter.h> 

#include <qslider.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qlistview.h> 
#include <qpixmap.h>
#include <qpainter.h>
#include <qfontmetrics.h>
#include <qevent.h>
#include <qapplication.h>
#include <qheader.h>
#include <qtoolbutton.h>
#include <qpalette.h> 
#include <qcolor.h> 
#include <qcolordialog.h> 
#include <qbrush.h> 

#include <algorithm>
#include <sstream>



#include "itkImage.h"
#include <itkCastImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkTranslationTransform.h>
#include "itkImageDuplicator.h"
#include "itkVTKImageToImageFilter.h"
#include <itkVTKImageImport.h>
#include <itkCurvatureFlowImageFilter.h>

#include "DataTreeViewItem.h"
#include "mitkImageAccessByItk.h"
#include "mitkDataStorage.h"
#include "mitkStringProperty.h"
#include "mitkDataTreeNode.h"
#include "mitkProperties.h"
#include <mitkPlaneCutFilter.h>
#include <mitkPlaneGeometry.h>
#include <mitkVector.h>
#include <mitkImage.h>
#include <mitkPlaneCutFilter.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include <mitkDataTreeFilterFunctions.h>
#include <mitkVolumeCalculator.h>
#include <mitkLabeledImageToSurfaceFilter.h>
#include <mitkRenderingManager.h>
#include <mitkSurface.h>
#include "mitkImageAccessByItk.h"
#include "mitkITKImageImport.h"
#include <mitkColorProperty.h>
#include <mitkSurface.h>
#include <mitkManualSegmentationToSurfaceFilter.h>

#include <QmitkDataTreeComboBox.h>
#include "QmitkToolReferenceDataSelectionBox.h"
#include "QmitkToolWorkingDataSelectionBox.h"
#include <QmitkRenderWindow.h>

#define ROUND(a)     ((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))

namespace mitk
{
	class IsSegmentNode: public mitk::DataTreeFilterFunction
	{
	public:
		virtual ~IsSegmentNode() {}

		virtual bool NodeMatches(DataTreeNode* node) const
		{
			//bool result = node && node->GetData() && dynamic_cast<mitk::Image*>(node->GetData());
			bool result = node && node->GetData();
			if(result)
			{
				std::string name;
				node->GetName(name);
				//result = ! (name=="volume threshold overlay image");
				std::cout <<name << " : " <<name.find("Segment")<<std::endl;
				result =!(name.find("Segment")==name.npos);
			}
			return result;
		} 

		virtual DataTreeFilterFunction* Clone() const
		{
			return new IsSegmentNode();
		}
	};

	class vtkPickCallback : public vtkCommand
	{
	public:
		static vtkPickCallback *New() 
		{
			return new vtkPickCallback; 
		}

		virtual void Execute(vtkObject *caller, unsigned long, void*)
		{
			vtkCellPicker *picker = reinterpret_cast<vtkCellPicker *>(caller);
			if (picker->GetCellId() == -1 )
			{
				std::cout << "No cell selected." << std::endl;
			}
			else
			{
				double pickPos[3];
				//picker->GetPickPosition( pickPos );
				picker->GetPCoords(pickPos);
				double xp = pickPos[0];
				double yp = pickPos[1];
				double zp = pickPos[2];
				cellid = picker->GetCellId();
				std::cout << "Picked a cell cellid = "<<cellid<<std::endl;

			}
		}
		vtkPickCallback()
		{
			cellid = 0;
		}
		vtkIdType cellid;
	};
}
//=====================================================//


QmitkAnnotation* QmitkAnnotation::m_pInstance = 0;

QmitkAnnotation* QmitkAnnotation::Instance(QObject* parent, mitk::DataTreeIteratorBase* iterator,mitk::DataTreeNode* SelectedNode )
{
	DWORD dwScreenX,dwScreenY; 
	dwScreenX = ::GetSystemMetrics(SM_CXSCREEN); 
	dwScreenY = ::GetSystemMetrics(SM_CYSCREEN);

	std::cout << "dwScreenX = " << dwScreenX << std::endl;
	std::cout << "dwScreenY = " << dwScreenY << std::endl;

	DWORD dwSize = dwScreenX > dwScreenY ? dwScreenY : dwScreenX;
	DWORD renderWindowSize = (DWORD)(0.8*dwSize);

	DWORD innerDwScreenX,innerDwScreenY;
	innerDwScreenX = (DWORD)(0.75*dwScreenX);
	innerDwScreenY = (DWORD)(0.75*dwScreenY);
	
	if (!m_pInstance) 
	{
		m_pInstance = new QmitkAnnotation(parent);
		printf("iterator = %p\n", iterator);
		if(iterator)
		{
	//		m_pInstance->m_TreeNodeSelector->SetDataTreeNodeIterator(iterator);
			std::cout << "SetDataTreeNodeIterator" << std::endl;
 			while (m_pInstance->listView->firstChild())
 			{
 				delete m_pInstance->listView->firstChild();
 			}
 			mitk::DataTreeIteratorClone tempIt = iterator;
			m_pInstance->m_DataTreeIteratorBase  = iterator;
 			if (!tempIt->IsAtEnd())
 			{
				DataTreeViewItem* dataView = new DataTreeViewItem(m_pInstance->listView, "Loaded Data", "root", tempIt.GetPointer(), SelectedNode);	
 				++tempIt;
 			}

// 			int count = m_pInstance->listView->childCount();
// 			std::cout << "count = " <<count << std::endl;
// 			//m_pInstance->listView->takeItem(m_pInstance->listView->lastItem());
// 			QListViewItem *firstChild = m_pInstance->listView->firstChild();
// 		   QListViewItem *sibling = firstChild;
// 		   QListViewItem *temp;
// 		   while (sibling)
// 		   { 
// 			   sibling = sibling->nextSibling();
// 		   }
// 
// 			count = m_pInstance->listView->childCount();
// 			std::cout << "count = " <<count << std::endl;

          QListViewItem *firstChild = m_pInstance->listView->firstChild();
		  if(firstChild)
		  {
			//  QCheckListItem *p = new QCheckListItem (firstChild, QString("QCheckListItem"),QCheckListItem::CheckBox ); 
			//  p->setOn(true);
			//  firstChild->insertItem(p);

			  QListViewItem *item = m_pInstance->listView->findItem(QString("Widgets"), 0 );
			  if(item)
			  {
				  QString text = item->text(0);
				//  std::cout << text.ascii() << std::endl;
				  firstChild->takeItem(item);
			  }

		  }

		}
		
		m_pInstance->show();
		m_pInstance->SelectedNode = SelectedNode;
		
		m_pInstance->segListView->clear();
		/*m_pInstance->actorListBox->setColumnText(0,"Actor"); 
		while(m_pInstance->actorListBox->columns()>1)
			m_pInstance->actorListBox->removeColumn(1);
		m_pInstance->actorListBox->clear();*/
	}

	
	
	
	//m_pInstance->m_ToolReferenceDataSelectionBox
	return m_pInstance;
}
void QmitkAnnotation::Destroy()
{
	if (m_pInstance) delete m_pInstance;
}

void QmitkAnnotation::closeEvent(QCloseEvent* e)
{
    Destroy();
    e->accept();
}

mitk::PixelType originalPixelType;
// constructor
QmitkAnnotation::QmitkAnnotation(QObject* parent)
{
   CreateConnections();
   m_pRenderer = vtkRenderer::New();
   m_pRenderer->SetBackground(0.20, 0.27 ,0.32);
   m_pRenWin = RenderWindow->GetRenderer()->GetRenderWindow();
   m_pRenWin->AddRenderer(m_pRenderer); 
   picker  = vtkCellPicker::New();
   clipPolyData = clipImageData = NULL;

   m_MitkImage=NULL;
   m_MitkSurface=NULL;
   actor = NULL;
   polydata  = NULL;
   plane = NULL;
   pWidget = NULL;
   axes = NULL;
   vesselGraph = NULL;
   voxelNumberOfImage = 0;
   QString s;
   s.sprintf("%.3f",0.60);
   radiusRatio->setText(s);
   isInterpolate = false;
   isInitDataTree = false;
   Isaxes->setChecked(false);
   
   btnMiddleHepaticVein->setEnabled(false);
   btnLeftHepaticVein->setEnabled(false);
   btnRightHepaticVein->setEnabled(false);
   btnPortalVein->setEnabled(false);
   btnClip->setEnabled(false);
   

   numOfActors = 0;

   m_btnNewSubTree->setEnabled(false);
   m_btnDivideTree->setEnabled(false);
   m_btnManuLiverSegment->setEnabled(false);
   segListView->QListView::addColumn(QString("Segmentation"));
   segListView->setColumnText (0,QString(""));
   //segListView->setColumnText (1,QString("Segmentation"));

   activeActorProperty(false);
}

//Connections
void QmitkAnnotation::CreateConnections()
{
    QObject::connect( (QObject*)(RenderButton), SIGNAL(clicked()),(QObject*) this, SLOT(RenderButtonClicked()));
	QObject::connect( (QObject*)(opacity), SIGNAL(valueChanged(int)),(QObject*) this, SLOT(setOpacity()));
	QObject::connect((QObject*)(Interpolation), SIGNAL(toggled(bool)), (QObject*)this, SLOT(EnableInterpolation()) );
	QObject::connect((QObject*)(Isaxes), SIGNAL(toggled(bool)), (QObject*)this, SLOT(EnableAxes()) );
	QObject::connect((QObject*)(btnMiddleHepaticVein), SIGNAL(clicked()), (QObject*)this, SLOT(SegmentByMiddleHepaticVein()) );
	QObject::connect((QObject*)(btnLeftHepaticVein), SIGNAL(clicked()), (QObject*)this, SLOT(SegmentByLeftHepaticVein()) );
	QObject::connect((QObject*)(btnRightHepaticVein), SIGNAL(clicked()), (QObject*)this, SLOT(SegmentByRightHepaticVein()) );
	QObject::connect((QObject*)(btnPortalVein), SIGNAL(clicked()), (QObject*)this, SLOT(SegmentByPortalVein()) );
	QObject::connect((QObject*)(btnISeg), SIGNAL(clicked()), (QObject*)this, SLOT(SegmentISeg()) );
	QObject::connect((QObject*)(btnClip), SIGNAL(clicked()), (QObject*)this, SLOT(Clip()) );
	QObject::connect((QObject*)(btnCacl), SIGNAL(clicked()), (QObject*)this, SLOT(CaclwithAnnotation()) );
	
	QObject::connect((QObject*)(m_btnActorColor), SIGNAL(clicked()), (QObject*)this, SLOT(SetActorColor()) );
	
	QObject::connect( (QObject*)(listView), SIGNAL( clicked( QListViewItem* ) ), this, SLOT( actorItemSelected( QListViewItem* ) ) );

	QObject::connect( (QObject*)(m_TreeNodeSelector), SIGNAL(Activated(mitk::DataTreeIteratorClone)),(QObject*) this, SLOT(ImageSelected(mitk::DataTreeIteratorClone)) );
	QObject::connect( (QObject*)(segListView), SIGNAL( clicked( QListViewItem* ) ), this, SLOT( itemSelected( QListViewItem* ) ) );
	QObject::connect((QObject*)( m_btnDisplayVesselTree ), SIGNAL( clicked() ), (QObject*)this, SLOT(DoDisplayVesselTree()) );
	QObject::connect((QObject*)( m_btnNewSubTree ), SIGNAL( clicked() ), (QObject*)this, SLOT(CreateNewSubTree()) );
	QObject::connect((QObject*)( m_btnDivideTree ), SIGNAL( clicked() ), (QObject*)this, SLOT(DoVesselDivision()) );
	QObject::connect((QObject*)( m_btnAutoLiverSegment ), SIGNAL( clicked() ), (QObject*)this, SLOT(AutomaticLiverSegment()) );
	QObject::connect((QObject*)( m_btnManuLiverSegment ), SIGNAL( clicked() ), (QObject*)this, SLOT(ManualLiverSegment()) );
	//connect((QObject*)NULL, SIGNAL(statusChange(DataTreeViewItem*, bool)), this, SLOT(SetDataNodeChecked(DataTreeViewItem*, bool)) );
	//QObject::connect( (QObject*)())
	QObject::connect((QObject*)( actorVisible ), SIGNAL( toggled(bool) ), (QObject*)this, SLOT(EnableVisible()) );
}


void QmitkAnnotation::RenderButtonClicked() 
{
	std::cout << "RenderButtonClicked! " << std::endl;
	m_pRenderer->SetBackground(0.0, 0.0, 0.0);
	int *origin = m_pRenderer->GetOrigin();
	double originOfRender[3];
	for (int index=0;index!=3;index++)
	{
		std::cout << "2 : Origin of render is"<<origin[index]<<" "<<std::endl;
		originOfRender[index] = (double)origin[index];
	}
	
	if (axes == NULL)
	{
		vtkSmartPointer<vtkAxesActor> axes2 = vtkSmartPointer<vtkAxesActor>::New();
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		transform->Identity();
		axes2->SetShaftTypeToCylinder();
		axes2->SetUserTransform( transform );
		axes2->SetXAxisLabelText( "x" );
		axes2->SetYAxisLabelText( "y" );
		axes2->SetZAxisLabelText( "z" );
		axes2->SetTotalLength( 210.5, 210.5, 210.5 );
		axes2->SetCylinderRadius( 0.1 * axes2->GetCylinderRadius() );
		axes2->SetConeRadius    ( 0.1 * axes2->GetConeRadius() );
		axes2->SetSphereRadius  ( 0.1 * axes2->GetSphereRadius() );
		vtkTextProperty* tprop = axes2->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
		tprop->ItalicOn();
		tprop->ShadowOn();
		tprop->SetFontFamilyToTimes();
		tprop->SetColor(0.0,1.0,1.0);
		tprop->SetFontSize(32);
		tprop->SetBold(false);
		axes2->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy( tprop );
		axes2->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy( tprop );		
		axes = axes2;
		axes->SetVisibility(false);
		m_pRenderer->AddViewProp( axes2 );

	}
	
	double standardDiviation;
	QListViewItemIterator it(m_pInstance->listView);
	it++;
	while (it.current()) 
	{
		DataTreeViewItem *item = (DataTreeViewItem *)it.current();
		std::cout<<item->text(0)<<std::endl;
		mitk::Image::Pointer copyImage;
		mitk::Image::Pointer castImage;
		mitk::Image::Pointer interpolateImage,curvatureFlowImage;
		if (item->isOn())
		{
			std::cout<<"Selected and begin render!"<<std::endl;
			standardDiviation = 0.35;
			m_MitkImage = dynamic_cast<mitk::Image*> (item->GetDataTreeNode()->GetData());
			//m_MitkImage = dynamic_cast<mitk::Image*> (mitk::DataStorage::GetInstance()->GetNamedNode("segmented_liver")->GetData());
			if (m_MitkImage==NULL)
			{
				std::cout<<"Surprise: is not image data."<<std::endl;
				m_MitkSurface = dynamic_cast<mitk::Surface*> (item->GetDataTreeNode()->GetData());
				if (m_MitkSurface!=NULL)
				{
					std::cout<<"Surprise: is surface data."<<std::endl;
					vtkTriangleFilter *triangleFilter = vtkTriangleFilter::New();
					triangleFilter->SetInput(m_MitkSurface->GetVtkPolyData());
					vtkDecimatePro *decimatePro = vtkDecimatePro::New();
					decimatePro->SetInputConnection(triangleFilter->GetOutputPort());
					decimatePro->SetTargetReduction(0.5);
					decimatePro->PreserveTopologyOn();
					vtkSmoothPolyDataFilter *smoother = vtkSmoothPolyDataFilter::New();
					smoother->SetInputConnection(decimatePro->GetOutputPort());
					smoother->SetNumberOfIterations( 50 );
					smoother->FeatureEdgeSmoothingOn();

					vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
					mapper->SetInput(smoother->GetOutput());
					vtkActor *iactor = vtkActor::New();
					iactor->SetMapper(mapper);
					mitk::Color color = m_RainbowColor.GetNextColor();
					iactor->GetProperty()->SetColor(color.GetRed(),color.GetGreen(),color.GetBlue());
					m_pRenderer->AddActor(iactor);
					m_pRenderer->ResetCamera();
					mitk::RenderingManager::GetInstance()->AddRenderWindow(m_pRenWin);
					m_pRenWin->Render();
					return;
				}
			}

			if (m_MitkImage==NULL && m_MitkSurface==NULL)
			{
				std::cout<<"Surprise: what a image."<<std::endl;
				return;
			}
			
			std::cout << "Item :image pixel type is :"<<m_MitkImage->GetPixelType().GetItkTypeAsString() << std::endl;
			int maxIntensity, minIntensity;
			AccessByItk_2(m_MitkImage, GetMaxMinPixelIntensity, maxIntensity, minIntensity);
			std::cout << "maxIntensity: " << maxIntensity << std::endl;
			std::cout << "minIntensity: " << minIntensity << std::endl;
			if (isInterpolate)
			{
				standardDiviation = 3.0;
			}
			if(actor == NULL)
			{
				AccessByItk_1(m_MitkImage, CopyImage, &copyImage);
				AccessByItk_1(m_MitkImage, MyCastPixelType, &castImage);
				//AccessByItk_1(castImage, CurvatureFlow, &curvatureFlowImage);
			}
			else
			{
				castImage = m_MitkImage;	
			}

			mitk::ManualSegmentationToSurfaceFilter::Pointer filter = mitk::ManualSegmentationToSurfaceFilter::New();
			if (filter.IsNull())
			{
				std::cout<<"NULL Pointer for ManualSegmentationToSurfaceFilter"<<std::endl;
				return;
			}
			filter->SetInput( castImage );  
			filter->UseGaussianImageSmoothOn();
			filter->SetGaussianStandardDeviation( standardDiviation );
			filter->SetUseGaussianImageSmooth( true );
			filter->SetThreshold( 1.0 ); 
			filter->SetTargetReduction( 0.5 );
			

			
			//vtkImageChangeInformation *indexCoordinatesImageFilter = vtkImageChangeInformation::New();
			//indexCoordinatesImageFilter->SetInput(curvatureFlowImage->GetVtkImageData());
			//indexCoordinatesImageFilter->SetOutputOrigin(0.0,0.0,0.0);
			//vtkImageGaussianSmooth *gaussian = vtkImageGaussianSmooth::New();
			//gaussian->SetInputConnection(indexCoordinatesImageFilter->GetOutputPort());
			//gaussian->SetStandardDeviation(standardDiviation);
			//indexCoordinatesImageFilter->Delete();
			/*if(curvatureFlowImage->IsEmpty())
			{
				std::cout << "curvatureFlowImage==NULL" << std::endl;
				return;	
			}
			vtkMarchingCubes *skinExtractor = vtkMarchingCubes::New();
			skinExtractor->ComputeScalarsOff();
			skinExtractor->SetInput(curvatureFlowImage->GetVtkImageData());*/
			//skinExtractor->SetInputConnection(gaussian->GetOutputPort());
			//gaussian->Delete();
			//skinExtractor->SetValue(0, 1);
			
			//----------------------三角面片约简------------------------------
			//vtkTriangleFilter *triangle= vtkTriangleFilter::New();
			//triangle->SetInputConnection(skinExtractor->GetOutputPort());
			//vtkDecimatePro *decimate = vtkDecimatePro::New();
			//decimate->SetInputConnection(triangle->GetOutputPort());//RC++
			//decimate->SetTargetReduction(0.5);
			//decimate->PreserveTopologyOn();
			//polydata->Delete();//RC--
			//polydata = decimate->GetOutput();
			//polydata->Register(NULL);//RC++
			//decimate->Delete();
			//----------------------------------------------------------------
			


			vtkSmoothPolyDataFilter *smoother = vtkSmoothPolyDataFilter::New();
			smoother->SetInput(filter->GetOutput()->GetVtkPolyData());
			//smoother->SetInput(skinExtractor->GetOutput());
			smoother->SetNumberOfIterations( 50 );
			smoother->FeatureEdgeSmoothingOn();

			vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
			vtkActor* iactor = vtkActor::New();
			/*if(actor == NULL)
			{
				this->actor = iactor;
				
			}*/
			polydata = smoother->GetOutput();
			mapper->SetInput(polydata);
			iactor->SetMapper(mapper);
			iactor->SetOrigin(originOfRender);
			iactor->GetProperty()->SetOpacity(0.6);
			mitk::Color color = m_RainbowColor.GetNextColor();
			iactor->GetProperty()->SetColor(color.GetRed(),color.GetGreen(),color.GetBlue());
			m_pRenderer->AddActor(iactor);
			
			QString str;
			str.sprintf("Surface %d",++numOfActors);
		    QCheckListItem *actorItem = new QCheckListItem(item,str);
			QPixmap pixmap(12,8);
			std::cout << color.GetRed()<< " "<<color.GetGreen()<<" "<<color.GetBlue()<<std::endl;
			pixmap.fill(QColor(color.GetRed()*255,color.GetGreen()*255,color.GetBlue()*255));
			actorItem->setPixmap(0,pixmap);
			actorItem->setEnabled(true);
			actorMap.insert(std::make_pair(actorItem,iactor));
			std::cout<<"render over!"<<std::endl;
		
		}

	 it++;

	}

	//=================DataStorage添加肝脏PolyData节点=====================//
	if(mitk::DataStorage::GetInstance()->GetNamedNode("*Liver Polydata*")==NULL)
	{
		std::cout << "Liver Polydata DataNode is existed." << std::endl;
		mitk::DataTreeNode::Pointer liverNode = mitk::DataTreeNode::New();
		mitk::Surface::Pointer surface = mitk::Surface::New();
		surface->SetVtkPolyData(polydata);
		liverNode->SetData(surface);
		liverNode->SetProperty("name", mitk::StringProperty::New("*Liver Polydata*"));
		liverNode->SetProperty("opacity", mitk::FloatProperty::New(0.0));
		//rightLobe->SetColor( m_RainbowColor.GetNextColor());
		liverNode->SetVisibility(true);

		mitk::DataStorage::GetInstance()->Add( liverNode );
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	}
	//====================================================================//
	

	

	m_pRenderer->ResetCamera();
	//mitk::RenderingManager::GetInstance()->AddRenderWindow(m_pRenWin);
	m_pRenWin->Render();
}

//this function outdate
void QmitkAnnotation::HighlightButtonClicked() 
{
	std::cout << "Skeletonization Button clicked!" << std::endl;
	QString sPath = QFileDialog::getOpenFileName(
                   "/",
                   "All Files (*)",
                   this,
                   "open file dialog",
                   "Choose a file to open" );
}

QmitkAnnotation::~QmitkAnnotation()
{
   m_pInstance = NULL;
} 



//this function outdate just for testing
void QmitkAnnotation::SetDataNodeChecked(DataTreeViewItem* node, bool inflag)
{
	std::cout << "SetDataNodeChecked " << std::endl;
	std::cout << node->text().ascii() << std::endl;
	std::cout << inflag << std::endl;
}


void QmitkAnnotation::setOpacity()
{
	if(actor==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "No actor to set opacity!");
		return;
	}

	std::cout<<opacity->value()/100.0<<std::endl;
	this->actor->GetProperty()->SetOpacity(opacity->value()/100.0);
	m_pRenWin->Render();
	
}

void QmitkAnnotation::EnableInterpolation()
{
	isInterpolate = isInterpolate?false:true;
	std::cout<<"isInterpolate is "<<isInterpolate<<std::endl;
}
 
template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::CastPixelType(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	/**duplicate a image to clip 8 segment**/
	typedef itk::Image<TPixel, VImageDimension> TImageType;
	typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage( itkImage );
	duplicator->Update();
	typename TImageType::Pointer clonedImage = duplicator->GetOutput();
	mitk::Image::Pointer resultImage = mitk::ImportItkImage( clonedImage );
	mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
	newNode->SetData(resultImage);
	newNode->SetProperty("name", mitk::StringProperty::New("Liver Image"));
	newNode->SetProperty("opacity", mitk::FloatProperty::New(0.0));
	mitk::DataStorage::GetInstance()->Add( newNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	
	typedef itk::Image<TPixel, VImageDimension> InputImageType;
	originalPixelType = m_MitkImage->GetPixelType();
	//originalPixelType = pointer->GetPixelType();
	if(originalPixelType == typeid(unsigned char))
	{
		//std::cout << "CastPixelType: unsigned char!" << std::endl;
		typedef itk::Image<unsigned char,  VImageDimension> OutputImageType;
		typedef itk::CastImageFilter<InputImageType, OutputImageType> CasterType;

		CasterType::Pointer caster = CasterType::New();
		caster->SetInput(itkImage);
		try
		{
			caster->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "CastPixelType exception thrown." 
				<< std::endl;
			std::cerr << excp << std::endl;
			return ;
		}

		OutputImageType::Pointer output = caster->GetOutput();
		*pointer = mitk::ImportItkImage(output);
	}
	else if(originalPixelType == typeid(short))
	{
		//std::cout << "CastPixelType: short!" << std::endl;
		typedef itk::Image<short,  VImageDimension> OutputImageType;
		typedef itk::CastImageFilter<InputImageType, OutputImageType> CasterType;

		CasterType::Pointer caster = CasterType::New();
		caster->SetInput(itkImage);
		try
		{
			caster->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "CastPixelType exception thrown." 
				<< std::endl;
			std::cerr << excp << std::endl;
			return ;
		}

		OutputImageType::Pointer output = caster->GetOutput();
		*pointer = mitk::ImportItkImage(output);
	} 
	else 
	{   
		char buf[128];
		sprintf(buf, "%s not supported by this method yet!", originalPixelType.GetItkTypeAsString().c_str());
		QMessageBox::information( NULL, "Virtual Surgery functionality", buf);
	}

}


template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::MyCastPixelType(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	std::cout<< "Function MyCastPixelType :Begin."<<std::endl;

	///***************duplicate a image to clip 8 segment******************/
	//typedef itk::Image<TPixel, VImageDimension> TImageType;
	//typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	//DuplicatorType::Pointer duplicator = DuplicatorType::New();
	//duplicator->SetInputImage( itkImage );
	//duplicator->Update();
	//typename TImageType::Pointer clonedImage = duplicator->GetOutput();
	//mitk::Image::Pointer resultImage = mitk::ImportItkImage( clonedImage );
	//mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
	//newNode->SetData(resultImage);
	//newNode->SetProperty("name", mitk::StringProperty::New("Liver Image"));
	//newNode->SetProperty("opacity", mitk::FloatProperty::New(0.0));
	//mitk::DataStorage::GetInstance()->Add( newNode );
	//mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	///*********************************************************************/

	originalPixelType = m_MitkImage->GetPixelType();
	if (originalPixelType == typeid(float))
	{
		std::cout<< "Function MyCastPixelType : Origin image pixel type is float."<<std::endl;
		*pointer = mitk::ImportItkImage(itkImage);
	}
	else
	{
		typedef itk::Image<TPixel, VImageDimension> InputImageType;
		typedef float OutputImagePixelType;
		typedef itk::Image<OutputImagePixelType,  VImageDimension> OutputImageType;
		typedef itk::CastImageFilter<InputImageType, OutputImageType> CasterType;

		CasterType::Pointer caster = CasterType::New();
		caster->SetInput(itkImage);
		try
		{
			caster->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "CastPixelType exception thrown." 
				<< std::endl;
			std::cerr << excp << std::endl;
			return ;
		}

		OutputImageType::Pointer output = caster->GetOutput();
		*pointer = mitk::ImportItkImage(output);
	}
	
	std::cout<< "Function MyCastPixelType End."<<std::endl;
	
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::Interpolate(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	originalPixelType = m_MitkImage->GetPixelType();
	if (originalPixelType == typeid(float))
	{
		//originalPixelType = short;
	}
	
	typedef itk::Image< TPixel, VImageDimension >  InputImageType;
	InputImageType::RegionType region = itkImage->GetLargestPossibleRegion();
	InputImageType::SpacingType space = itkImage->GetSpacing();
	InputImageType::DirectionType direction = itkImage->GetDirection();
	InputImageType::PointType origin = itkImage->GetOrigin();
	InputImageType::SizeType size = region.GetSize();
	
	std::cout << "region: " << region << std::endl;
	std::cout << "space:" << space << std::endl;

	space[2] = space[2]/4.0;
	size[2] <<=2;

	std::cout << "direction:" << direction << std::endl; 

	if(originalPixelType == typeid(unsigned char))
	{
		typedef itk::Image< unsigned char, VImageDimension > OuputImageType;
		typedef itk::ResampleImageFilter< InputImageType,OuputImageType,double>  FilterType;
		FilterType::Pointer resampler = FilterType::New();
		resampler->SetOutputSpacing( space );
		resampler->SetOutputOrigin( origin );
		resampler->SetSize( size );
		resampler->SetDefaultPixelValue( 100 );
		resampler->SetInput( itkImage );
		//typedef itk::LinearInterpolateImageFunction<InputImageType,double>InterpolatorType;
		typedef itk::BSplineInterpolateImageFunction<InputImageType>InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		typedef itk::TranslationTransform<double, VImageDimension>  TransformType;

		TransformType::Pointer transform = TransformType::New();
		resampler->SetInterpolator( interpolator );
		resampler->SetTransform( transform );
		try{
			resampler->Update();
			std::cout << "type is unsiged char and Updated!" << std::endl;
		}catch(itk::ExceptionObject & e)
		{
			std::cout << "exception occured!" << e.what() << std::endl;
			system("pause");
		}
		OuputImageType * outputImage = resampler->GetOutput();
		*pointer = mitk::ImportItkImage(outputImage);
		mitk::Image::Pointer resultImage = mitk::ImportItkImage(outputImage);
		mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
		newNode->SetData(resultImage);
		newNode->SetProperty("name", mitk::StringProperty::New("interpolate image"));
		newNode->SetProperty("opacity", mitk::FloatProperty::New(0.5));
		mitk::DataStorage::GetInstance()->Add( newNode );
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	}
	else if(originalPixelType == typeid(short))
	{
		typedef itk::Image< short, VImageDimension > OuputImageType;
		typedef itk::ResampleImageFilter< InputImageType,OuputImageType,double>  FilterType;
		FilterType::Pointer resampler = FilterType::New();
		resampler->SetOutputSpacing( space );
		resampler->SetOutputOrigin( origin );
		resampler->SetSize( size );
		resampler->SetDefaultPixelValue( 100 );
		resampler->SetInput( itkImage );
		//typedef itk::LinearInterpolateImageFunction<InputImageType,double>InterpolatorType;
		typedef itk::BSplineInterpolateImageFunction<InputImageType>InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		typedef itk::TranslationTransform<double, VImageDimension>  TransformType;

		TransformType::Pointer transform = TransformType::New();
		resampler->SetInterpolator( interpolator );
		resampler->SetTransform( transform );
		try{
			resampler->Update();
			std::cout << "pixel type is short and Updated!" << std::endl;
		}catch(itk::ExceptionObject & e)
		{
			std::cout << "exception occured!" << e.what() << std::endl;
			system("pause");
		}
		OuputImageType * outputImage = resampler->GetOutput();
		*pointer = mitk::ImportItkImage(outputImage);
		mitk::Image::Pointer resultImage = mitk::ImportItkImage(outputImage);
		mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
		newNode->SetData(resultImage);
		newNode->SetProperty("name", mitk::StringProperty::New("interpolate image"));
		newNode->SetProperty("opacity", mitk::FloatProperty::New(0.5));
		mitk::DataStorage::GetInstance()->Add( newNode );
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();

	}
	else 
	{   
		char buf[128];
		sprintf(buf, "%s not supported by this method yet!", originalPixelType.GetItkTypeAsString().c_str());
		QMessageBox::information( NULL, "Virtual Surgery functionality", buf);
	}


	
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::CurvatureFlow(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	std::cout<< "Function Curvature Flow Image Start..." <<std::endl;
	typedef itk::Image< TPixel, VImageDimension >  InternalImageType;
	typedef itk::CurvatureFlowImageFilter< InternalImageType, InternalImageType >   CurvatureFlowImageFilterType;
	CurvatureFlowImageFilterType::Pointer smoothing =  CurvatureFlowImageFilterType::New();
	smoothing->SetInput(itkImage);
	smoothing->SetNumberOfIterations(60);
	smoothing->SetTimeStep(0.25);
	InternalImageType * outputImage = smoothing->GetOutput();
	*pointer = mitk::ImportItkImage(outputImage);
	std::cout<< "Function Curvature Flow Image End..." <<std::endl;
}

void QmitkAnnotation::SegmentByMiddleHepaticVein()
{
	std::cout<<"SegmentByMiddleHepaticVein Begin"<<std::endl;
	m_RainbowColor.GoToBegin();
	clipPolyData = "*Main Part*";
	clipImageData = "Main Image";
	mitk::DataTreeNode *node = mitk::DataStorage::GetInstance()->GetNamedNode(clipPolyData);
	vtkPolyData *source = ((mitk::Surface *)node->GetData())->GetVtkPolyData();
	if (source==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "No Polydata to Clip!");
		return;
	} 
	btnClip->setEnabled(true);
	btnMiddleHepaticVein->setEnabled(false);
	if(pWidget ==NULL)
	{

		QMessageBox::information( NULL, "Virtual Surgery functionality", "No ImplicitPlaneWidget.");
		return;	
	}

	std::cout<<"SegmentByMiddleHepaticVein Begin"<<std::endl;

}


void QmitkAnnotation::SegmentByLeftHepaticVein()
{
	
	std::cout<<"SegmentByLeftHepaticVein Begin!"<<std::endl;
	clipPolyData = "*Left Lobe*";
	clipImageData  = "Left Lobe";
	mitk::DataTreeNode *node = mitk::DataStorage::GetInstance()->GetNamedNode(clipPolyData);
	vtkPolyData *source = ((mitk::Surface *)node->GetData())->GetVtkPolyData();
	if (source==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "No Polydata to Clip!");
		return;
	} 

	btnClip->setEnabled(true);
	btnLeftHepaticVein->setEnabled(false);

	if(pWidget ==NULL)
	{

		QMessageBox::information( NULL, "Virtual Surgery functionality", "No ImplicitPlaneWidget.");
		return;	
	}
	std::cout<<"SegmentByLeftHepaticVein End!"<<std::endl;
}


void QmitkAnnotation::SegmentByRightHepaticVein()
{
	std::cout<<"SegmentByRightHepaticVein Begin"<<std::endl;
	clipPolyData = "*Right Lobe*";
	clipImageData = "Right Lobe";
	mitk::DataTreeNode *node = mitk::DataStorage::GetInstance()->GetNamedNode(clipPolyData);
	vtkPolyData *source = ((mitk::Surface *)node->GetData())->GetVtkPolyData();
	if (source==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "No Polydata to Clip!");
		return;
	} 
	
	btnClip->setEnabled(true);
	btnRightHepaticVein->setEnabled(false);


	if(pWidget ==NULL)
	{

		QMessageBox::information( NULL, "Virtual Surgery functionality", "No ImplicitPlaneWidget.");
		return;	
	}
	std::cout<<"SegmentByRightHepaticVein End!"<<std::endl;
}

void QmitkAnnotation::SegmentISeg()
{
	std::cout<<"SegmentISeg Begin"<<std::endl;

	this->isInitDataTree = true;


	clipPolyData = "*Liver Polydata*";
	clipImageData = "Liver Image";
	mitk::DataTreeNode *node = mitk::DataStorage::GetInstance()->GetNamedNode("*Liver Polydata*");
	vtkPolyData *source = ((mitk::Surface *)node->GetData())->GetVtkPolyData();
	if (source==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "No Polydata to Clip!");
		return;
	} 

	btnClip->setEnabled(true);
	btnISeg->setEnabled(false);

	if(pWidget == NULL)
	{
		vtkImplicitPlaneWidget *planeWidget = vtkImplicitPlaneWidget::New();
		planeWidget->SetInteractor(m_pRenWin->GetInteractor());
		planeWidget->SetPlaceFactor(1.25);
		planeWidget->GetPlaneProperty()->SetOpacity(1.0);
		planeWidget->GetOutlineProperty()->SetColor(0,0,1);
		planeWidget->GetPlaneProperty()->SetColor(0, 1.0, 0.8);
		double *origin = actor->GetCenter();
		planeWidget->SetOrigin(origin[0],origin[1],origin[2]);
		
		if(plane == NULL)
		{
			plane = vtkPlane::New();
			plane->SetOrigin(origin[0],origin[1],origin[2]);
			plane->SetNormal(1,1,1);
		}
		vtkTIPWCallback *pCaller = new vtkTIPWCallback();
		pCaller->plane = plane; 
		planeWidget->AddObserver(vtkCommand::InteractionEvent,pCaller);

		pWidget = planeWidget;
	}

	if (source!=NULL)
	{
		pWidget->SetInput(source);
	} 
	else
	{
		std::cout<<"No polydata!"<<std::endl;
		return;
	}
	pWidget->PlaceWidget();
	pWidget->On();

	std::cout<<"SegmentISeg End!"<<std::endl;
}

void QmitkAnnotation::SegmentByPortalVein()
{
	std::cout<<"SegmentByPortalVein Begin"<<std::endl;
	clipPolyData = "*Portal Vein*";
	clipImageData = "Portal Vein";
	btnPortalVein->setEnabled(false);
	btnClip->setEnabled(true);

	if(pWidget ==NULL)
	{

		QMessageBox::information( NULL, "Virtual Surgery functionality", "No ImplicitPlaneWidget.");
		return;	
	}
	this->isInitDataTree = true;
	std::cout<<"SegmentByPortalVein End!"<<std::endl;
}

void QmitkAnnotation::Clip()
{
	std::cout<<"---------------Clip Begin---------------"<<std::endl;
	if (plane==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "Not Specific Clipping Plane!");
		return;
	} 
	double normal[3] = {plane->GetNormal()[0],plane->GetNormal()[1],plane->GetNormal()[2]};
	double origin[3] = {plane->GetOrigin()[0],plane->GetOrigin()[1],plane->GetOrigin()[2]};

	for (int i=0;i<3;i++)
	{
		std::cout << normal[i] <<" "<< origin[i] <<std::endl;
	}
	//======================================clip poly data================================================//
	char *nodeName1,*nodeName2;
	char *nodeName3,*nodeName4;
	if(clipPolyData=="*Portal Vein*")
	{
		mitk::DataTreeNode *anteriorLobeNode = mitk::DataStorage::GetInstance()->GetNamedNode("*Anterior Lobe*");
		mitk::DataTreeNode *posteriorLobeNode = mitk::DataStorage::GetInstance()->GetNamedNode("*Posterior Lobe*");
		mitk::DataTreeNode *siphonalLobeNode = mitk::DataStorage::GetInstance()->GetNamedNode("*Siphonal Lobe*");
		vtkPolyData *anteriorLobeSource = ((mitk::Surface *)anteriorLobeNode->GetData())->GetVtkPolyData();
		vtkPolyData *posteriorLobeSource = ((mitk::Surface *)posteriorLobeNode->GetData())->GetVtkPolyData();
		vtkPolyData *siphonalLobeSource = ((mitk::Surface *)siphonalLobeNode->GetData())->GetVtkPolyData();
		vtkPolyData *part1 = vtkPolyData::New();
		vtkPolyData *part2 = vtkPolyData::New();
		ClipPolydata(normal,origin,anteriorLobeSource,part1,part2);
		setNode(part1,anteriorLobeNode,"*V Segment*");
		setNode(part2,anteriorLobeNode,"*VIII Segment*");

		vtkPolyData *part3 = vtkPolyData::New();
		vtkPolyData *part4 = vtkPolyData::New();
		ClipPolydata(normal,origin,posteriorLobeSource,part3,part4);
		setNode(part3,posteriorLobeNode,"*VI Segment*");
		setNode(part4,posteriorLobeNode,"*VII Segment*");

		vtkPolyData *part5 = vtkPolyData::New();
		vtkPolyData *part6 = vtkPolyData::New();
		ClipPolydata(normal,origin,siphonalLobeSource,part5,part6);
		setNode(part5,siphonalLobeNode,"*III Segment*");
		setNode(part6,siphonalLobeNode,"*II Segment*");

		/****clip image****/
		mitk::DataTreeNode *anteriorLobeImageNode = mitk::DataStorage::GetInstance()->GetNamedNode("Anterior Lobe");
		mitk::DataTreeNode *posteriorLobeImageNode = mitk::DataStorage::GetInstance()->GetNamedNode("Posterior Lobe");
		mitk::DataTreeNode *siphonalLobeImageNode = mitk::DataStorage::GetInstance()->GetNamedNode("Siphonal Lobe");
		mitk::Image *anteriorLobeImageSource = dynamic_cast<mitk::Image*>(anteriorLobeImageNode->GetData());//前叶
		mitk::Image *posteriorLobeImageSource = dynamic_cast<mitk::Image*>(posteriorLobeImageNode->GetData());//后页
		mitk::Image *siphonalLobeImageSource = dynamic_cast<mitk::Image*>(siphonalLobeImageNode->GetData());//外页
	
		ClipImagedata(normal,origin,anteriorLobeImageSource,anteriorLobeImageNode,"V Segment","VIII Segment");
		ClipImagedata(normal,origin,posteriorLobeImageSource,posteriorLobeImageNode,"VI Segment","VII Segment");
		ClipImagedata(normal,origin,siphonalLobeImageSource,siphonalLobeImageNode,"III Segment","II Segment");
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	}
	else 
	{
		mitk::DataTreeNode *node = mitk::DataStorage::GetInstance()->GetNamedNode(clipPolyData);
		vtkPolyData *source = ((mitk::Surface *)node->GetData())->GetVtkPolyData();
		mitk::DataTreeNode *clipImageNode = mitk::DataStorage::GetInstance()->GetNamedNode(clipImageData);
		mitk::Image  *clipImage = dynamic_cast<mitk::Image*>(clipImageNode->GetData());
		if ( source==NULL || clipImage == NULL )
		{
			QMessageBox::information( NULL, "Virtual Surgery functionality", "No data to Clip!");
			return;
		} 
		vtkPolyData *part1 = vtkPolyData::New();
		vtkPolyData *part2 = vtkPolyData::New();
		mitk::Image *imagePart1 = mitk::Image::New();
		mitk::Image *imagePart2 = mitk::Image::New();
		//planeWidget->GetPlane(plane);
		ClipPolydata(normal,origin,source,part1,part2);
		if (clipPolyData == "*Main Part*"){nodeName1="*Right Lobe*";nodeName2="*Left Lobe*";btnLeftHepaticVein->setEnabled(true);}
		else if(clipPolyData == "*Right Lobe*"){nodeName1="*Anterior Lobe*";nodeName2="*Posterior Lobe*";btnPortalVein->setEnabled(true);}
		else if(clipPolyData == "*Left Lobe*"){nodeName1="*IV Segment*";nodeName2="*Siphonal Lobe*";btnRightHepaticVein->setEnabled(true);}
		else if(clipPolyData=="*Liver Polydata*"){nodeName1="*I Segment*";nodeName2="*Main Part*";btnMiddleHepaticVein->setEnabled(true);}
		setNode(part1,node,nodeName1);
		setNode(part2,node,nodeName2);
		
		
		if (clipImageData == "Main Image"){nodeName3="Right Lobe";nodeName4="Left Lobe";}
		else if(clipImageData == "Right Lobe"){nodeName3="Anterior Lobe";nodeName4="Posterior Lobe";}
		else if(clipImageData == "Left Lobe"){nodeName3="IV Segment";nodeName4="Siphonal Lobe";}
		else if(clipImageData=="Liver Image"){nodeName3="I Segment";nodeName4="Main Image";}
		ClipImagedata(normal,origin,clipImage,clipImageNode,nodeName3,nodeName4);
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	}
	//==========================================================================================================================//

	btnClip->setEnabled(false);
	
	//========================Init Segemnt DataTree======================//
	if(isInitDataTree)
	{
		InitDataTree();
	}
	//===================================================================//

	std::cout<<"---------------Clip End-----------------"<<std::endl;
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::GetMaxMinPixelIntensity(itk::Image<TPixel, VImageDimension> *itkImage, int &maxIntensity, int &minIntensity)
{
	typedef itk::Image<TPixel, VImageDimension> InputImageType;
	typedef itk::ImageRegionIterator<InputImageType> RegionIteratorType;

	InputImageType::RegionType region    = itkImage->GetLargestPossibleRegion();

	maxIntensity = -3024;
	minIntensity =  3024;

	RegionIteratorType it(itkImage, itkImage->GetLargestPossibleRegion() );
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > maxIntensity)
			maxIntensity = it.Get();
		if(it.Get() < minIntensity)
			minIntensity = it.Get();
	}
}

void QmitkAnnotation::ClipImagedata(double (&normal)[3],double (&origin)[3], mitk::Image *image,mitk::DataTreeNode *parent,char *name1,char *name2)
{
	std::cout << "Clip Image Data begin..." <<std::endl;
	mitk::PlaneCutFilter::Pointer cutter = mitk::PlaneCutFilter::New();
	cutter->SetFillMode(mitk::PlaneCutFilter::FillMode::FILL);
	cutter->SetBackgroundLevel(-3024);
	mitk::PlaneGeometry::Pointer cutPlane = mitk::PlaneGeometry::New();

	float nor[3],ori[3];
	for (int i=0;i<3;i++)
	{
		nor[i] = -float(normal[i]);//strange
		ori[i] = float(origin[i]);
	}
	mitk::Vector3D planeNormal(nor);
	mitk::Point3D planeOrigin(ori);

	cutPlane->InitializePlane(planeOrigin,planeNormal);
	cutter->SetPlane(cutPlane);
	cutter->SetInput( image );
	try
	{
		cutter->UpdateLargestPossibleRegion();
	}
	catch(itk::ExceptionObject&)
	{
		QMessageBox::warning ( NULL,
			tr("Cutting not possible"),
			tr("Sorry, the bounding box has to be completely inside the image.\n\n"
			"The possibility to drag it larger than the image a bug and has to be fixed."),
			QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
		return;
	}
	mitk::Image *one = cutter->GetOutput();
	addImageNode(one,parent,name1);

	mitk::PlaneCutFilter::Pointer cutter1 = mitk::PlaneCutFilter::New();
	cutter1->SetFillMode(mitk::PlaneCutFilter::FillMode::FILL);
	cutter1->SetBackgroundLevel(-3024);
	mitk::PlaneGeometry::Pointer cutPlane1 = mitk::PlaneGeometry::New();
	for (int i=0;i<3;i++)
	{
		nor[i] = float(normal[i]);
		ori[i] = float(origin[i]);
	}
	mitk::Vector3D planeNormal1(nor);
	mitk::Point3D planeOrigin1(ori);
	cutPlane1->InitializePlane(planeOrigin1,planeNormal1);
	cutter1->SetPlane(cutPlane1);
	cutter1->SetInput( image );
	try
	{
		cutter1->UpdateLargestPossibleRegion();
	}
	catch(itk::ExceptionObject&)
	{
		QMessageBox::warning ( NULL,
			tr("Cutting not possible"),
			tr("Sorry, the bounding box has to be completely inside the image.\n\n"
			"The possibility to drag it larger than the image a bug and has to be fixed."),
			QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
		return;
	}
	mitk::Image *another = cutter1->GetOutput();
	addImageNode(another,parent,name2);
	
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	std::cout << "Clip Image Data end..." <<std::endl;
}

void QmitkAnnotation::setNode(vtkPolyData *polydata,mitk::DataTreeNode *parent,char *nodeName)
{
	mitk::DataTreeNode::Pointer node = mitk::DataTreeNode::New();
	mitk::Surface::Pointer surface = mitk::Surface::New();
	surface->SetVtkPolyData(polydata);
	node->SetData(surface);
	node->SetProperty("name", mitk::StringProperty::New(nodeName));
	node->SetProperty("opacity", mitk::FloatProperty::New(1.0));
	node->SetProperty("Surface", mitk::BoolProperty::New(true));
	std::string name(nodeName);
	
	
	if(name.find("Segment")!=name.npos)
	{
		node->SetProperty("color",mitk::ColorProperty::New(m_RainbowColor.GetNextColor()));
	}

	node->SetVisibility(true);
	mitk::DataStorage::GetInstance()->Add( node,parent );
}


void QmitkAnnotation::ClipPolydata(double (&normal)[3],double (&origin)[3],vtkPolyData *source,vtkPolyData *one,vtkPolyData *another)
{
	std::cout<<"Clip Begin."<<std::endl;
	if (source==NULL)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "No Polydata to Clip!");
		return;
	} 
	vtkPlane *plane = vtkPlane::New();
	plane->SetNormal(normal);
	plane->SetOrigin(origin);
	vtkPlane *plane1 = vtkPlane::New();
	plane1->SetNormal(0-normal[0],0-normal[1],0-normal[2]);
	plane1->SetOrigin(origin);
	vtkClipPolyData *clipper = vtkClipPolyData::New();
	clipper->SetInput(source);
	clipper->SetClipFunction(plane);
	clipper->SetValue(0.2);
	clipper->Update();
	one->DeepCopy(clipper->GetOutput());

	clipper->SetClipFunction(plane1);
	clipper->SetValue(0.2);
	clipper->Update();
	another->DeepCopy(clipper->GetOutput());

	std::cout<<"Clip End."<<std::endl;

}

void QmitkAnnotation::addImageNode(mitk::Image *data,mitk::DataTreeNode *parent,char *nodeName)
{
	mitk::DataTreeNode::Pointer node = mitk::DataTreeNode::New();
	node->SetData(data);
	node->SetProperty("name", mitk::StringProperty::New(nodeName));
	node->SetProperty("opacity", mitk::FloatProperty::New(0.0));
	node->SetProperty("Surface", mitk::BoolProperty::New(true));
	mitk::DataStorage::GetInstance()->Add( node,parent );
}


void QmitkAnnotation::EnableAxes()
{
	bool isAxes = Isaxes->isChecked();
	std::cout << "Axes : " << isAxes <<std::endl;
	if (axes == NULL && isAxes == true)
	{
		QMessageBox::information( NULL, "Virtual Surgery functionality", "Have not render yet!");
		return;
	} 

	if(axes!=NULL)
	{
		axes->SetVisibility(isAxes);
	}
}

void QmitkAnnotation::CaclwithAnnotation()
{
	if (voxelNumberOfImage == 0)
	{
		mitk::DataTreeNode *imageNode = mitk::DataStorage::GetInstance()->GetNamedNode("Liver Image");
		mitk::Image  *image = dynamic_cast<mitk::Image*>(imageNode->GetData());
		if ( image == NULL )
		{
			QMessageBox::information( NULL, "Virtual Surgery functionality", "No data to Clip!");
			return;
		} 
		long temp;
		AccessByItk_2(image, VoxelNumberOfImage, voxelNumberOfImage, temp);
	}


	std::cout << "CaclwithAnnotation begin ..." << std::endl;
	if(selectedImage.IsNull())
	{
		QMessageBox::warning ( NULL,
			tr("erode vessel not possible"),
			tr("Sorry, you must select a image.\n\n"
			"The possibility to not load an image."),
			QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
		return;
	}
	mitk::DataTreeNode* node = selectedImage->Get();
	if ( node )
	{
		mitk::BaseData* data = node->GetData();
		if (data)
		{
			mitk::Image* image = dynamic_cast<mitk::Image*>( data );
			if (image)
			{
				long param1 = 0;
				long param2 = 0;
				float ratio = 0;
				AccessByItk_2(image, VoxelNumberOfImage, param1, param2);
				ratio = (float)param1/voxelNumberOfImage;
				std::cout <<param1<<" "<<param2<<" "<<ratio<< std::endl;
				QString numberStr,ratioStr;
				numOfPixels->setText(numberStr.sprintf("%ld",param1));
				percentInfo->setText(ratioStr.sprintf("%0.5f",(float)param1/voxelNumberOfImage));

				mitk::VolumeCalculator::Pointer volCalc = mitk::VolumeCalculator::New();
				volCalc->SetImage(image);
				volCalc->ComputeVolume();
				std::stringstream vs;
				vs << volCalc->GetVolume() << " ml";
				volumeInfo->setText(vs.str().c_str() );
				

				//===================input to file====================//
				ofstream outfile("D:\\annotaionInformation.txt",ofstream::app);
				outfile <<setiosflags(ios::left)<< "Image:"<<node->GetName()<<'\t'<<"Voxel Number:"<<param1<<'\t'<<"Volume:"<<vs.str()
					<<'\t'<<"Percent Ratio:"<<setprecision(5)<<ratio<<endl;
				outfile.flush();
				outfile.close();
				//====================================================//
			}
		}
	}

	
	std::cout << "CaclwithAnnotation end ..." << std::endl;
}

void QmitkAnnotation::InitDataTree()
{
	std::cout << "Init DataTree Selector start..." << std::endl;
	m_TreeNodeSelector->SetDataTreeNodeIterator(m_DataTreeIteratorBase);
	//m_TreeNodeSelector->GetFilter()->SetFilter(mitk::IsSegmentNode());
	std::cout << "Init DataTree Selector over..." << std::endl;
}


template <typename ImageType>
void QmitkAnnotation::ConnectVTKToITK(vtkImageExport* in, itk::VTKImageImport<ImageType>* out)
{
	out->SetUpdateInformationCallback(in->GetUpdateInformationCallback());
	out->SetPipelineModifiedCallback(in->GetPipelineModifiedCallback());
	out->SetWholeExtentCallback(in->GetWholeExtentCallback());
	out->SetSpacingCallback(in->GetSpacingCallback());
	out->SetOriginCallback(in->GetOriginCallback());
	out->SetScalarTypeCallback(in->GetScalarTypeCallback());
	out->SetNumberOfComponentsCallback(in->GetNumberOfComponentsCallback());
	out->SetPropagateUpdateExtentCallback(in->GetPropagateUpdateExtentCallback());
	out->SetUpdateDataCallback(in->GetUpdateDataCallback());
	out->SetDataExtentCallback(in->GetDataExtentCallback());
	out->SetBufferPointerCallback(in->GetBufferPointerCallback());
	out->SetCallbackUserData(in->GetCallbackUserData());
}

void QmitkAnnotation::ImageSelected(mitk::DataTreeIteratorClone imageIt)
{
	assert( imageIt.IsNotNull() ); // should never fail, the selection widget cares for that
	std::cout << " Image Selected " << std::endl;
	mitk::DataTreeNode* node = imageIt->Get();
	selectedImage = imageIt;
	if ( node )
	{
		std::string name;
		if (node->GetName(name))
		{
			std::cout << "Tree node selected with name '" << name << "'" << std::endl;
		}
		mitk::BaseData* data = node->GetData();
		if (data)
		{
			mitk::Image* image = dynamic_cast<mitk::Image*>( data );
			if (image)
			{
				std::cout << "Surprise: this node contains a real image dataset." << std::endl;
				originalPixelType = image->GetPixelType();
			}
		}
	}
}

void QmitkAnnotation::TreeChanged()
{
	m_TreeNodeSelector->SetDataTreeNodeIterator(m_DataTreeIteratorBase);
}

template <typename TPixel, unsigned int VImageDimension>
long QmitkAnnotation::VoxelNumberOfImage(itk::Image<TPixel, VImageDimension> *itkImage, long &param1, long &param2)
{
	long numberOfVoxel  = 0;
	typedef itk::Image<TPixel, VImageDimension> ImageType; 
	typedef itk::ImageRegionIterator<ImageType> RegionIteratorType;
	//ImageType::RegionType region    = itkImage->GetLargestPossibleRegion();
	RegionIteratorType it(itkImage, itkImage->GetLargestPossibleRegion() );
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if(it.Get() > 0)
			param1++;
		else
			param2++;

	}
	return numberOfVoxel;
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::CopyImage(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	/***************duplicate a image to clip 8 segment******************/
	typedef itk::Image<TPixel, VImageDimension> TImageType;
	typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage( itkImage );
	duplicator->Update();
	typename TImageType::Pointer clonedImage = duplicator->GetOutput();
	mitk::Image::Pointer resultImage = mitk::ImportItkImage( clonedImage );
	mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
	newNode->SetData(resultImage);
	newNode->SetProperty("name", mitk::StringProperty::New("Liver Image"));
	newNode->SetProperty("opacity", mitk::FloatProperty::New(0.0));
	mitk::DataStorage::GetInstance()->Add( newNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	/*********************************************************************/
}


//use vertex set and edge set to display vascular tree
void QmitkAnnotation::DisplayVesselTree(std::vector<Vertex> vs,std::vector<Edge> es)
{
	
	vtkPoints *points = vtkPoints::New();
	vtkFloatArray *scalar  = vtkFloatArray::New();
	grid = vtkUnstructuredGrid::New();
	std::vector<Vertex>::iterator beg;
	//convert point index to graph index
	for (size_t index= 0;index<vs.size();index++)
	{
		//pts->InsertPoint(vertexIndexMap[*beg],vesselGraph->txMap[*beg],vesselGraph->tyMap[*beg],vesselGraph->tzMap[*beg]);
		points->InsertPoint(index,vesselGraph->txMap[vs.at(index)],vesselGraph->tyMap[vs.at(index)],vesselGraph->tzMap[vs.at(index)]);
		std::cout << index<<" : (" <<vesselGraph->txMap[vs.at(index)]<<","
			<<vesselGraph->tyMap[vs.at(index)]<<","<<vesselGraph->tzMap[vs.at(index)]<<")"<<std::endl;
		vtkVertex *vertex = vtkVertex::New();
		vertex->GetPointIds()->SetId(0,index);
		scalar->InsertNextValue(1.0);
		grid->InsertNextCell(vertex->GetCellType(),vertex->GetPointIds());
	}
	grid->SetPoints(points);
	grid->GetPointData()->SetScalars(scalar);
	vtkCellArray *lines=vtkCellArray::New();
	std::vector<Edge>::iterator e_beg;
	int src,tar;
	for (e_beg = es.begin();e_beg!= es.end();e_beg++)
	{
		src = tar = -1;
		vtkLine *line = vtkLine::New();
		for(size_t i= 0;i<vs.size();i++)
		{
			if(vesselGraph->GetVertexIndex(vs.at(i))==vesselGraph->GetEdgeSourceIndex(*e_beg))
			{
				src = i;
			}
			if(vesselGraph->GetVertexIndex(vs.at(i))==vesselGraph->GetEdgeTargetIndex(*e_beg))
			{
				tar = i;	
			}
		}
		std::cout << "("<<src<<","<<tar<<")"<<std::endl;
		if(src>=0 && tar>=0)
		{
			line->GetPointIds()->SetId(0,src);
			line->GetPointIds()->SetId(1,tar);
			lines->InsertNextCell(line);
		}
	}

	vtkPolyData *treeData = vtkPolyData::New();
	treeData->SetPoints(points);
	treeData->SetLines(lines);

	vtkPolyDataMapper *treeMapper = vtkPolyDataMapper::New();
	treeMapper->ScalarVisibilityOff();
	treeMapper->SetInput(treeData);
	vtkActor *treeActor = vtkActor::New();
	treeActor->PickableOff();
	treeActor->SetMapper(treeMapper);
	treeActor->GetProperty()->SetColor(1.0,1.0,1.0);

	vtkSphereSource *sphere = vtkSphereSource::New();
	sphere->SetRadius(0.70);
	sphere->SetPhiResolution(20);
	sphere->SetThetaResolution(20);

	vtkThresholdPoints *thresholdIn = vtkThresholdPoints::New();
	thresholdIn->SetInput(grid);
	thresholdIn->ThresholdByUpper(0.5);
	vtkGlyph3D *vertices = vtkGlyph3D::New();
	vertices->SetInputConnection(thresholdIn->GetOutputPort());
	vertices->SetSource(sphere->GetOutput());
	vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
	sphereMapper->SetInputConnection(vertices->GetOutputPort());
	sphereMapper->ScalarVisibilityOff();
	vtkActor *cubeActor = vtkActor::New();
	cubeActor->SetMapper(sphereMapper);
	cubeActor->GetProperty()->SetColor(0.0,1.0,0.0);
	cubeActor->PickableOff();

	vtkDataSetMapper *verMapper = vtkDataSetMapper::New();
	verMapper->SetInput(grid);
	vtkActor *verActor = vtkActor::New();
	verActor->SetMapper(verMapper);
	verActor->GetProperty()->SetColor(0.0,1.0,1.0);

	m_pRenderer->AddActor(treeActor);
	m_pRenderer->AddActor(verActor);
	m_pRenderer->AddActor(cubeActor);
	m_pRenderer->ResetCamera();
	mitk::RenderingManager::GetInstance()->AddRenderWindow(m_pRenWin);
	m_pRenWin->Render();

	mitk::vtkPickCallback *pickObsever  = mitk::vtkPickCallback::New();
	picker->SetTolerance(0.02);
	picker->AddObserver(vtkCommand::EndPickEvent,pickObsever);
	m_pRenWin->GetInteractor()->SetPicker(picker);
	m_pRenWin->Render();

}

// for automatic vascular branch segemnt
void QmitkAnnotation::DoDisplayVesselTree()
{
	std::cout << "DoVesselDivision begin..." << std::endl;;
	QFileDialog* fd = new QFileDialog( this, "File Dialog", true );
	fd->setMode(QFileDialog::ExistingFile);
	fd->addFilter( "Xmls (*.xml)" );
	fd->show();

	QString fileName;
	if ( fd->exec() == QDialog::Accepted )
		fileName = fd->selectedFile();
	std::cout << fileName.ascii() << std::endl;

	if(!vesselGraph)
		vesselGraph = new VesselGraph(fileName.ascii());
	vesselGraph->SetSubGraphCount(1);
	m_RainbowColor.GoToBegin();
	//DisplayVesselTree(vesselGraph->vs,vesselGraph->es);
	/*
	Function code of DisplayVesselTree.
	*/
	vesselGraph->RadiusFilter(2.30);
	DisplayVesselTree(vesselGraph->vs,vesselGraph->es);
	
	vesselGraph->KMeans(8,0.5);
	DisplaySubtreeNodes();


	m_btnNewSubTree->setEnabled(true);
	m_btnDisplayVesselTree->setEnabled(false);
	segListView->setEnabled(true);
}


void QmitkAnnotation::CreateNewSubTree()
{
	int numOfSeg = segListView->childCount();
	if(numOfSeg == segExp)
	{
		QMessageBox::information( NULL, "Vessel Tree segmentation", "Eight subtrees are already created, please set seeds!");
		m_btnNewSubTree->setEnabled(false);
		m_btnManuLiverSegment->setEnabled(true);
		return;
	}
	else
	{
		int index = vesselGraph->GetSubGraphCount();
		std::cout <<"index:"<< index << std::endl;
		if ((index-1)%7==0)
		{
			m_RainbowColor.GoToBegin();
		}
		color = m_RainbowColor.GetNextColor();
		float rgb[3]; rgb[0] = color.GetRed(); rgb[1] = color.GetGreen(); rgb[2] = color.GetBlue();;
		QRgb qrgb = qRgb( ROUND(rgb[0]*255.0), ROUND(rgb[1]*255.0), ROUND(rgb[2]*255.0) );

		QColor qtcolor(qrgb);
		QString text = QString("%1").arg( index );
		QPixmap pixmap(25,18);
		pixmap.fill(qrgb);
		QPainter painter( &pixmap );
		QPen pen = painter.pen();
		int ha,es,vau;
		qtcolor.getHsv(ha,es,vau);
		if ( vau < 160  )
			pen.setColor ( Qt::white );
		else
			pen.setColor ( Qt::black );
		painter.setPen( pen );
		QFontMetrics fm = painter.fontMetrics();
		QRect bb = fm.boundingRect( text );
		//painter.drawText( (25 - bb.width()) / 2,(18 + bb.height()) / 2,text, 0, -1 );
		
		QListViewItem* lastItem = segListView->lastItem();
		QListViewItem* newItem = new QListViewItem(segListView);
		newItem->setPixmap(0,pixmap);
		QString name = QString("%1 subtree").arg( index );
		newItem->setText(1,name);
		std::cout << segListView->columns()<< std::endl;

		m_btnNewSubTree->setEnabled(false);
		m_btnDivideTree->setEnabled(true);
	}
	
}

void QmitkAnnotation::itemSelected(QListViewItem* item)
{
	if(!item)
	{
		return;
	}
	int pos = atoi(item->text(1).left(1).ascii());
	std::cout << pos << std::endl;
	if(item->isSelected())
	{
		if(segMap[pos]!=NULL)
		{
			segMap[pos]->VisibilityOn();
			m_pRenWin->Render();
		}
	}
	else
	{
		if(segMap[pos]!=NULL)
		{
			segMap[pos]->VisibilityOff();
			m_pRenWin->Render();
		}
	}
	
}

void QmitkAnnotation::actorItemSelected(QListViewItem* item)
{
	std::map<const QListViewItem*,vtkActor*>::iterator iter;
	iter = actorMap.find(item);
	if(iter!=actorMap.end())
	{
		this->actor = iter->second;
		activeActorProperty(true);
	}
	else
	{
		activeActorProperty(false);
	}
}

void QmitkAnnotation::DoVesselDivision()
{
	if(!vesselGraph)
	{
		return;
	}

	vtkProp *prop = picker->GetViewProp();
	if (prop != NULL)
	{
		int iSel = picker->GetCellId();//获得选择的细胞编号
		vtkActor* pA = (vtkActor*)picker->GetActor();
		vtkIdList *ptlds = vtkIdList::New();
		pA->GetMapper()->GetInput()->GetCellPoints(iSel,ptlds);
		vtkIdType start;
		for (int i=0;i<ptlds->GetNumberOfIds();i++)
		{
			start = ptlds->GetId(i);
		}
		vesselGraph->SetStartVertex(start);
		if(!vesselGraph->GetSubGraph())
		{
			return;	
		}

		//display another color for new subgraph
		OutEdgeIterator e_beg,e_end;
		std::list<Vertex>::iterator begin;
		vtkPoints *points = vtkPoints::New();	
		int index = vesselGraph->GetSubGraphCount();
		vtkCellArray *lines=vtkCellArray::New();
		for (begin = vesselGraph->QDiVertex[index].begin();begin!=vesselGraph->QDiVertex[index].end();begin++)
		{
			points->InsertPoint(vesselGraph->GetVertexIndex(*begin),vesselGraph->txMap[*begin],vesselGraph->tyMap[*begin],vesselGraph->tzMap[*begin]);
			for(tie(e_beg,e_end) = out_edges(*begin,vesselGraph->g);e_beg!=e_end;e_beg++)
			{
				vtkLine *line = vtkLine::New();
				line->GetPointIds()->SetId(0,vesselGraph->GetEdgeSourceIndex(*e_beg));
				line->GetPointIds()->SetId(1,vesselGraph->GetEdgeTargetIndex(*e_beg));
				lines->InsertNextCell(line);
			}
			
		}
		vtkPolyData *subTreeData = vtkPolyData::New();
		subTreeData->SetPoints(points);
		subTreeData->SetPolys(lines);
		subTreeData->Update();
		std::cout << subTreeData->GetNumberOfPoints() << " , " << subTreeData->GetNumberOfLines()<<std::endl;
		vtkExtractEdges *edges = vtkExtractEdges::New(); 
		edges->SetInput(subTreeData);
		edges->Update();
		vtkTubeFilter *tubes = vtkTubeFilter::New();
		tubes->SetInputConnection(edges->GetOutputPort());
		tubes->SetRadius(0.55);
		tubes->SetNumberOfSides(6);
		//tubes->UseDefaultNormalOn();
		//tubes->SetDefaultNormal(.577, .577, .577);
		tubes->Update();
		vtkPolyData *pdata = tubes->GetOutput();
		std::cout << pdata->GetNumberOfPoints() << " , " << pdata->GetNumberOfLines()<<std::endl;

		vtkPolyDataMapper *subTreeMapper = vtkPolyDataMapper::New();
		subTreeMapper->ScalarVisibilityOff();
		//subTreeMapper->SetInput(subTreeData);
		subTreeMapper->SetInput(pdata);
		vtkActor *subTreeActor = vtkActor::New();
		subTreeActor->PickableOff();
		subTreeActor->SetMapper(subTreeMapper);
		subTreeActor->GetProperty()->SetColor(color.GetRed(),color.GetGreen(),color.GetBlue());
		m_pRenderer->AddActor(subTreeActor);
		m_pRenWin->Render();

		m_btnNewSubTree->setEnabled(true);
		m_btnDivideTree->setEnabled(false);

		if(segListView->childCount()==segExp)
		{
			m_btnNewSubTree->setEnabled(false);
			m_btnManuLiverSegment->setEnabled(true);
		}
	}
}

void QmitkAnnotation::AutomaticLiverSegment()
{
	if(!vesselGraph)
	    vesselGraph = new VesselGraph("D:\\VesseltreeGraph.xml");
	double arg = atof(radiusRatio->text());
	std::cout << arg << std::endl;
	vesselGraph->SetSeparateArg(arg);
	vesselGraph->DivideVesselTree();
	vesselGraph->OutputGraph();
	
	DisplaySubtreeNodes();

	mitk::Image* image;
	mitk::Image::Pointer resultImage;
	QListViewItemIterator it(m_pInstance->listView);
	it++;
	while (it.current()) 
	{
		DataTreeViewItem *item = (DataTreeViewItem *)it.current();
		if (item->isOn())
		{
			std::cout<<item->text(0)<<std::endl;
			image = dynamic_cast<mitk::Image*> (item->GetDataTreeNode()->GetData());
			if (image)
			{
				if(vesselGraph)
				{	
					resultImage = vesselGraph->VoxelDivision(image);
				}
				else
				{	
					QMessageBox::warning ( NULL,
						tr("Sorry,should first divide veeesl"),tr(""),
						QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
				}
			}
			break;
		}
		it++;
	}

	mitk::LabeledImageToSurfaceFilter::Pointer filter = mitk::LabeledImageToSurfaceFilter::New();
	if (filter.IsNull())
	{
		std::cout<<"[FAILED]"<<std::endl;
		return;
	}
	if (resultImage.IsNotNull())
	{
		filter->SetInput(resultImage);
		filter->GenerateAllLabelsOn();
		filter->SetGaussianStandardDeviation( 1.5 ); 
		filter->SetSmooth(true);
		filter->SetDecimate( mitk::ImageToSurfaceFilter::DecimatePro );
		filter->SetTargetReduction( 0.2 );
		filter->SetBackgroundLabel(-10000);
		filter->Update();

		int numOfSurfaces = filter->GetNumberOfOutputs();
		std::cout << numOfSurfaces <<std::endl;
		mitk::ColorSequenceRainbow m_RainbowColor;
		if (numOfSurfaces>=1)
		{
			for (int index=0;index<numOfSurfaces;index++)
			{
				if (index%8==0)
				{
					m_RainbowColor.GoToBegin();
				}
				mitk::Surface::Pointer pSurface = filter->GetOutput(index);
				vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
				mapper->ScalarVisibilityOff();
				mapper->SetInput(pSurface->GetVtkPolyData());
				vtkActor *actor = vtkActor::New();
				actor->SetMapper(mapper);
				mitk::Color color =  m_RainbowColor.GetNextColor();
				actor->GetProperty()->SetColor(color.GetRed(),color.GetGreen(),color.GetBlue());
				m_pRenderer->AddActor(actor);
				m_pRenderer->ResetCamera();
				mitk::RenderingManager::GetInstance()->AddRenderWindow(m_pRenWin);
				m_pRenWin->Render();
			}
		}

		
	}
	else
	{
		QMessageBox::warning ( NULL,
			tr("Sorry,No label result image to Surface."),tr(""),
			QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
	}

}

void QmitkAnnotation::ManualLiverSegment()
{
	mitk::Image::Pointer resultImage;
	mitk::Image* image;
	mitk::Image::Pointer temp;
	QListViewItemIterator it(m_pInstance->listView);
	it++;
	while (it.current()) 
	{
		DataTreeViewItem *item = (DataTreeViewItem *)it.current();
		if (item->isOn())
		{
			image = dynamic_cast<mitk::Image*> (item->GetDataTreeNode()->GetData());
			if (image)
			{
				if(vesselGraph)
				{	
					resultImage = vesselGraph->VoxelDivision(image);
					AccessByItk_1(resultImage, SeparateLabelImage, &temp);
				}
				else
				{	
					QMessageBox::warning ( NULL,
						tr("Sorry,should first divide veeesl"),tr(""),
						QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
				}
			}
			break;
		}
		it++;
	}
	
	mitk::LabeledImageToSurfaceFilter::Pointer filter = mitk::LabeledImageToSurfaceFilter::New();
	if (filter.IsNull())
	{
		std::cout<<"[FAILED]"<<std::endl;
		return;
	}
	if (resultImage.IsNotNull())
	{
		filter->SetInput(resultImage);
		filter->GenerateAllLabelsOn();
		filter->SetGaussianStandardDeviation( 2.0 ); 
		filter->SetSmooth(true);
		filter->SetDecimate( mitk::ImageToSurfaceFilter::DecimatePro );
		filter->SetTargetReduction( 0.2 );
		filter->SetBackgroundLabel(-10000);
		filter->Update();

		int numOfSurfaces = filter->GetNumberOfOutputs();
		std::cout << numOfSurfaces <<std::endl;
		mitk::ColorSequenceRainbow m_RainbowColor;
		m_RainbowColor.GoToBegin();
		if (numOfSurfaces==vesselGraph->GetSubGraphCount()-1)
		{
			for (int index=0;index<numOfSurfaces;index++)
			{
				mitk::LabeledImageToSurfaceFilter::LabelType label = filter->GetLabelForNthOutput(index);
				std::cout << (label/10)-1 << std::endl;
				
				if (index%7==0)
				{
					m_RainbowColor.GoToBegin();
				}
				mitk::Surface::Pointer pSurface = filter->GetOutput(index);
				vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
				mapper->ScalarVisibilityOff();
				mapper->SetInput(pSurface->GetVtkPolyData());
				vtkActor *actor = vtkActor::New();
				actor->SetMapper(mapper);
				mitk::Color color =  m_RainbowColor.GetNextColor();
				actor->GetProperty()->SetColor(color.GetRed(),color.GetGreen(),color.GetBlue());
				m_pRenderer->AddActor(actor);
				actor->VisibilityOff();
				segMap.insert(std::make_pair((label/10)-1,actor));
				m_pRenderer->ResetCamera();
				mitk::RenderingManager::GetInstance()->AddRenderWindow(m_pRenWin);
				m_pRenWin->Render();

				QString str;
				str.sprintf("%d Segment",index+1);
				QCheckListItem *actorItem = new QCheckListItem(listView->firstChild(),str);
				QPixmap pixmap(12,8);
				pixmap.fill(QColor(color.GetRed()*255,color.GetGreen()*255,color.GetBlue()*255));
				actorItem->setPixmap(0,pixmap);
				actorItem->setEnabled(true);
				actorMap.insert(std::make_pair(actorItem,actor));

			}

			m_btnManuLiverSegment->setEnabled(false);
		}
	}
	else
	{
		QMessageBox::warning ( NULL,
			tr("Sorry,No label result image to Surface."),tr(""),
			QMessageBox::Ok,  QMessageBox::NoButton,  QMessageBox::NoButton );
	}

	
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkAnnotation::SeparateLabelImage(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	typedef itk::Image<TPixel, VImageDimension> TImageType;
	typedef typename TImageType::RegionType TRegionType;
	typedef TImageType::SizeType Size;
	typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	typedef itk::ImageRegionIterator<TImageType> RegionIteratorType;
	
	for (int i=2;i<=vesselGraph->subGraphCount;i++)
	{
		TPixel pv = TPixel(i*10);
		DuplicatorType::Pointer duplicator = DuplicatorType::New();
		duplicator->SetInputImage(itkImage);
		duplicator->Update();
		typename TImageType::Pointer clonedImage = duplicator->GetOutput();
		RegionIteratorType iterator1(itkImage, itkImage->GetLargestPossibleRegion());
		RegionIteratorType iterator2(clonedImage, clonedImage->GetLargestPossibleRegion());
		for (iterator1.GoToBegin(),iterator2.Begin(); !iterator1.IsAtEnd(); ++iterator1,++iterator2)
		{
			if(iterator1.Get()==pv)
			{
				iterator2.Set(pv);		 
			}
			else
			{
				iterator2.Set(-10000);		
			}
		}

		mitk::DataTreeNode::Pointer node  = mitk::DataTreeNode::New();
		node->SetData(mitk::ImportItkImage(clonedImage));
		QString s;
		s.sprintf("Segment %d",i-1);
		node->SetName(s.ascii());
		mitk::DataStorage::GetInstance()->Add( node );
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();

		InitDataTree();
	}
	
}

void QmitkAnnotation::SetActorColor()
{
	assert(this->actor!=NULL);
	QPalette palette = m_btnActorColor->palette();
	const QColor & color =
		QColorDialog::getColor();
	if(color.isValid())
	{
		std::cout << color.red() << " " << color.green() << " " << color.blue() << std::endl;
		palette.setColor(QColorGroup::Button,color);
		palette.setBrush(QColorGroup::Background,QBrush(color));
		m_btnActorColor->setPalette(palette);
	}
	this->actor->GetProperty()->SetColor(color.red()/255.0,color.green()/255.0,color.blue()/255.0);

	QListViewItemIterator iter(listView);
	iter++;
	while(iter.current())
	{
		QCheckListItem *item = dynamic_cast<QCheckListItem *>(iter.current());
		if(item && item->isSelected())
		{
			QPixmap pixmap(12,8);
			pixmap.fill(QColor(color.red(),color.green(),color.blue()));
			item->setPixmap(0,pixmap);
			m_btnActorColor->setPixmap(pixmap);
			break;
		}
		iter++;
	}
}

void QmitkAnnotation::activeActorProperty(bool isActive)
{
	actorVisible->setEnabled(isActive);
	m_btnActorColor->setEnabled(isActive);
	opacity->setEnabled(isActive);

	if(isActive)
	{
		
		if (this->actor->GetVisibility())
		{
			actorVisible->setChecked(true);
		}
		else
		{
			actorVisible->setChecked(false);	
		}

		std::cout << "ACTOR OPACITY :" <<this->actor->GetProperty()->GetOpacity()<< std::endl;
		opacity->setValue(this->actor->GetProperty()->GetOpacity()*100);
		double color[3];
		this->actor->GetProperty()->GetColor(color);
		QPixmap pixmap(12,8);
		pixmap.fill(QColor(color[0]*255,color[1]*255,color[2]*255));
		m_btnActorColor->setPixmap(pixmap);
	}

	
}

void QmitkAnnotation::EnableVisible()
{
	assert(this->actor!=NULL);
	if(actorVisible->isChecked())
	{
		this->actor->VisibilityOn();
		
	}
	else
	{
		this->actor->VisibilityOff();
	}
	m_pRenWin->Render();
}

void QmitkAnnotation::DisplaySubtreeNodes()
{
	assert(vesselGraph!=NULL && vesselGraph->subGraphCount>1);
	std::list<Vertex>::iterator v_i;
	mitk::Point3D point3D;
	for (int i = 2;i<=vesselGraph->subGraphCount;i++)
	{
		mitk::PointSet::Pointer points = mitk::PointSet::New();
		int index = 0;
		for (v_i=vesselGraph->QDiVertex[i].begin();v_i!=vesselGraph->QDiVertex[i].end();++v_i)
		{
			point3D[0] = vesselGraph->txMap[*v_i];
			point3D[1] = vesselGraph->tyMap[*v_i];
			point3D[2] = vesselGraph->tzMap[*v_i];
			points->GetPointSet()->GetPoints()->InsertElement( index++, point3D );
		}
		mitk::DataTreeNode::Pointer pNode = mitk::DataTreeNode::New();
		pNode->SetData( points );
		char name[20];
		sprintf(name,"No.%d SubGraph",i);
		pNode->SetProperty("name", mitk::StringProperty::New(name));
		pNode->SetProperty("layer", mitk::IntProperty::New(1));
		if (vesselGraph->subGraphCount%8==0)
		{
			m_RainbowColor.GoToBegin();
		}
		mitk::Color color = m_RainbowColor.GetNextColor();
		pNode->SetProperty("color", mitk::ColorProperty::New(color.GetRed(),color.GetGreen(),color.GetBlue()));
		pNode->SetProperty("pointsize", mitk::FloatProperty::New(1));
		mitk::DataStorage::GetInstance()->Add( pNode );
		mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	}
}