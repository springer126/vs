/*=========================================================================
 
Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date: 2008-09-04 19:15:46 +0800 (鏄熸湡鍥? 04 涔濇湀 2008) $
Version:   $Revision: 15159 $
 
Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.
 
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.
 
=========================================================================*/

#include "QmitkVirtualSurgery.h"
#include "QmitkVirtualSurgeryControls.h"
#include <qaction.h>
#include "icon.xpm"
#include "QmitkTreeNodeSelector.h"
#include "QmitkStdMultiWidget.h"
#include "mitkStatusBar.h"
#include "mitkProgressBar.h"
#include "QmitkAnnotation.h"
#include "qmessagebox.h"
#include "mitkTransferFunction.h"
#include "mitkColorSequenceRainbow.h"

/***add 2012-09-20**/
#include "mitkImageAccessByItk.h"
#include <iostream>
#include "itkBinaryBallStructuringElement.h"
#include <itkImageRegionIterator.h>
#include <itkImageDuplicator.h>
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "mitkDataStorage.h"
#include "mitkStringProperty.h"
#include "mitkDataTreeNode.h"
#include "mitkProperties.h"
#include "mitkITKImageImport.h"
//**add 2012-09-21**//
#include <vtkImageGaussianSmooth.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyData.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageData.h>
#include <vtkDataObject.h>
#include <itkImage.h>
#include <itkImageDuplicator.h>
#include <itkCurvatureFlowImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkTranslationTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <string>

#include <vtkImageChangeInformation.h>
#include <vtkLinearTransform.h>
#include <vtkMarchingContourFilter.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>


//**add 2012-09-25**//
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>

//add 2012-10-18
#include "mitkGlobalInteraction.h"
#include <itkConnectedThresholdImageFilter.h>
#include <itkCommand.h>
#include <limits>

#include "VesselGraph.h"


//add 2012-11-20
#include <mitkLabeledImageToSurfaceFilter.h>

QmitkVirtualSurgery::QmitkVirtualSurgery(QObject *parent, const char *name, QmitkStdMultiWidget *mitkStdMultiWidget, mitk::DataTreeIteratorBase* it)
    : QmitkFunctionality(parent, name, it), m_MultiWidget(mitkStdMultiWidget), m_Controls(NULL)
{
  SetAvailability(true);
  vesselGraph = NULL;
}


QmitkVirtualSurgery::~QmitkVirtualSurgery()
{
	if(vesselGraph!=NULL)
		delete vesselGraph;
}


QWidget * QmitkVirtualSurgery::CreateMainWidget(QWidget *parent)
{
  if ( m_MultiWidget == NULL )
  {
    m_MultiWidget = new QmitkStdMultiWidget( parent );
  }
  return m_MultiWidget;
}


QWidget * QmitkVirtualSurgery::CreateControlWidget(QWidget *parent)
{
  if (m_Controls == NULL)
  {
    m_Controls = new QmitkVirtualSurgeryControls(parent);
  }
  return m_Controls;
}


void QmitkVirtualSurgery::CreateConnections()
{
  if ( m_Controls )
  {
    connect( (QObject*)(m_Controls->m_TreeNodeSelector), SIGNAL(Activated(mitk::DataTreeIteratorClone)),(QObject*) this, SLOT(ImageSelected(mitk::DataTreeIteratorClone)) );
    connect( (QObject*)(m_Controls->m_StartButton), SIGNAL(clicked()),(QObject*) this, SLOT(StartButtonClicked()));
	connect( (QObject*)(m_Controls->m_btnManualErode), SIGNAL(clicked()),(QObject*) this, SLOT(ManualErodeVesselClicked()));
	connect( (QObject*)(m_Controls->m_btnPortalVeinSegment), SIGNAL(clicked()),(QObject*) this, SLOT(PortalVeinSegmentButtonClicked()));
	
  }
}


QAction * QmitkVirtualSurgery::CreateAction(QActionGroup *parent)
{
  QAction* action;
  action = new QAction(
	  tr( "Virtual Surgery" ), 
	  QPixmap((const char**)icon_xpm), 
	  tr( "QmitkVirtualSurgery menu" ), 
	  0, 
	  parent, 
	  "QmitkVirtualSurgery" );
  return action;
}


void QmitkVirtualSurgery::TreeChanged()
{
   m_Controls->m_TreeNodeSelector->SetDataTreeNodeIterator(this->GetDataTreeIterator());
}


void QmitkVirtualSurgery::Activated()
{
  QmitkFunctionality::Activated();
  if (m_PointSetNode.IsNull()) // only once create a new DataTreeNode containing a PointSet with some interaction
  {
	  // new node and data item
	  m_PointSetNode = mitk::DataTreeNode::New();
	  m_PointSetNode->SetName("Seedpoints for region growing");
	  m_PointSet = mitk::PointSet::New();
	  m_PointSetNode->SetData( m_PointSet );
	  m_Interactor = mitk::PointSetInteractor::New("pointsetinteractor", m_PointSetNode);
	  // add the pointset to the data tree (for rendering)
	  GetDataTreeIterator()->Add( m_PointSetNode );
  }
  // new behavior/interaction for the pointset node

  mitk::GlobalInteraction::GetInstance()->AddInteractor( m_Interactor );
}


void QmitkVirtualSurgery::ImageSelected(mitk::DataTreeIteratorClone imageIt)
{
  assert( imageIt.IsNotNull() ); // should never fail, the selection widget cares for that
  std::cout << " Image Selected " << std::endl;
  mitk::DataTreeNode* node = imageIt->Get();
  selectedImage = imageIt;
  if ( node )
  {
    // here we have a valid mitk::DataTreeNode
    std::string name;
    if (node->GetName(name))
    {
      // a property called "name" was found for this DataTreeNode
      std::cout << "Tree node selected with name '" << name << "'" << std::endl;
    }

    // a node itself is not very useful, we need its data item
    mitk::BaseData* data = node->GetData();
    if (data)
    {
      // test if this data item is an image or not (could also be a surface or something totally different)
      mitk::Image* image = dynamic_cast<mitk::Image*>( data );
      if (image)
      {
        std::cout << "Surprise: this node contains a real image dataset." << std::endl;
		originalPixelType = image->GetPixelType();
      }
    }
  }
}

void QmitkVirtualSurgery::StartButtonClicked() 
{
	if(selectedImage.IsNotNull())
	{
		QmitkAnnotation *AnnotationWindow = QmitkAnnotation::Instance(this, this->GetDataTreeIterator(), selectedImage->Get() );
	}

	mitk::TransferFunction::Pointer tf = mitk::TransferFunction::New();
	//tf->SetMyTransferFunction();
    std::cout << "Start Button clicked!" << std::endl;
    WaitCursorOn(); // always good to show the user that the application is processing and will not react to user input for a while
    mitk::StatusBar::GetInstance()->DisplayText("QmitkVirtualSurgery is doing something...", 4000);  // tell the user what you are doing
    mitk::ProgressBar::GetInstance()->AddStepsToDo(2);  // use progress bar to show that the application is doing something
  
    // Do something here, add progress while you do it:
  
    mitk::ProgressBar::GetInstance()->Progress();
    mitk::ProgressBar::GetInstance()->Progress(); // last step, we're finished!
    mitk::StatusBar::GetInstance()->DisplayText("QmitkVirtualSurgery has finished...", 2000);
    WaitCursorOff();  // restore normal mouse cursor after you finished
}


void QmitkVirtualSurgery::PortalVeinSegmentButtonClicked()
{
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
				std::cout << "Erode process begin..." << std::endl;
				//AccessByItk_2(image, ErodeVessel, &erodeImage,image->GetGeometry());
				AccessByItk_1(image,DoPortalVeinSegment,image->GetGeometry());
				std::cout << "segment process end..." << std::endl;
			}
		}
	}
}


void QmitkVirtualSurgery::ManualErodeVesselClicked()
{
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
				std::cout << "Erode process begin..." << std::endl;
				AccessByItk_1(image, ErodeVessel,image->GetGeometry());
				std::cout << "segment process end..." << std::endl;
			}
		}
	}
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkVirtualSurgery::ErodeVessel(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Geometry3D* imageGeometry)
{
	std::cout << "Function ErodeVessel begin..." <<std::endl;
	typedef itk::Image<TPixel, VImageDimension> TImageType;
	typedef TImageType::IndexType  IndexType;
	typedef TImageType::SizeType   SizeType;

	typedef itk::ImageDuplicator< typename TImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage( itkImage );
	duplicator->Update();
	typename TImageType::Pointer clonedImage = duplicator->GetOutput();
	mitk::Image::Pointer resultImage = mitk::ImportItkImage( clonedImage );
	mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
	newNode->SetData(resultImage);
	newNode->SetProperty("name", mitk::StringProperty::New("Erode Vessel Image"));
	newNode->SetProperty("opacity", mitk::FloatProperty::New(0.0));
	mitk::DataStorage::GetInstance()->Add( newNode );
	mitk::RenderingManager::GetInstance()->RequestUpdateAll();

	//typedef itk::ImageRegionIterator<TImageType> RegionIteratorType;
	//typedef itk::ImageRegionIterator<TImageType> ImageIteratorType;
	//typedef itk::BinaryBallStructuringElement<float,VImageDimension> StructuringElementType;
	//typedef itk::BinaryDilateImageFilter<TImageType,
	//	TImageType, 
	//	StructuringElementType > DilateFilterType;
	//typedef itk::BinaryErodeImageFilter<TImageType, 
	//	TImageType, 
	//	StructuringElementType > ErodeFilterType;
	

	/*
	StructuringElementType structuringElement;
	structuringElement.SetRadius(1);
	structuringElement.CreateStructuringElement();

	ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();
	DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
	binaryErode->SetKernel( structuringElement );
	binaryDilate->SetKernel( structuringElement );
	binaryErode->SetErodeValue( 1 );
	binaryErode->SetBackgroundValue(0);
	binaryDilate->SetDilateValue( 1 );
	binaryDilate->SetBackgroundValue(0);

	//闭操作：先膨胀再腐蚀
	binaryDilate->SetInput(itkImage);
	binaryErode->SetInput(binaryDilate->GetOutput());
	binaryErode->Update();
	//Update，完成闭操作

	//只腐蚀
	binaryErode->SetInput(itkImage);
	binaryDilate->SetInput(binaryErode->GetOutput());
	binaryDilate->Update();

	TImageType::Pointer outputImage = binaryDilate->GetOutput();
	*pointer = mitk::ImportItkImage(outputImage);
	mitk::Image::Pointer resultImage = mitk::ImportItkImage(outputImage);

	mitk::DataTreeNode::Pointer node = mitk::DataTreeNode::New();
	node->SetData(resultImage);
	node->SetProperty("name", mitk::StringProperty::New("Result Vessel"));
	node->SetProperty("opacity", mitk::FloatProperty::New(1.0));
	node->SetProperty("Surface", mitk::BoolProperty::New(true));
	mitk::DataStorage::GetInstance()->Add( node );
	*/

	//add 2012-10-15
	/*
	typedef itk::Image<TPixel, 2>  SliceImageType;
	typedef SliceImageType::SizeType   SliceSizeType;
	typedef itk::ImageRegionIterator<SliceImageType> SliceIteratorType;
	SizeType ImageSize = itkImage->GetLargestPossibleRegion().GetSize();
	SliceSizeType  SliceSize;
	SliceSize[0] = ImageSize[0];
	SliceSize[1] = ImageSize[1];

	SliceImageType::Pointer slicer = SliceImageType::New();
	slicer->SetRegions(SliceSize);
	slicer->Allocate();
	SliceImageType::RegionType sliceRegion = slicer->GetLargestPossibleRegion();
	
	TImageType::SizeType   RequestSize;
	TImageType::IndexType  RequestIndex;
	TImageType::RegionType RequestRegion;

	RequestSize[0] = ImageSize[0];
	RequestSize[1] = ImageSize[1];
	RequestSize[2] = 1;

	RequestIndex[0] = 0;
	RequestIndex[1] = 0;

	RequestRegion.SetSize(RequestSize);

	typedef itk::BinaryBallStructuringElement<float,2> StructuringElementType2;
	typedef itk::BinaryDilateImageFilter<SliceImageType,SliceImageType, StructuringElementType2 > DilateFilterType2;
	typedef itk::BinaryErodeImageFilter<SliceImageType, SliceImageType, StructuringElementType2 > ErodeFilterType2;
	StructuringElementType2 structuringElement2;
	structuringElement2.SetRadius(2);
	structuringElement2.CreateStructuringElement();
	
	int slices = ImageSize[2];
	for (int i=0; i<slices; i++)
	{
		RequestIndex[2] = i;
		RequestRegion.SetIndex(RequestIndex);
		
		RegionIteratorType it(itkImage, RequestRegion);
		SliceIteratorType  sliceIt(slicer, sliceRegion);

		for(it.GoToBegin(), sliceIt.GoToBegin(); !it.IsAtEnd(); ++it, ++sliceIt)
		{
			sliceIt.Set(it.Get() );
		}

		ErodeFilterType2::Pointer binaryErode1 = ErodeFilterType2::New();
		DilateFilterType2::Pointer binaryDilate1 = DilateFilterType2::New();
		binaryErode1->SetKernel( structuringElement2 );
		binaryDilate1->SetKernel( structuringElement2 );
		binaryErode1->SetErodeValue( 1 );
		binaryDilate1->SetDilateValue( 1 );

		ErodeFilterType2::Pointer binaryErode2 = ErodeFilterType2::New();
		DilateFilterType2::Pointer binaryDilate2 = DilateFilterType2::New();
		binaryErode2->SetKernel( structuringElement2 );
		binaryDilate2->SetKernel( structuringElement2 );
		binaryErode2->SetErodeValue( 1 );
		binaryDilate2->SetDilateValue( 1 );

		//1.开操作：先腐蚀后膨胀
		binaryErode1->SetInput(slicer);
		binaryDilate1->SetInput(binaryErode1->GetOutput());
		binaryDilate1->Update();   //先Update，完成第一步开操作
		//2.闭操作：先膨胀再腐蚀
		//binaryDilate2->SetInput( binaryDilate1->GetOutput() );
		//binaryErode2->SetInput(binaryDilate2->GetOutput());
		//binaryErode2->Update();    //Update，完成闭操作

		SliceImageType::Pointer mask = binaryDilate1->GetOutput();

		SliceIteratorType it3(mask, sliceRegion);
		for(it.GoToBegin(), it3.GoToBegin(); !it.IsAtEnd(); ++it, ++it3)
		{
			it.Set(it3.Get());
		}
		
		
	}
	*/
	
	IndexType seedIndex;
	for ( mitk::PointSet::PointsConstIterator pointsIterator = m_PointSet->GetPointSet()->GetPoints()->Begin(); // really nice syntax to get an interator for all points
		 pointsIterator != m_PointSet->GetPointSet()->GetPoints()->End();
		 ++pointsIterator ) 
	{
		// first test if this point is inside the image at all
		if ( !imageGeometry->IsInside( pointsIterator.Value()) ) 
			continue;

		// convert world coordinates to image indices
		imageGeometry->WorldToIndex( pointsIterator.Value(), seedIndex);
		for (int i=-1;i<2;i++)
		{
			for (int j=-1;j<2;j++)
			{
				for (int k=-1;k<2;k++)
				{
					IndexType index;
					index[0]=seedIndex[0]+i;
					index[1]=seedIndex[1]+j;
					index[2]=seedIndex[2]+k;
					clonedImage->SetPixel(index,0);
				}
				
			}
		}
	}	

	//*pointer = mitk::ImportItkImage(processImage);
	//mitk::Image::Pointer resultImage = mitk::ImportItkImage(processImage);
	//mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
	//newNode->SetData(resultImage);
	//newNode->SetProperty("name", mitk::StringProperty::New("erode vessel image"));
	//newNode->SetProperty("opacity", mitk::FloatProperty::New(0.5));
	//mitk::DataStorage::GetInstance()->Add( newNode );
	//mitk::RenderingManager::GetInstance()->RequestUpdateAll();
	std::cout << "Function ErodeVessel end..." <<std::endl;
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkVirtualSurgery::DoPortalVeinSegment(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Geometry3D* imageGeometry)
{
	typedef itk::Image< TPixel, VImageDimension > InputImageType;
	typedef typename InputImageType::IndexType    IndexType;
	
	// instantiate an ITK region growing filter, set its parameters
	typedef itk::ConnectedThresholdImageFilter<InputImageType, InputImageType> RegionGrowingFilterType;
	typename RegionGrowingFilterType::Pointer regionGrower = RegionGrowingFilterType::New();
	regionGrower->SetInput( itkImage ); // don't forget this

	// determine a thresholding interval
	IndexType seedIndex;
	TPixel min( std::numeric_limits<TPixel>::max() );
	TPixel max( std::numeric_limits<TPixel>::min() );
	for ( mitk::PointSet::PointsConstIterator pointsIterator = m_PointSet->GetPointSet()->GetPoints()->Begin(); // really nice syntax to get an interator for all points
		pointsIterator != m_PointSet->GetPointSet()->GetPoints()->End();
		++pointsIterator ) 
	{
		// first test if this point is inside the image at all
		if ( !imageGeometry->IsInside( pointsIterator.Value()) ) 
			continue;

		// convert world coordinates to image indices
		imageGeometry->WorldToIndex( pointsIterator.Value(), seedIndex);

		// get the pixel value at this point
		TPixel currentPixelValue = itkImage->GetPixel( seedIndex );

		// adjust minimum and maximum values
		if (currentPixelValue > max)
			max = currentPixelValue;

		if (currentPixelValue < min)
			min = currentPixelValue;

		regionGrower->AddSeed( seedIndex );
	}

	std::cout << "Values between " << min << " and " << max << std::endl;

	//min -= 30;
	//max += 30;

	// set thresholds and execute filter
	regionGrower->SetLower( min );
	regionGrower->SetUpper( max );

	regionGrower->Update();

	mitk::Image::Pointer resultImage = mitk::ImportItkImage( regionGrower->GetOutput() );
	mitk::DataTreeNode::Pointer newNode = mitk::DataTreeNode::New();
	newNode->SetData( resultImage );

	// set some properties
	newNode->SetProperty("binary", mitk::BoolProperty::New(true));
	newNode->SetProperty("name", mitk::StringProperty::New("Portal Vein Image"));
	newNode->SetProperty("color", mitk::ColorProperty::New(1.0,0.0,0.0));
	//newNode->SetProperty("volumerendering", mitk::BoolProperty::New(true));
	newNode->SetProperty("layer", mitk::IntProperty::New(1));
	newNode->SetProperty("opacity", mitk::FloatProperty::New(0.5));

	// add result to data tree
	mitk::DataStorage::GetInstance()->Add( newNode );

	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

template <typename TPixel, unsigned int VImageDimension>
void QmitkVirtualSurgery::MyCastPixelType(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer)
{
	std::cout<< "Function MyCastPixelType :Begin."<<std::endl;

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

	//originalPixelType = (dynamic_cast<mitk::Image *>(selectedImage->Get()->GetData()))->GetPixelType();
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


