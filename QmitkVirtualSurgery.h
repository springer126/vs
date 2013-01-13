/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date: 2008-06-10 23:03:03 +0800 (星期�? 10 六月 2008) $
Version:   $Revision: 14578 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#if !defined(QmitkVirtualSurgery_H__INCLUDED)
#define QmitkVirtualSurgery_H__INCLUDED

#include "QmitkFunctionality.h"
#include "mitkTestingConfig.h"
/** add 2012-09-20 **/
#include <mitkImage.h>
#include <itkImage.h>

#include "mitkPointSet.h"
#include "mitkPointSetInteractor.h"
#include "VesselGraph.h"
class QmitkStdMultiWidget;
class QmitkVirtualSurgeryControls;

/*!
\brief QmitkVirtualSurgery 

One needs to reimplement the methods CreateControlWidget(..), CreateMainWidget(..) 
and CreateAction(..) from QmitkFunctionality. 

\sa QmitkFunctionality
\ingroup Functionalities
*/
class QmitkVirtualSurgery : public QmitkFunctionality
{  
  Q_OBJECT
  
  public:  
  /*!  
  \brief default constructor  
  */  
  QmitkVirtualSurgery(QObject *parent=0, const char *name=0, QmitkStdMultiWidget *mitkStdMultiWidget = NULL, mitk::DataTreeIteratorBase* dataIt = NULL);

  /*!  
  \brief default destructor  
  */  
  virtual ~QmitkVirtualSurgery();

  /*!  
  \brief method for creating the widget containing the application controls, like sliders, buttons etc.  
  */  
  virtual QWidget * CreateControlWidget(QWidget *parent);

  /*!  
  \brief method for creating the applications main widget  
  */  
  virtual QWidget * CreateMainWidget(QWidget * parent);

  /*!  
  \brief method for creating the connections of main and control widget  
  */  
  virtual void CreateConnections();

  /*!  
  \brief method for creating an QAction object, i.e. button & menu entry  @param parent the parent QWidget  
  */  
  virtual QAction * CreateAction(QActionGroup *parent);

  virtual void Activated();


private:
	template <typename TPixel, unsigned int VImageDimension>
	void ErodeVessel(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Geometry3D* imageGeometry);

	template <typename TPixel, unsigned int VImageDimension>
	void DoPortalVeinSegment(itk::Image<TPixel, VImageDimension> *itkImage,mitk::Geometry3D* imageGeometry);

	template <typename TPixel, unsigned int VImageDimension>
	void MyCastPixelType(itk::Image<TPixel, VImageDimension> *itkImage, mitk::Image::Pointer *pointer);

	template <typename TPixel, unsigned int VImageDimension>
	void GetMaxMinPixelIntensity(itk::Image<TPixel, VImageDimension> *itkImage, int &maxIntensity, int &minIntensity);

#ifdef BUILD_TESTING
  /**
  \brief This method performs an automated functionality test. 
  
  Will be called by testing subsystem. Do not call it yourself!
  The method should be implemented in its own cpp file 'QmitkMyFunctionalityTesting.cpp'.
  */
  virtual int TestYourself();
#endif

protected slots:  
  void TreeChanged();
  /*
   * just an example slot for the example TreeNodeSelector widget
   */
  void ImageSelected(mitk::DataTreeIteratorClone imageIt);
  
  /*
   * just an example slot for the example "Start!" button
   */
   void StartButtonClicked();
   void PortalVeinSegmentButtonClicked();
   void ManualErodeVesselClicked();


protected:  
  /*!  
  * default main widget containing 4 windows showing 3   
  * orthogonal slices of the volume and a 3d render window  
  */  
  QmitkStdMultiWidget * m_MultiWidget;

  /*!  
  * controls for the functionality. Sliders, Buttons, TreenodeSelectors,...  
  */  
  QmitkVirtualSurgeryControls * m_Controls;

  mitk::DataTreeIteratorClone selectedImage;
  mitk::PixelType originalPixelType;

  //brief This node is created once and used for storing seed points
  mitk::DataTreeNode::Pointer m_PointSetNode;

  //brief This is the actual seed point data object
  mitk::PointSet::Pointer m_PointSet;

  mitk::PointSetInteractor::Pointer m_Interactor;
  VesselGraph* vesselGraph;
};
#endif // !defined(QmitkVirtualSurgery_H__INCLUDED)
