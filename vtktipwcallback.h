
#pragma once
#include "vtkcommand.h"
#include <vtkPlane.h>
#include <vtkImageData.h>
#include <vtkImplicitPlaneWidget.h>

class vtkTIPWCallback :
	public vtkCommand
{
public:
	vtkPlane *plane;
	vtkImageData *data;

public:
	//static vtkTIPWCallback *New()
	//{
	//	return new vtkTIPWCallback;
	//}

	vtkTIPWCallback(void)
	{

	}

	virtual void Execute(vtkObject *caller,unsigned long,void *)
	{
		vtkImplicitPlaneWidget *planeWidget = reinterpret_cast<vtkImplicitPlaneWidget *>(caller);
		planeWidget->GetPlane(this->plane);
		double normal[3] = {plane->GetNormal()[0],plane->GetNormal()[1],plane->GetNormal()[2]};
		double origin[3] = {plane->GetOrigin()[0],plane->GetOrigin()[1],plane->GetOrigin()[2]};


		std::cout << "Normal:" ;
		for (int i=0;i<3;i++)
		{
			std::cout << normal[i] <<" ";
		}
		std::cout << endl;
		std::cout << "Origin:" ;
		for (int i=0;i<3;i++)
		{
			std::cout << origin[i] <<" ";
		}
		std::cout << endl;
	}

public:
	~vtkTIPWCallback(void)
	{
		plane->Delete();
		data->Delete();
	}
};
