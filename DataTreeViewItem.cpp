/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date: 2008-09-04 19:54:36 +0800 (鏄熸湡鍥? 04 涔濇湀 2008) $
Version:   $Revision: 15161 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#include <sstream>
#include "DataTreeViewItem.h"
#include "mitkPropertyList.h"
#include "mitkStringProperty.h"

/*
DataTreeViewItem::DataTreeViewItem( QListView *parent, const QString &s1 , const QString &s2 , mitk::DataTreeIteratorBase* nodeIt, mitk::DataTreeNode* SelectedNode)
: QListViewItem( parent, s1, s2 ), m_DataTreeIterator(NULL)
{
	this->SelectedNode = SelectedNode;
	std::cout << "SelectedNode = " << SelectedNode << std::endl;
	if(SelectedNode)
	{
		assert(nodeIt!=NULL);
		m_DataTreeIterator = nodeIt;

		m_TreeNode = nodeIt->Get();
		QListViewItem::setText(0,QString("All Datasets"));

		mitk::DataTreeChildIterator it(*nodeIt);
		while (!it.IsAtEnd())
		{
			//if(SelectedNode == m_TreeNode)
				new DataTreeViewItem(this,&it);
			++it;
		}
		setOpen(true);
	}
}
*/


DataTreeViewItem::DataTreeViewItem( QListView *parent, const QString &s1 , const QString &s2 , mitk::DataTreeIteratorBase* nodeIt, mitk::DataTreeNode* SelectedNode)
: QCheckListItem( parent, s1, QCheckListItem::CheckBox), m_DataTreeIterator(NULL)
{
	this->SelectedNode = SelectedNode;
	std::cout << "SelectedNode = " << SelectedNode << std::endl;
	if(SelectedNode)
	{
		assert(nodeIt!=NULL);
		m_DataTreeIterator = nodeIt;

		m_TreeNode = nodeIt->Get();
		QListViewItem::setText(0,QString("All Datasets"));

		mitk::DataTreeChildIterator it(*nodeIt);
		while (!it.IsAtEnd())
		{
			//if(SelectedNode == m_TreeNode)
			new DataTreeViewItem(this,&it);
			++it;
		}
		setOpen(true);
	}
}



DataTreeViewItem::DataTreeViewItem( DataTreeViewItem * parent, mitk::DataTreeIteratorBase * nodeIt )
: QCheckListItem(parent, QString(""), QCheckListItem::CheckBox ), m_DataTreeIterator(NULL) 
{
	assert(nodeIt!=NULL);
	m_DataTreeIterator = nodeIt;

	m_TreeNode = nodeIt->Get();
	//char name[256];
	std::string name;

	if(m_TreeNode.IsNotNull())
	{
		if (m_TreeNode->GetName(name, NULL)) {
			std::cout << "name = " << name << std::endl;
			//if(name=="widget1Plane" || name=="widget2Plane" || name=="widget3Plane")
			//	;
			//else
			//{
			//std::cout << name << std::endl;
			QListViewItem::setText(0, QString(name.c_str()));
			//}
		} else {
			QListViewItem::setText(0,"No name!"); 
		}

		std::stringstream ss;
		ss << m_TreeNode->GetReferenceCount() << "/";
		mitk::BaseData::ConstPointer bd = m_TreeNode->GetData();
		if (bd.IsNotNull()) {
			//std::cout << bd->GetNameOfClass() << std::endl;
			QListViewItem::setText(1,QString(bd->GetNameOfClass())) ;
			ss << m_TreeNode->GetData()->GetReferenceCount();
		} else {
			QListViewItem::setText(1,QString("empty DataTreeNode"));
			ss << "-";
		}
		QListViewItem::setText(2,QString(ss.str().c_str()));
	}
	else
	{
		QListViewItem::setText(0,"NULL node!");
	}

	mitk::DataTreeChildIterator it(*nodeIt);
	while (!it.IsAtEnd())
	{
		new DataTreeViewItem(this,&it);
		++it;
	}
	setOpen(true);
}


//如果选了某一个结点，则同样选中它的所有子结点
void DataTreeViewItem::stateChange(bool value)
{
	//if(value)
	//	std::cout << "true" <<std::endl;
	//else
	//	std::cout << "false" <<std::endl;

	//得到第一个孩子
	QListViewItem *child = firstChild();
	int count = childCount();
	//std::cout << "childCount = " << count << std::endl;


	QListViewItem *p = child;
	for (int i=0; p&&i<count; i++)
	{
		//std::cout << "value = " << value << std::endl;
		DataTreeViewItem *item = dynamic_cast<DataTreeViewItem*>(p);
		//std::cout << item->text(1).ascii() << std::endl;
		//std::cout << "item = " << item << std::endl;
		p = p->nextSibling();
		//statusChange(item, value);
		if(item)
		{
			item->setOn(value);
		}
	}
}

/*
DataTreeViewItem::DataTreeViewItem( DataTreeViewItem * parent, mitk::DataTreeIteratorBase * nodeIt )
: QListViewItem(parent),m_DataTreeIterator(NULL) 
{
	assert(nodeIt!=NULL);
	m_DataTreeIterator = nodeIt;

	m_TreeNode = nodeIt->Get();
	//char name[256];
	std::string name;

	if(m_TreeNode.IsNotNull())
	{
		if (m_TreeNode->GetName(name, NULL)) {
			//if(name=="widget1Plane" || name=="widget2Plane" || name=="widget3Plane")
			//	;
			//else
			//{
				std::cout << name << std::endl;
				QListViewItem::setText(0, QString(name.c_str()));
			//}
		} else {
			QListViewItem::setText(0,"No name!");
		}

		std::stringstream ss;
		ss << m_TreeNode->GetReferenceCount() << "/";
		mitk::BaseData::ConstPointer bd = m_TreeNode->GetData();
		if (bd.IsNotNull()) {
			std::cout << bd->GetNameOfClass() << std::endl;
			QListViewItem::setText(1,QString(bd->GetNameOfClass())) ;
			ss << m_TreeNode->GetData()->GetReferenceCount();
		} else {
			QListViewItem::setText(1,QString("empty DataTreeNode"));
			ss << "-";
		}
		QListViewItem::setText(2,QString(ss.str().c_str()));
	}
	else
	{
		QListViewItem::setText(0,"NULL node!");
	}

	mitk::DataTreeChildIterator it(*nodeIt);
	while (!it.IsAtEnd())
	{
		new DataTreeViewItem(this,&it);
		++it;
	}
	setOpen(true);
}

*/






/*
DataTreeViewItem::DataTreeViewItem( DataTreeViewItem * parent, mitk::DataTreeIteratorBase * nodeIt )
: QListViewItem(parent),m_DataTreeIterator(NULL) 
{
	assert(nodeIt!=NULL);
	m_DataTreeIterator = nodeIt;

	m_TreeNode = nodeIt->Get();
	//char name[256];
	std::string name;

	if(m_TreeNode.IsNotNull())
	{
		if (m_TreeNode->GetName(name, NULL)) {
			if(name=="" || name=="Widgets" || name=="widget1Plane" || name=="widget2Plane" || name=="widget3Plane")
				return;
			else
			{
				std::cout << name << std::endl;
				QListViewItem::setText(0, QString(name.c_str()));
			}
		}
// 		else {
// 			QListViewItem::setText(0,"No name!");
// 		}

// 		std::stringstream ss;
// 		ss << m_TreeNode->GetReferenceCount() << "/";
// 		mitk::BaseData::ConstPointer bd = m_TreeNode->GetData();
// 		if (bd.IsNotNull()) {
// 			std::cout << bd->GetNameOfClass() << std::endl;
// 			QListViewItem::setText(1,QString(bd->GetNameOfClass())) ;
// 			ss << m_TreeNode->GetData()->GetReferenceCount();
// 		} else {
// 			QListViewItem::setText(1,QString("empty DataTreeNode"));
// 			ss << "-";
// 		}
// 		QListViewItem::setText(2,QString(ss.str().c_str()));

		mitk::DataTreeChildIterator it(*nodeIt);
		while (!it.IsAtEnd())
		{
			new DataTreeViewItem(this,&it);
			++it;
		}
		//setOpen(true);
	}
	else
	{
		QListViewItem::setText(0,"NULL node!");
	}

// 	mitk::DataTreeChildIterator it(*nodeIt);
// 	while (!it.IsAtEnd())
// 	{
// 		new DataTreeViewItem(this,&it);
// 		++it;
// 	}

//	setOpen(true);
}

*/
