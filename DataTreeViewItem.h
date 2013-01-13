/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date: 2009-04-23 19:50:34 +0800 (星期�? 23 四月 2009) $
Version:   $Revision: 16947 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#ifndef DATA_TREE_VIEW_H
#define DATA_TREE_VIEW_H
#include <qlistview.h>
#include <qstring.h>

class QWidget;
#include <mitkDataTree.h>

/*
class DataTreeViewItem : public QListViewItem
{

public:
	DataTreeViewItem( QListView *parent, const QString &s1, const QString &s2, mitk::DataTreeIteratorBase * nodeIt, mitk::DataTreeNode* SelectedNode);
	DataTreeViewItem( DataTreeViewItem * parent, mitk::DataTreeIteratorBase * nodeIt );
	virtual ~DataTreeViewItem() {}

	mitk::DataTreeNode::Pointer GetDataTreeNode() const {return m_TreeNode;}
	mitk::DataTreeIteratorBase* GetDataTreeIterator() {return m_DataTreeIterator.GetPointer();}

protected:
	mitk::DataTreeNode::Pointer m_TreeNode;
	mitk::DataTreeIteratorClone m_DataTreeIterator;
	mitk::DataTreeNode* SelectedNode;
};
*/


class DataTreeViewItem : public QObject, public QCheckListItem
{
	Q_OBJECT
public:
	DataTreeViewItem( QListView *parent, const QString &s1, const QString &s2, mitk::DataTreeIteratorBase * nodeIt, mitk::DataTreeNode* SelectedNode);
	DataTreeViewItem( DataTreeViewItem * parent, mitk::DataTreeIteratorBase * nodeIt );
	virtual ~DataTreeViewItem() {}

	mitk::DataTreeNode::Pointer GetDataTreeNode() const {return m_TreeNode;}
	mitk::DataTreeIteratorBase* GetDataTreeIterator() {return m_DataTreeIterator.GetPointer();}

public slots:
	void stateChange (bool);
//signals:
//	void statusChange (DataTreeViewItem*, bool);


protected:
	mitk::DataTreeNode::Pointer m_TreeNode;
	mitk::DataTreeIteratorClone m_DataTreeIterator;
	mitk::DataTreeNode* SelectedNode;
};




#endif //QMITK_DATA_TREE_VIEW_H
