/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#include "NewickNode.h"

BaseNode::BaseNode(BaseNode const &other)
{
  copy(other);
}


BaseNode const &BaseNode::operator=(BaseNode const &other)
{
  if (this != &other){
    destroy();
    copy(other);
  }

  return *this;
}


BaseNode::~BaseNode()
{
  destroy();
}


BaseNode *const BaseNode::makeSubTreeCopy() const
{
  return newCopy(this);
}


BaseNode *BaseNode::popChild()
{
  if (getChildrenNumber() == 0){
    cerr << "attempt to popChild from node with no children" << endl;
    exit (1);
  }

  BaseNode *np;
  np = getChild(getChildrenNumber() - 1);
  children.pop_back();

  return np;
}


void BaseNode::destroy()
{
  int n = getChildrenNumber();
  for(int i = 0; i < n; i++)
    delete children[i];
  children.clear();
}

  
BaseNode *const BaseNode::newCopy(BaseNode const *other) const
{
  return other->copyNode();
}


void BaseNode::copy(BaseNode const &other)
{
  parent = other.getParent();

  for(int i = 0; i < other.getChildrenNumber(); i++) {
    children.push_back(newCopy(other.children[i]) );
    children[i]->parent = this;
  }
}

