/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#ifndef __NewickNode_h
#define __NewickNode_h

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>

using namespace std;

/** 
 * Base class for nodes
 */
class BaseNode {
public:

  BaseNode(BaseNode *p) : parent(p) {}

  BaseNode(BaseNode const &other);

  BaseNode const &operator=(BaseNode const &other);

  virtual ~BaseNode();

  BaseNode *const makeSubTreeCopy() const;

  int getChildrenNumber() const {return children.size();};

  virtual BaseNode * const getChild(unsigned int n) const {return children[n];};

  virtual BaseNode * const getParent() const {return parent;};

  void addChild(BaseNode *c) {children.push_back(c);}

  BaseNode *popChild();

  void setParent(BaseNode *p) {parent = p;}

protected:
  //methods for handling construction, destruction etc.
  virtual void destroy();

  /** Return a pointer to a newly allocated copy of other (used by copy() ) */
  virtual BaseNode *const newCopy(BaseNode const *other) const;

  /** Overwite *this with a copy of other (copies all children recursively) */
  virtual void copy(BaseNode const &other);

  /** Must be defined in all derived classes in order to make copy
      funtion properly.  Should return Node of the relevant type. */
  virtual BaseNode * const copyNode() const {return new BaseNode(*this);}

  BaseNode *parent;
  vector<BaseNode *> children;
};

/** 
 * NewickNode class can hold the information specified in the
 * newick format.
 */
class NewickNode : public BaseNode{
public:
  NewickNode(BaseNode* p, string nm = "", double bl = 0) 
    : BaseNode(p),
      name(nm),
      branchLength(bl) {}
  virtual ~NewickNode() {};

  virtual NewickNode * const getChild(unsigned int n) const {return dynamic_cast<NewickNode *>(children[n]);}
  virtual NewickNode * const getParent() const {return dynamic_cast<NewickNode *>(parent);}
  string getName() const {return name;}
  double getBranchLength() const {return branchLength;}

  void setName(string n) {name = n;}
  void setBranchLength(double bl) {branchLength = bl;}

protected:
  /** Must be defined in all derived classes in order to make
      copy(const other) (inherited from BaseNode) funtion properly.
      Should return a pointer to a copy of *this. */
  BaseNode * const copyNode() const {  return new NewickNode(*this);}

  string name;
  /** Specify branch length to parent.  
   *  Set at -1 if undefined.   */
  double branchLength;
};

#endif  //__NewickNode_h
