%name NewickParserBase

%header{
#include <iostream>
#include <sstream>
#undef  yyFlexLexer
#define yyFlexLexer NewickFlexLexer
#include <FlexLexer.h>
#include "phy/NewickNode.h"
%}

%union {
  NewickNode *node;
  vector<NewickNode *> *nodePointerList;
  string *label;
  double value;
  };

%define LEX_BODY  {return lexer.yylex(); }
%define ERROR_BODY                                        \
{cerr << "syntax error in newick tree format before: \""  \
      << lexer.YYText()                                   \
      << "\"."                                            \
      << endl;}

%define MEMBERS                                           \
  virtual ~NewickParserBase() {}                          \
  NewickNode *parseTree(string const &tree);              \
  protected:                                              \
  virtual NewickNode *makeNode(string nm,                 \
		     double bl,				  \
		     NewickNode *p,			  \
		     vector<NewickNode *> *dl) = 0;	  \
  void setParentPointer(NewickNode *np, NewickNode *pp);  \
  NewickNode *theTree;	                                  \
  NewickFlexLexer lexer;                                  \

%token SIGNED_NUMBER, UNSIGNED_NUMBER, UNQUOTED_LABEL, QUOTED_LABEL
%type <node> tree subtree  
%type <nodePointerList> subtree_list descendant_list
%type <label> grammar_label root_label internal_node_label leaf_label
%type <value> branch_length

%start tree


%%
tree                       : descendant_list ';'
			   {theTree = makeNode("", 0, NULL, $1);
			   delete $1;}

                           | descendant_list root_label ';'
			   {theTree = makeNode(*$2, 0, NULL, $1);
			   delete $1;
			   delete $2;}

		           | descendant_list branch_length ';'
			   {theTree = makeNode("", $2, NULL, $1);
			   delete $1;}

			   | descendant_list root_label branch_length ';'
			   {theTree = makeNode(*$2, $3, NULL, $1);
			   delete $1;
			   delete $2;};

descendant_list            : '(' subtree_list ')'
			   {$$ = $2;};


subtree_list               : subtree
			   {$$ = new vector<NewickNode *>;
			    $$->push_back($1) };			   

                           | subtree_list ',' subtree
			   { $1->push_back($3);
			   $$ =	$1;};


subtree                    : descendant_list 
                           { $$ = makeNode("", 0, NULL, $1);
			   delete $1;}

                           | descendant_list internal_node_label
			   { $$ = makeNode(*$2, 0, NULL, $1);
			   delete $1;
			   delete $2;}

                           | descendant_list branch_length 
			   { $$ = makeNode("", $2, NULL, $1);
			   delete $1;}

                           | descendant_list internal_node_label branch_length
			   { $$ = makeNode(*$2, $3, NULL, $1);
			   delete $1;
			   delete $2;}

                           | leaf_label 
			   { $$ = makeNode(*$1, 0, NULL, NULL);
			   delete $1;}

                           | leaf_label branch_length
			   { $$ = makeNode(*$1, $2, NULL, NULL);
                           delete $1;};

root_label                 : grammar_label
			   {$$ = $1;};


internal_node_label        : grammar_label
			   {$$ = $1;};


leaf_label                 : grammar_label
			   {$$ = $1;};


grammar_label              : UNQUOTED_LABEL 
			   {$$ = new string(lexer.YYText() );}

			   | QUOTED_LABEL
			   {$$ = new string(lexer.YYText() );};


branch_length              : SIGNED_NUMBER 
			   {$$ = atof(lexer.YYText() ); }

			   | UNSIGNED_NUMBER
   			   {$$ = atof(lexer.YYText() ); };
