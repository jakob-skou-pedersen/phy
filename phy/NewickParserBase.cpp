/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define YY_NewickParserBase_h_included

/*  A Bison++ parser, made from newick.y  */

 /* with Bison++ version bison++ version 1.21-45, adapted from GNU Bison by coetmeur@icdc.fr
  */


#line 1 "/usr/local/lib/bison.cc"
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Bob Corbett and Richard Stallman

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* HEADER SECTION */
#if defined( _MSDOS ) || defined(MSDOS) || defined(__MSDOS__) 
#define __MSDOS_AND_ALIKE
#endif
#if defined(_WINDOWS) && defined(_MSC_VER)
#define __HAVE_NO_ALLOCA
#define __MSDOS_AND_ALIKE
#endif

#ifndef alloca
#if defined( __GNUC__)
#define alloca __builtin_alloca

#elif (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc)  || defined (__sgi)
#include <alloca.h>

#elif defined (__MSDOS_AND_ALIKE)
#include <malloc.h>
#ifndef __TURBOC__
/* MS C runtime lib */
#define alloca _alloca
#endif

#elif defined(_AIX)
#include <malloc.h>
#pragma alloca

#elif defined(__hpux)
#ifdef __cplusplus
extern "C" {
void *alloca (unsigned int);
};
#else /* not __cplusplus */
void *alloca ();
#endif /* not __cplusplus */

#endif /* not _AIX  not MSDOS, or __TURBOC__ or _AIX, not sparc.  */
#endif /* alloca not defined.  */
#ifdef c_plusplus
#ifndef __cplusplus
#define __cplusplus
#endif
#endif
#ifdef __cplusplus
#ifndef YY_USE_CLASS
#define YY_USE_CLASS
#endif
#else
#ifndef __STDC__
#define const
#endif
#endif
#include <stdio.h>
#define YYBISON 1  

/* #line 73 "/usr/local/lib/bison.cc" */
#line 85 "NewickParserBase.cpp"
#line 3 "newick.y"

#include <iostream>
#include <sstream>
#undef  yyFlexLexer
#define yyFlexLexer NewickFlexLexer
#include <FlexLexer.h>
#include "phy/NewickNode.h"

#line 12 "newick.y"
typedef union {
  NewickNode *node;
  vector<NewickNode *> *nodePointerList;
  string *label;
  double value;
  } yy_NewickParserBase_stype;
#define YY_NewickParserBase_STYPE yy_NewickParserBase_stype
#define YY_NewickParserBase_LEX_BODY   {return lexer.yylex(); }
#define YY_NewickParserBase_ERROR_BODY                                         \
{cerr << "syntax error in newick tree format before: \""  \
      << lexer.YYText()                                   \
      << "\"."                                            \
      << endl;}
#define YY_NewickParserBase_MEMBERS                                            \
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


#line 73 "/usr/local/lib/bison.cc"
/* %{ and %header{ and %union, during decl */
#define YY_NewickParserBase_BISON 1
#ifndef YY_NewickParserBase_COMPATIBILITY
#ifndef YY_USE_CLASS
#define  YY_NewickParserBase_COMPATIBILITY 1
#else
#define  YY_NewickParserBase_COMPATIBILITY 0
#endif
#endif

#if YY_NewickParserBase_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YYLTYPE
#ifndef YY_NewickParserBase_LTYPE
#define YY_NewickParserBase_LTYPE YYLTYPE
#endif
#endif
#ifdef YYSTYPE
#ifndef YY_NewickParserBase_STYPE 
#define YY_NewickParserBase_STYPE YYSTYPE
#endif
#endif
#ifdef YYDEBUG
#ifndef YY_NewickParserBase_DEBUG
#define  YY_NewickParserBase_DEBUG YYDEBUG
#endif
#endif
#ifdef YY_NewickParserBase_STYPE
#ifndef yystype
#define yystype YY_NewickParserBase_STYPE
#endif
#endif
/* use goto to be compatible */
#ifndef YY_NewickParserBase_USE_GOTO
#define YY_NewickParserBase_USE_GOTO 1
#endif
#endif

/* use no goto to be clean in C++ */
#ifndef YY_NewickParserBase_USE_GOTO
#define YY_NewickParserBase_USE_GOTO 0
#endif

#ifndef YY_NewickParserBase_PURE

/* #line 117 "/usr/local/lib/bison.cc" */
#line 169 "NewickParserBase.cpp"

#line 117 "/usr/local/lib/bison.cc"
/*  YY_NewickParserBase_PURE */
#endif

/* section apres lecture def, avant lecture grammaire S2 */

/* #line 121 "/usr/local/lib/bison.cc" */
#line 178 "NewickParserBase.cpp"

#line 121 "/usr/local/lib/bison.cc"
/* prefix */
#ifndef YY_NewickParserBase_DEBUG

/* #line 123 "/usr/local/lib/bison.cc" */
#line 185 "NewickParserBase.cpp"
#define YY_NewickParserBase_DEBUG 1

#line 123 "/usr/local/lib/bison.cc"
/* YY_NewickParserBase_DEBUG */
#endif


#ifndef YY_NewickParserBase_LSP_NEEDED

/* #line 128 "/usr/local/lib/bison.cc" */
#line 196 "NewickParserBase.cpp"

#line 128 "/usr/local/lib/bison.cc"
 /* YY_NewickParserBase_LSP_NEEDED*/
#endif



/* DEFAULT LTYPE*/
#ifdef YY_NewickParserBase_LSP_NEEDED
#ifndef YY_NewickParserBase_LTYPE
typedef
  struct yyltype
    {
      int timestamp;
      int first_line;
      int first_column;
      int last_line;
      int last_column;
      char *text;
   }
  yyltype;

#define YY_NewickParserBase_LTYPE yyltype
#endif
#endif
/* DEFAULT STYPE*/
      /* We used to use `unsigned long' as YY_NewickParserBase_STYPE on MSDOS,
	 but it seems better to be consistent.
	 Most programs should declare their own type anyway.  */

#ifndef YY_NewickParserBase_STYPE
#define YY_NewickParserBase_STYPE int
#endif
/* DEFAULT MISCELANEOUS */
#ifndef YY_NewickParserBase_PARSE
#define YY_NewickParserBase_PARSE yyparse
#endif
#ifndef YY_NewickParserBase_LEX
#define YY_NewickParserBase_LEX yylex
#endif
#ifndef YY_NewickParserBase_LVAL
#define YY_NewickParserBase_LVAL yylval
#endif
#ifndef YY_NewickParserBase_LLOC
#define YY_NewickParserBase_LLOC yylloc
#endif
#ifndef YY_NewickParserBase_CHAR
#define YY_NewickParserBase_CHAR yychar
#endif
#ifndef YY_NewickParserBase_NERRS
#define YY_NewickParserBase_NERRS yynerrs
#endif
#ifndef YY_NewickParserBase_DEBUG_FLAG
#define YY_NewickParserBase_DEBUG_FLAG yydebug
#endif
#ifndef YY_NewickParserBase_ERROR
#define YY_NewickParserBase_ERROR yyerror
#endif
#ifndef YY_NewickParserBase_PARSE_PARAM
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
#define YY_NewickParserBase_PARSE_PARAM
#ifndef YY_NewickParserBase_PARSE_PARAM_DEF
#define YY_NewickParserBase_PARSE_PARAM_DEF
#endif
#endif
#endif
#endif
#ifndef YY_NewickParserBase_PARSE_PARAM
#define YY_NewickParserBase_PARSE_PARAM void
#endif
#endif
#if YY_NewickParserBase_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YY_NewickParserBase_LTYPE
#ifndef YYLTYPE
#define YYLTYPE YY_NewickParserBase_LTYPE
#else
/* WARNING obsolete !!! user defined YYLTYPE not reported into generated header */
#endif
#endif
#ifndef YYSTYPE
#define YYSTYPE YY_NewickParserBase_STYPE
#else
/* WARNING obsolete !!! user defined YYSTYPE not reported into generated header */
#endif
#ifdef YY_NewickParserBase_PURE
#ifndef YYPURE
#define YYPURE YY_NewickParserBase_PURE
#endif
#endif
#ifdef YY_NewickParserBase_DEBUG
#ifndef YYDEBUG
#define YYDEBUG YY_NewickParserBase_DEBUG 
#endif
#endif
#ifndef YY_NewickParserBase_ERROR_VERBOSE
#ifdef YYERROR_VERBOSE
#define YY_NewickParserBase_ERROR_VERBOSE YYERROR_VERBOSE
#endif
#endif
#ifndef YY_NewickParserBase_LSP_NEEDED
#ifdef YYLSP_NEEDED
#define YY_NewickParserBase_LSP_NEEDED YYLSP_NEEDED
#endif
#endif
#endif
#ifndef YY_USE_CLASS
/* TOKEN C */

/* #line 236 "/usr/local/lib/bison.cc" */
#line 309 "NewickParserBase.cpp"
#define	SIGNED_NUMBER	258
#define	UNSIGNED_NUMBER	259
#define	UNQUOTED_LABEL	260
#define	QUOTED_LABEL	261


#line 236 "/usr/local/lib/bison.cc"
 /* #defines tokens */
#else
/* CLASS */
#ifndef YY_NewickParserBase_CLASS
#define YY_NewickParserBase_CLASS NewickParserBase
#endif
#ifndef YY_NewickParserBase_INHERIT
#define YY_NewickParserBase_INHERIT
#endif
#ifndef YY_NewickParserBase_MEMBERS
#define YY_NewickParserBase_MEMBERS 
#endif
#ifndef YY_NewickParserBase_LEX_BODY
#define YY_NewickParserBase_LEX_BODY  
#endif
#ifndef YY_NewickParserBase_ERROR_BODY
#define YY_NewickParserBase_ERROR_BODY  
#endif
#ifndef YY_NewickParserBase_CONSTRUCTOR_PARAM
#define YY_NewickParserBase_CONSTRUCTOR_PARAM
#endif
#ifndef YY_NewickParserBase_CONSTRUCTOR_CODE
#define YY_NewickParserBase_CONSTRUCTOR_CODE
#endif
#ifndef YY_NewickParserBase_CONSTRUCTOR_INIT
#define YY_NewickParserBase_CONSTRUCTOR_INIT
#endif
/* choose between enum and const */
#ifndef YY_NewickParserBase_USE_CONST_TOKEN
#define YY_NewickParserBase_USE_CONST_TOKEN 0
/* yes enum is more compatible with flex,  */
/* so by default we use it */ 
#endif
#if YY_NewickParserBase_USE_CONST_TOKEN != 0
#ifndef YY_NewickParserBase_ENUM_TOKEN
#define YY_NewickParserBase_ENUM_TOKEN yy_NewickParserBase_enum_token
#endif
#endif

class YY_NewickParserBase_CLASS YY_NewickParserBase_INHERIT
{
public: 
#if YY_NewickParserBase_USE_CONST_TOKEN != 0
/* static const int token ... */

/* #line 280 "/usr/local/lib/bison.cc" */
#line 363 "NewickParserBase.cpp"
static const int SIGNED_NUMBER;
static const int UNSIGNED_NUMBER;
static const int UNQUOTED_LABEL;
static const int QUOTED_LABEL;


#line 280 "/usr/local/lib/bison.cc"
 /* decl const */
#else
enum YY_NewickParserBase_ENUM_TOKEN { YY_NewickParserBase_NULL_TOKEN=0

/* #line 283 "/usr/local/lib/bison.cc" */
#line 376 "NewickParserBase.cpp"
	,SIGNED_NUMBER=258
	,UNSIGNED_NUMBER=259
	,UNQUOTED_LABEL=260
	,QUOTED_LABEL=261


#line 283 "/usr/local/lib/bison.cc"
 /* enum token */
     }; /* end of enum declaration */
#endif
public:
 int YY_NewickParserBase_PARSE (YY_NewickParserBase_PARSE_PARAM);
 virtual void YY_NewickParserBase_ERROR(const char *msg) YY_NewickParserBase_ERROR_BODY;
#ifdef YY_NewickParserBase_PURE
#ifdef YY_NewickParserBase_LSP_NEEDED
 virtual int  YY_NewickParserBase_LEX (YY_NewickParserBase_STYPE *YY_NewickParserBase_LVAL,YY_NewickParserBase_LTYPE *YY_NewickParserBase_LLOC) YY_NewickParserBase_LEX_BODY;
#else
 virtual int  YY_NewickParserBase_LEX (YY_NewickParserBase_STYPE *YY_NewickParserBase_LVAL) YY_NewickParserBase_LEX_BODY;
#endif
#else
 virtual int YY_NewickParserBase_LEX() YY_NewickParserBase_LEX_BODY;
 YY_NewickParserBase_STYPE YY_NewickParserBase_LVAL;
#ifdef YY_NewickParserBase_LSP_NEEDED
 YY_NewickParserBase_LTYPE YY_NewickParserBase_LLOC;
#endif
 int   YY_NewickParserBase_NERRS;
 int    YY_NewickParserBase_CHAR;
#endif
#if YY_NewickParserBase_DEBUG != 0
 int YY_NewickParserBase_DEBUG_FLAG;   /*  nonzero means print parse trace     */
#endif
public:
 YY_NewickParserBase_CLASS(YY_NewickParserBase_CONSTRUCTOR_PARAM);
public:
 YY_NewickParserBase_MEMBERS 
};
/* other declare folow */
#if YY_NewickParserBase_USE_CONST_TOKEN != 0

/* #line 314 "/usr/local/lib/bison.cc" */
#line 417 "NewickParserBase.cpp"
const int YY_NewickParserBase_CLASS::SIGNED_NUMBER=258;
const int YY_NewickParserBase_CLASS::UNSIGNED_NUMBER=259;
const int YY_NewickParserBase_CLASS::UNQUOTED_LABEL=260;
const int YY_NewickParserBase_CLASS::QUOTED_LABEL=261;


#line 314 "/usr/local/lib/bison.cc"
 /* const YY_NewickParserBase_CLASS::token */
#endif
/*apres const  */
YY_NewickParserBase_CLASS::YY_NewickParserBase_CLASS(YY_NewickParserBase_CONSTRUCTOR_PARAM) YY_NewickParserBase_CONSTRUCTOR_INIT
{
#if YY_NewickParserBase_DEBUG != 0
YY_NewickParserBase_DEBUG_FLAG=0;
#endif
YY_NewickParserBase_CONSTRUCTOR_CODE;
};
#endif

/* #line 325 "/usr/local/lib/bison.cc" */
#line 438 "NewickParserBase.cpp"


#define	YYFINAL		30
#define	YYFLAG		-32768
#define	YYNTBASE	11

#define YYTRANSLATE(x) ((unsigned)(x) <= 261 ? yytranslate[x] : 20)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     8,
     9,     2,     2,    10,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     7,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
     6
};

#if YY_NewickParserBase_DEBUG != 0
static const short yyprhs[] = {     0,
     0,     3,     7,    11,    16,    20,    22,    26,    28,    31,
    34,    38,    40,    43,    45,    47,    49,    51,    53,    55
};

static const short yyrhs[] = {    12,
     7,     0,    12,    15,     7,     0,    12,    19,     7,     0,
    12,    15,    19,     7,     0,     8,    13,     9,     0,    14,
     0,    13,    10,    14,     0,    12,     0,    12,    16,     0,
    12,    19,     0,    12,    16,    19,     0,    17,     0,    17,
    19,     0,    18,     0,    18,     0,    18,     0,     5,     0,
     6,     0,     3,     0,     4,     0
};

#endif

#if YY_NewickParserBase_DEBUG != 0
static const short yyrline[] = { 0,
    48,    52,    57,    61,    66,    70,    74,    79,    83,    88,
    92,    97,   101,   105,   109,   113,   117,   120,   124,   127
};

static const char * const yytname[] = {   "$","error","$illegal.","SIGNED_NUMBER",
"UNSIGNED_NUMBER","UNQUOTED_LABEL","QUOTED_LABEL","';'","'('","')'","','","tree",
"descendant_list","subtree_list","subtree","root_label","internal_node_label",
"leaf_label","grammar_label","branch_length",""
};
#endif

static const short yyr1[] = {     0,
    11,    11,    11,    11,    12,    13,    13,    14,    14,    14,
    14,    14,    14,    15,    16,    17,    18,    18,    19,    19
};

static const short yyr2[] = {     0,
     2,     3,     3,     4,     3,     1,     3,     1,     2,     2,
     3,     1,     2,     1,     1,     1,     1,     1,     1,     1
};

static const short yydefact[] = {     0,
     0,     0,    17,    18,     8,     0,     6,    12,    16,    19,
    20,     1,     0,    14,     0,     9,    15,    10,     5,     0,
    13,     2,     0,     3,    11,     7,     4,     0,     0,     0
};

static const short yydefgoto[] = {    28,
     5,     6,     7,    13,    16,     8,     9,    15
};

static const short yypact[] = {    -4,
     1,     9,-32768,-32768,    14,    13,-32768,    22,-32768,-32768,
-32768,-32768,    -2,-32768,     3,    22,-32768,-32768,-32768,     1,
-32768,-32768,    20,-32768,-32768,-32768,-32768,    28,    29,-32768
};

static const short yypgoto[] = {-32768,
    30,-32768,    11,-32768,-32768,-32768,    19,    -5
};


#define	YYLAST		31


static const short yytable[] = {    18,
    10,    11,    21,     1,    22,     3,     4,    23,     1,    24,
    25,    10,    11,     3,     4,    12,    10,    11,     3,     4,
    14,    19,    20,    17,    10,    11,    27,    29,    30,     2,
    26
};

static const short yycheck[] = {     5,
     3,     4,     8,     8,     7,     5,     6,    13,     8,     7,
    16,     3,     4,     5,     6,     7,     3,     4,     5,     6,
     2,     9,    10,     5,     3,     4,     7,     0,     0,     0,
    20
};

#line 325 "/usr/local/lib/bison.cc"
 /* fattrs + tables */

/* parser code folow  */


/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: dollar marks section change
   the next  is replaced by the list of actions, each action
   as one case of the switch.  */ 

#if YY_NewickParserBase_USE_GOTO != 0
/* 
 SUPRESSION OF GOTO : on some C++ compiler (sun c++)
  the goto is strictly forbidden if any constructor/destructor
  is used in the whole function (very stupid isn't it ?)
 so goto are to be replaced with a 'while/switch/case construct'
 here are the macro to keep some apparent compatibility
*/
#define YYGOTO(lb) {yy_gotostate=lb;continue;}
#define YYBEGINGOTO  enum yy_labels yy_gotostate=yygotostart; \
                     for(;;) switch(yy_gotostate) { case yygotostart: {
#define YYLABEL(lb) } case lb: {
#define YYENDGOTO } } 
#define YYBEGINDECLARELABEL enum yy_labels {yygotostart
#define YYDECLARELABEL(lb) ,lb
#define YYENDDECLARELABEL  };
#else
/* macro to keep goto */
#define YYGOTO(lb) goto lb
#define YYBEGINGOTO 
#define YYLABEL(lb) lb:
#define YYENDGOTO
#define YYBEGINDECLARELABEL 
#define YYDECLARELABEL(lb)
#define YYENDDECLARELABEL 
#endif
/* LABEL DECLARATION */
YYBEGINDECLARELABEL
  YYDECLARELABEL(yynewstate)
  YYDECLARELABEL(yybackup)
/* YYDECLARELABEL(yyresume) */
  YYDECLARELABEL(yydefault)
  YYDECLARELABEL(yyreduce)
  YYDECLARELABEL(yyerrlab)   /* here on detecting error */
  YYDECLARELABEL(yyerrlab1)   /* here on error raised explicitly by an action */
  YYDECLARELABEL(yyerrdefault)  /* current state does not do anything special for the error token. */
  YYDECLARELABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */
  YYDECLARELABEL(yyerrhandle)  
YYENDDECLARELABEL
/* ALLOCA SIMULATION */
/* __HAVE_NO_ALLOCA */
#ifdef __HAVE_NO_ALLOCA
int __alloca_free_ptr(char *ptr,char *ref)
{if(ptr!=ref) free(ptr);
 return 0;}

#define __ALLOCA_alloca(size) malloc(size)
#define __ALLOCA_free(ptr,ref) __alloca_free_ptr((char *)ptr,(char *)ref)

#ifdef YY_NewickParserBase_LSP_NEEDED
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		    __ALLOCA_free(yyls,yylsa)+\
		   (num))
#else
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		   (num))
#endif
#else
#define __ALLOCA_return(num) return(num)
#define __ALLOCA_alloca(size) alloca(size)
#define __ALLOCA_free(ptr,ref) 
#endif

/* ENDALLOCA SIMULATION */

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (YY_NewickParserBase_CHAR = YYEMPTY)
#define YYEMPTY         -2
#define YYEOF           0
#define YYACCEPT        __ALLOCA_return(0)
#define YYABORT         __ALLOCA_return(1)
#define YYERROR         YYGOTO(yyerrlab1)
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL          YYGOTO(yyerrlab)
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do                                                              \
  if (YY_NewickParserBase_CHAR == YYEMPTY && yylen == 1)                               \
    { YY_NewickParserBase_CHAR = (token), YY_NewickParserBase_LVAL = (value);                 \
      yychar1 = YYTRANSLATE (YY_NewickParserBase_CHAR);                                \
      YYPOPSTACK;                                               \
      YYGOTO(yybackup);                                            \
    }                                                           \
  else                                                          \
    { YY_NewickParserBase_ERROR ("syntax error: cannot back up"); YYERROR; }   \
while (0)

#define YYTERROR        1
#define YYERRCODE       256

#ifndef YY_NewickParserBase_PURE
/* UNPURE */
#define YYLEX           YY_NewickParserBase_LEX()
#ifndef YY_USE_CLASS
/* If nonreentrant, and not class , generate the variables here */
int     YY_NewickParserBase_CHAR;                      /*  the lookahead symbol        */
YY_NewickParserBase_STYPE      YY_NewickParserBase_LVAL;              /*  the semantic value of the */
				/*  lookahead symbol    */
int YY_NewickParserBase_NERRS;                 /*  number of parse errors so far */
#ifdef YY_NewickParserBase_LSP_NEEDED
YY_NewickParserBase_LTYPE YY_NewickParserBase_LLOC;   /*  location data for the lookahead     */
			/*  symbol                              */
#endif
#endif


#else
/* PURE */
#ifdef YY_NewickParserBase_LSP_NEEDED
#define YYLEX           YY_NewickParserBase_LEX(&YY_NewickParserBase_LVAL, &YY_NewickParserBase_LLOC)
#else
#define YYLEX           YY_NewickParserBase_LEX(&YY_NewickParserBase_LVAL)
#endif
#endif
#ifndef YY_USE_CLASS
#if YY_NewickParserBase_DEBUG != 0
int YY_NewickParserBase_DEBUG_FLAG;                    /*  nonzero means print parse trace     */
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif
#endif



/*  YYINITDEPTH indicates the initial size of the parser's stacks       */

#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif


#if __GNUC__ > 1                /* GNU C and GNU C++ define this.  */
#define __yy_bcopy(FROM,TO,COUNT)       __builtin_memcpy(TO,FROM,COUNT)
#else                           /* not GNU C or C++ */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */

#ifdef __cplusplus
static void __yy_bcopy (char *from, char *to, int count)
#else
#ifdef __STDC__
static void __yy_bcopy (char *from, char *to, int count)
#else
static void __yy_bcopy (from, to, count)
     char *from;
     char *to;
     int count;
#endif
#endif
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}
#endif

int
#ifdef YY_USE_CLASS
 YY_NewickParserBase_CLASS::
#endif
     YY_NewickParserBase_PARSE(YY_NewickParserBase_PARSE_PARAM)
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
/* parameter definition without protypes */
YY_NewickParserBase_PARSE_PARAM_DEF
#endif
#endif
#endif
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YY_NewickParserBase_STYPE *yyvsp;
  int yyerrstatus;      /*  number of tokens to shift before error messages enabled */
  int yychar1=0;          /*  lookahead token as an internal (translated) token number */

  short yyssa[YYINITDEPTH];     /*  the state stack                     */
  YY_NewickParserBase_STYPE yyvsa[YYINITDEPTH];        /*  the semantic value stack            */

  short *yyss = yyssa;          /*  refer to the stacks thru separate pointers */
  YY_NewickParserBase_STYPE *yyvs = yyvsa;     /*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YY_NewickParserBase_LSP_NEEDED
  YY_NewickParserBase_LTYPE yylsa[YYINITDEPTH];        /*  the location stack                  */
  YY_NewickParserBase_LTYPE *yyls = yylsa;
  YY_NewickParserBase_LTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;

#ifdef YY_NewickParserBase_PURE
  int YY_NewickParserBase_CHAR;
  YY_NewickParserBase_STYPE YY_NewickParserBase_LVAL;
  int YY_NewickParserBase_NERRS;
#ifdef YY_NewickParserBase_LSP_NEEDED
  YY_NewickParserBase_LTYPE YY_NewickParserBase_LLOC;
#endif
#endif

  YY_NewickParserBase_STYPE yyval;             /*  the variable used to return         */
				/*  semantic values from the action     */
				/*  routines                            */

  int yylen;
/* start loop, in which YYGOTO may be used. */
YYBEGINGOTO

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    fprintf(stderr, "Starting parse\n");
#endif
  yystate = 0;
  yyerrstatus = 0;
  YY_NewickParserBase_NERRS = 0;
  YY_NewickParserBase_CHAR = YYEMPTY;          /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YY_NewickParserBase_LSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
YYLABEL(yynewstate)

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YY_NewickParserBase_STYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YY_NewickParserBase_LSP_NEEDED
      YY_NewickParserBase_LTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YY_NewickParserBase_LSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YY_NewickParserBase_LSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  YY_NewickParserBase_ERROR("parser stack overflow");
	  __ALLOCA_return(2);
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
      yyss = (short *) __ALLOCA_alloca (yystacksize * sizeof (*yyssp));
      __yy_bcopy ((char *)yyss1, (char *)yyss, size * sizeof (*yyssp));
      __ALLOCA_free(yyss1,yyssa);
      yyvs = (YY_NewickParserBase_STYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yyvsp));
      __yy_bcopy ((char *)yyvs1, (char *)yyvs, size * sizeof (*yyvsp));
      __ALLOCA_free(yyvs1,yyvsa);
#ifdef YY_NewickParserBase_LSP_NEEDED
      yyls = (YY_NewickParserBase_LTYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yylsp));
      __yy_bcopy ((char *)yyls1, (char *)yyls, size * sizeof (*yylsp));
      __ALLOCA_free(yyls1,yylsa);
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YY_NewickParserBase_LSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YY_NewickParserBase_DEBUG != 0
      if (YY_NewickParserBase_DEBUG_FLAG)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  YYGOTO(yybackup);
YYLABEL(yybackup)

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* YYLABEL(yyresume) */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yydefault);

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (YY_NewickParserBase_CHAR == YYEMPTY)
    {
#if YY_NewickParserBase_DEBUG != 0
      if (YY_NewickParserBase_DEBUG_FLAG)
	fprintf(stderr, "Reading a token: ");
#endif
      YY_NewickParserBase_CHAR = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (YY_NewickParserBase_CHAR <= 0)           /* This means end of input. */
    {
      yychar1 = 0;
      YY_NewickParserBase_CHAR = YYEOF;                /* Don't call YYLEX any more */

#if YY_NewickParserBase_DEBUG != 0
      if (YY_NewickParserBase_DEBUG_FLAG)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(YY_NewickParserBase_CHAR);

#if YY_NewickParserBase_DEBUG != 0
      if (YY_NewickParserBase_DEBUG_FLAG)
	{
	  fprintf (stderr, "Next token is %d (%s", YY_NewickParserBase_CHAR, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, YY_NewickParserBase_CHAR, YY_NewickParserBase_LVAL);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    YYGOTO(yydefault);

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrlab);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrlab);

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    fprintf(stderr, "Shifting token %d (%s), ", YY_NewickParserBase_CHAR, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (YY_NewickParserBase_CHAR != YYEOF)
    YY_NewickParserBase_CHAR = YYEMPTY;

  *++yyvsp = YY_NewickParserBase_LVAL;
#ifdef YY_NewickParserBase_LSP_NEEDED
  *++yylsp = YY_NewickParserBase_LLOC;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  YYGOTO(yynewstate);

/* Do the default action for the current state.  */
YYLABEL(yydefault)

  yyn = yydefact[yystate];
  if (yyn == 0)
    YYGOTO(yyerrlab);

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
YYLABEL(yyreduce)
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


/* #line 811 "/usr/local/lib/bison.cc" */
#line 1044 "NewickParserBase.cpp"

  switch (yyn) {

case 1:
#line 49 "newick.y"
{theTree = makeNode("", 0, NULL, yyvsp[-1].nodePointerList);
			   delete yyvsp[-1].nodePointerList;;
    break;}
case 2:
#line 53 "newick.y"
{theTree = makeNode(*yyvsp[-1].label, 0, NULL, yyvsp[-2].nodePointerList);
			   delete yyvsp[-2].nodePointerList;
			   delete yyvsp[-1].label;;
    break;}
case 3:
#line 58 "newick.y"
{theTree = makeNode("", yyvsp[-1].value, NULL, yyvsp[-2].nodePointerList);
			   delete yyvsp[-2].nodePointerList;;
    break;}
case 4:
#line 62 "newick.y"
{theTree = makeNode(*yyvsp[-2].label, yyvsp[-1].value, NULL, yyvsp[-3].nodePointerList);
			   delete yyvsp[-3].nodePointerList;
			   delete yyvsp[-2].label;;
    break;}
case 5:
#line 67 "newick.y"
{yyval.nodePointerList = yyvsp[-1].nodePointerList;;
    break;}
case 6:
#line 71 "newick.y"
{yyval.nodePointerList = new vector<NewickNode *>;
			    yyval.nodePointerList->push_back(yyvsp[0].node) ;
    break;}
case 7:
#line 75 "newick.y"
{ yyvsp[-2].nodePointerList->push_back(yyvsp[0].node);
			   yyval.nodePointerList =	yyvsp[-2].nodePointerList;;
    break;}
case 8:
#line 80 "newick.y"
{ yyval.node = makeNode("", 0, NULL, yyvsp[0].nodePointerList);
			   delete yyvsp[0].nodePointerList;;
    break;}
case 9:
#line 84 "newick.y"
{ yyval.node = makeNode(*yyvsp[0].label, 0, NULL, yyvsp[-1].nodePointerList);
			   delete yyvsp[-1].nodePointerList;
			   delete yyvsp[0].label;;
    break;}
case 10:
#line 89 "newick.y"
{ yyval.node = makeNode("", yyvsp[0].value, NULL, yyvsp[-1].nodePointerList);
			   delete yyvsp[-1].nodePointerList;;
    break;}
case 11:
#line 93 "newick.y"
{ yyval.node = makeNode(*yyvsp[-1].label, yyvsp[0].value, NULL, yyvsp[-2].nodePointerList);
			   delete yyvsp[-2].nodePointerList;
			   delete yyvsp[-1].label;;
    break;}
case 12:
#line 98 "newick.y"
{ yyval.node = makeNode(*yyvsp[0].label, 0, NULL, NULL);
			   delete yyvsp[0].label;;
    break;}
case 13:
#line 102 "newick.y"
{ yyval.node = makeNode(*yyvsp[-1].label, yyvsp[0].value, NULL, NULL);
                           delete yyvsp[-1].label;;
    break;}
case 14:
#line 106 "newick.y"
{yyval.label = yyvsp[0].label;;
    break;}
case 15:
#line 110 "newick.y"
{yyval.label = yyvsp[0].label;;
    break;}
case 16:
#line 114 "newick.y"
{yyval.label = yyvsp[0].label;;
    break;}
case 17:
#line 118 "newick.y"
{yyval.label = new string(lexer.YYText() );;
    break;}
case 18:
#line 121 "newick.y"
{yyval.label = new string(lexer.YYText() );;
    break;}
case 19:
#line 125 "newick.y"
{yyval.value = atof(lexer.YYText() ); ;
    break;}
case 20:
#line 128 "newick.y"
{yyval.value = atof(lexer.YYText() ); ;
    break;}
}

#line 811 "/usr/local/lib/bison.cc"
   /* the action file gets copied in in place of this dollarsign  */
  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YY_NewickParserBase_LSP_NEEDED
  yylsp -= yylen;
#endif

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YY_NewickParserBase_LSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = YY_NewickParserBase_LLOC.first_line;
      yylsp->first_column = YY_NewickParserBase_LLOC.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  YYGOTO(yynewstate);

YYLABEL(yyerrlab)   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++YY_NewickParserBase_NERRS;

#ifdef YY_NewickParserBase_ERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      YY_NewickParserBase_ERROR(msg);
	      free(msg);
	    }
	  else
	    YY_NewickParserBase_ERROR ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YY_NewickParserBase_ERROR_VERBOSE */
	YY_NewickParserBase_ERROR("parse error");
    }

  YYGOTO(yyerrlab1);
YYLABEL(yyerrlab1)   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (YY_NewickParserBase_CHAR == YYEOF)
	YYABORT;

#if YY_NewickParserBase_DEBUG != 0
      if (YY_NewickParserBase_DEBUG_FLAG)
	fprintf(stderr, "Discarding token %d (%s).\n", YY_NewickParserBase_CHAR, yytname[yychar1]);
#endif

      YY_NewickParserBase_CHAR = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;              /* Each real token shifted decrements this */

  YYGOTO(yyerrhandle);

YYLABEL(yyerrdefault)  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) YYGOTO(yydefault);
#endif

YYLABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YY_NewickParserBase_LSP_NEEDED
  yylsp--;
#endif

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

YYLABEL(yyerrhandle)

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yyerrdefault);

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    YYGOTO(yyerrdefault);

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrpop);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrpop);

  if (yyn == YYFINAL)
    YYACCEPT;

#if YY_NewickParserBase_DEBUG != 0
  if (YY_NewickParserBase_DEBUG_FLAG)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = YY_NewickParserBase_LVAL;
#ifdef YY_NewickParserBase_LSP_NEEDED
  *++yylsp = YY_NewickParserBase_LLOC;
#endif

  yystate = yyn;
  YYGOTO(yynewstate);
/* end loop, in which YYGOTO may be used. */
  YYENDGOTO
}

/* END */

/* #line 1010 "/usr/local/lib/bison.cc" */
#line 1348 "NewickParserBase.cpp"
#line 129 "newick.y"
