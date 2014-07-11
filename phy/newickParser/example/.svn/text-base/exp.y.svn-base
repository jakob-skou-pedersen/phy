%name Parser 
%header{
#include <iostream>
#include <FlexLexer.h>
#include "stdlib.h"
%}

%union 
{
  int i;
  double d;
};

%token      INT DOUBLE DONE
%type   <i> intExpr
%type   <d> doubleExpr
    
%left   '+'
%left   '*'
%right  UnaryMinus
    
%define MEMBERS virtual ~Parser(){} \
 private:                           \
 yyFlexLexer lexer;
                  
%define LEX_BODY {return lexer.yylex();}
    
%define ERROR_BODY { cerr << "error encountered\n"; }

%%
lines:
lines
line
|
line
;
    
line:
intExpr
'\n'
{
  cerr << "int: " << $1 << endl;
}
|
doubleExpr
'\n'
{
  cerr << "double: " << $1 << endl;
}
|
'\n'
{
  cout << "Good bye\n";
  YYACCEPT;
}
|
error
'\n'
;
    
intExpr:
intExpr '*' intExpr
{
  $$ = $1 * $3;
}
|
intExpr '+' intExpr
{
  $$ = $1 + $3;
}
|
'(' intExpr ')'
{
  $$ = $2;
}
|
'-' intExpr         %prec UnaryMinus
{
  $$ = -$2;
}
|
INT
{
  $$ = atoi(lexer.YYText());
}
;
    
doubleExpr:
doubleExpr '*' doubleExpr
{
  $$ = $1 * $3;
}
|
doubleExpr '+' doubleExpr
{
  $$ = $1 + $3;
}
|
doubleExpr '*' intExpr
{
  $$ = $1 * $3;
}
|
doubleExpr '+' intExpr
{
  $$ = $1 + $3;
}
|
intExpr '*' doubleExpr
{
  $$ = $1 * $3;
}
|
intExpr '+' doubleExpr
{
  $$ = $1 + $3;
}
|
'(' doubleExpr ')'
{
  $$ = $2;
}
|
'-' doubleExpr         %prec UnaryMinus
{
  $$ = -$2;
}
|
DOUBLE
{
  $$ = atof(lexer.YYText());
}
;
%%
