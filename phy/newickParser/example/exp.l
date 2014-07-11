%{ 
#include "parser.h"
%}
%option yylineno
%%

[ \t]                       ;
[0-9]+                      return(Parser::INT);

"."[0-9]*                   |                                        
[0-9]+("."[0-9]*)?          return(Parser::DOUBLE);

exit                        |
quit                        return (Parser::DONE);

.|\n                        return (*yytext);

