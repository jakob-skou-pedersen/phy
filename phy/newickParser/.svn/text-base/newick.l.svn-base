%{
#include "NewickParserBase.h"
%}


%%
[ \t]+          /* white space */;

\n                 /* newlines */;

\[[^\]]*\]         /* comment  */;

"("		     return '(';

")"                  return ')';

","     	     return ',';

";"     	     return ';';

'[^'\n]+'	     {yytext[yyleng - 1] = 0; 
		     yytext = yytext + 1;
		     return NewickParserBase::QUOTED_LABEL;}

[^[\]()':;,\n ]+     return NewickParserBase::UNQUOTED_LABEL;

:[ \t\n]*[-+][0-9]+\.?[0-9]* {yytext = yytext+1; return NewickParserBase::SIGNED_NUMBER;}

:[ \t\n]*[0-9]+\.?[0-9]*     {yytext = yytext+1; 
                     return NewickParserBase::UNSIGNED_NUMBER;}
%%

int NewickFlexLexer::yywrap() {return 1;}
