#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mathematica.h"

struct mathlink_env_t *init_mathlink(void)
{
	struct mathlink_env_t *met;
	int err;

	met=(struct mathlink_env_t *)(malloc(sizeof(struct mathlink_env_t)));
	if(!met)
		return NULL;

	met->gEnvironment=NULL;
	met->theLink=NULL;

	if(NULL==(met->gEnvironment=MLInitialize(0)))
		return NULL;

	met->theLink=MLOpenString(met->gEnvironment,"-linkmode launch -linkname '\"/Applications/Mathematica.app/Contents/MacOS/MathKernel\" -mathlink'",&err);

	return met;
}

void fini_mathlink(struct mathlink_env_t *met)
{
	MLPutFunction(met->theLink,"Exit",0);

	
	if(met)
	{
		if(met->theLink!=NULL)
			MLClose(met->theLink);

		if(met->gEnvironment!=NULL)
			MLDeinitialize(met->gEnvironment);
		
		free(met);
	}
}

const char *mathlink_eval(struct mathlink_env_t *met,const char *fmt,...)
{
	int thePacket;
	va_list ap;
	char *theExpression;

	if(!met)
		return NULL;

	va_start(ap,fmt);
	vasprintf(&theExpression,fmt,ap);
	va_end(ap);

	if(!theExpression)
		return NULL;

	MLNewPacket(met->theLink);
	MLPutFunction(met->theLink,"EvaluatePacket",1);
	MLPutFunction(met->theLink,"ToString",2);
	MLPutFunction(met->theLink,"ToExpression",1);
	MLPutString(met->theLink,theExpression);
	MLPutSymbol(met->theLink,"InputForm");
	MLEndPacket(met->theLink);

	while((thePacket=MLNextPacket(met->theLink))!=0)
	{
		if(thePacket==RETURNPKT)
		{
			const char *theResult;
			if(!MLGetString(met->theLink,&theResult))
				theResult=NULL;
			
			if(theExpression)
				free(theExpression);

			return theResult;
		}
	
		MLNewPacket(met->theLink);
	}

	if(theExpression)
		free(theExpression);

	return NULL;
}

void mathlink_disown_string(struct mathlink_env_t *met,const char *theString)
{
	MLDisownString(met->theLink,theString);
}
