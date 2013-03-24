#ifndef MATHEMATICA_H
#define	MATHEMATICA_H

#include <mathlink.h>

struct mathlink_env_t
{	
	MLEnvironment gEnvironment;
	MLINK theLink;
};

struct mathlink_env_t *init_mathlink(void);
void fini_mathlink(struct mathlink_env_t *met);
const char *mathlink_eval(struct mathlink_env_t *met,const char *fmt,...);
void mathlink_disown_string(struct mathlink_env_t *met,const char *theString);

#endif	/* MATHEMATICA_H */
