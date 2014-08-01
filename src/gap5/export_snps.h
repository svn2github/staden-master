#ifndef _EXPORT_SNPS_H_
#define _EXPORT_SNPS_H_

#include <tcl.h>

int tcl_export_snps(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[]);

#endif /* _EXPORT_SNPS_H_ */
