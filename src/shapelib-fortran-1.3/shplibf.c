/*
 Fortran 90 bindings for shapelib (http://shapelib.maptools.org/)
 Copyright (C) Davide Cesari, 2007, dcesari <at> arpa.emr.it

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <alloca.h>
#include <string.h>

#include <libshp/shapefil.h>
#include "config.h"

#define strf2c(strf, strc, strl) strc = alloca(strl+1); strncpy(strc, strf, strl); strc[strl] = '\0'


SHPHandle FC_FUNC(shpopen, SHPOPEN)
  (char *pszShapeFile, char *pszAccess, int sfl, int al) {
  char *lpszShapeFile, *lpszAccess;

  strf2c(pszShapeFile, lpszShapeFile, sfl);
  strf2c(pszAccess, lpszAccess, al);
  return SHPOpen(lpszShapeFile, lpszAccess);
}


void FC_FUNC(shpgetinfo, SHPGETINFO)
  (SHPHandle *hshp, int *entities, int *shapetype,
   double *minbound, double *maxbound) {

  SHPGetInfo(*hshp, entities, shapetype, minbound, maxbound);
}


int FC_FUNC_(shpreadobject_int, SHPREADOBJECT_INT)
  (SHPHandle *hshp, int *iShape, void *ftnobject) {
  SHPObject *psObject;

  psObject = SHPReadObject(*hshp, *iShape);

    if (psObject) {
    FC_FUNC_(shpset_object, SHPSET_OBJECT)
      (ftnobject, &psObject,
       &psObject->nSHPType, &psObject->nShapeId, &psObject->nParts,
       psObject->panPartStart, psObject->panPartType, &psObject->nVertices,
       psObject->padfX, psObject->padfY, psObject->padfZ, psObject->padfM,
       &psObject->dfXMin, &psObject->dfYMin, &psObject->dfZMin, &psObject->dfMMin,
       &psObject->dfXMax, &psObject->dfYMax, &psObject->dfZMax, &psObject->dfMMax);
    return 0;
  } else {
    return 1;
  }
}
						 

void FC_FUNC(shpclose, SHPCLOSE)(SHPHandle *hshp) {

  return SHPClose(*hshp);
}


SHPHandle FC_FUNC(shpcreate, SHPCREATE)
  (char *pszShapeFile, int *nShapeType, int sfl) {
  char *lpszShapeFile;

  strf2c(pszShapeFile, lpszShapeFile, sfl);
  return SHPCreate(lpszShapeFile, *nShapeType);
}


int FC_FUNC_(shpcreatesimpleobject_int, SHPCREATESIMPLEOBJECT_INT)
  (int *nSHPType, int *nVertices, double *padfX, double *padfY, double *padfZ,
   void *ftnobject) {
  SHPObject *psObject;

  psObject = SHPCreateSimpleObject(*nSHPType, *nVertices, padfX, padfY, padfZ);
  if (psObject) {
    FC_FUNC_(shpset_object, SHPSET_OBJECT)
      (ftnobject, &psObject,
       &psObject->nSHPType, &psObject->nShapeId, &psObject->nParts,
       psObject->panPartStart, psObject->panPartType, &psObject->nVertices,
       psObject->padfX, psObject->padfY, psObject->padfZ, psObject->padfM,
       &psObject->dfXMin, &psObject->dfYMin, &psObject->dfZMin, &psObject->dfMMin,
       &psObject->dfXMax, &psObject->dfYMax, &psObject->dfZMax, &psObject->dfMMax);
    return 0;
  } else {
    return 1;
  }
}
						 

int FC_FUNC_(shpcreateobject_int, SHPCREATEOBJECT_INT)
  (int *nSHPType, int *iShape, int *nParts, int *panPartStart, int *panPartType,
   int *nVertices, double *padfX, double *padfY, double *padfZ, double *padfM,
   void *ftnobject) {
  SHPObject *psObject;

  psObject = SHPCreateObject(*nSHPType, *iShape, *nParts, panPartStart,
			     panPartType, *nVertices, padfX, padfY, padfZ, padfM);
  if (psObject) {
    FC_FUNC_(shpset_object, SHPSET_OBJECT)
      (ftnobject, &psObject,
       &psObject->nSHPType, &psObject->nShapeId, &psObject->nParts,
       psObject->panPartStart, psObject->panPartType, &psObject->nVertices,
       psObject->padfX, psObject->padfY, psObject->padfZ, psObject->padfM,
       &psObject->dfXMin, &psObject->dfYMin, &psObject->dfZMin, &psObject->dfMMin,
       &psObject->dfXMax, &psObject->dfYMax, &psObject->dfZMax, &psObject->dfMMax);
    return 0;
  } else {
    return 1;
  }
}
						 

void FC_FUNC_(shpcomputeextents_int, SHPCOMPUTEEXTENTS_INT)
  (SHPObject **psObject, void *ftnobject) {
  SHPObject *lpsObject;

  lpsObject = *psObject;
  SHPComputeExtents(lpsObject);

  FC_FUNC_(shpset_object, SHPSET_OBJECT)
    (ftnobject, &lpsObject,
     &lpsObject->nSHPType, &lpsObject->nShapeId, &lpsObject->nParts,
     lpsObject->panPartStart, lpsObject->panPartType, &lpsObject->nVertices,
     lpsObject->padfX, lpsObject->padfY, lpsObject->padfZ, lpsObject->padfM,
     &lpsObject->dfXMin, &lpsObject->dfYMin, &lpsObject->dfZMin, &lpsObject->dfMMin,
     &lpsObject->dfXMax, &lpsObject->dfYMax, &lpsObject->dfZMax, &lpsObject->dfMMax);

}


int FC_FUNC_(shpwriteobject_int, SHPWRITEOBJECT_INT)
  (SHPHandle *hSHP, int *iShape, SHPObject **psObject ) {

  return SHPWriteObject(*hSHP, *iShape, *psObject);

}


void FC_FUNC_(shpdestroyobject_int, SHPDESTROYOBJECT_INT)
  (SHPObject **psObject) {

  SHPDestroyObject(*psObject);

}


int FC_FUNC_(shprewindobject_int, SHPREWINDOBJECT_INT)
  (SHPHandle *hSHP, SHPObject *psObject) {

  return SHPRewindObject(*hSHP, psObject);

}
