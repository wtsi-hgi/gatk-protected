#include "Python.h"
#include "arrayobject.h"
#include <math.h>
#include <string.h>
#include "grmutils.h"
#include <unistd.h>

/* Module initialization. Stolen from http://docs.python.org/py3k/howto/cporting.html */

struct module_state {
 PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject* error_out(PyObject *m){
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyObject* grm_cleanup(PyObject *self){
 // release the DOSAGES matrix (the variable is global, defined in the header)
 free_Carrayptrs(DOSAGES);
 // remember to decref
 Py_DECREF(DOSAGES_PYMAT_OBJ);
 // return true 
 Py_RETURN_TRUE;
}

static PyMethodDef grmutilsmethods[] = {
 {"error_out",(PyCFunction)error_out,METH_NOARGS,NULL},
 {"loadDosage",grm_loadDosage,METH_VARARGS,NULL},
 {"calculateKinship",grm_calcdistance,METH_VARARGS,NULL},
 {"calculateKinshipNoMissing",grm_calcdistance_nomissing,METH_VARARGS,NULL},
 {"cleanup",grm_cleanup,METH_NOARGS,NULL},
 {NULL,NULL}
};

#if PY_MAJOR_VERSION >= 3

static int grmutils_traverse(PyObject *m,visitproc visit, void *arg) {
 Py_VISIT(GETSTATE(m)->error);
 return 0;
}

static int grmutils_clear(PyObject *m) {
 Py_CLEAR(GETSTATE(m)->error);
 return 0;
}

static struct PyModuleDef grmutilsModule = {
 PyModuleDef_HEAD_INIT,
 "grmutils",
 NULL,
 sizeof(struct module_state),
 grmutilsmethods,
 NULL,
 grmutils_traverse,
 grmutils_clear,
 NULL
};

#define INITERROR return NULL

PyObject* PyInit_grmutils(void)

#else
#define INITERROR return

void init_grmutils(void)
#endif
{
 #if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&grmutilsModule);
 #else
  PyObject *module = Py_InitModule("grmutils",grmutilsmethods);
 #endif

 if ( module == NULL )
  INITERROR;
 struct module_state *st = GETSTATE(module);
 st->error = PyErr_NewException("grmutils.Error",NULL,NULL);
 if ( st->error == NULL ) {
  Py_DECREF(module);
  INITERROR;
 }
 import_array();
 #if PY_MAJOR_VERSION >= 3
  return module;
 #endif
}

static PyObject *grm_loadDosage(PyObject *self, PyObject *args) {
 /***
  Load the dosage matrix into global memory
  Global variables: DOSAGES_PYMAT - will be a pointer to the PyArrayObject of the dosage array (must incref this)
                    DOSAGES - will be a pointer to the double array underlying DOSAGES_PYMAT
 ***/
 PyArrayObject *dosageMatrixArr;
 if ( ! PyArg_ParseTuple(args,"O",&DOSAGES_PYMAT_OBJ) ) return NULL;
 if ( NULL == DOSAGES_PYMAT_OBJ ) return NULL;
 Py_INCREF(DOSAGES_PYMAT_OBJ);
 
 /* from PyObject to Python Array */
 dosageMatrixArr = pymatrix(DOSAGES_PYMAT_OBJ);
 /* get the row and col sizes */
 N_VARIANTS = dosageMatrixArr->dimensions[0];
 N_SAMPLES = dosageMatrixArr->dimensions[1];
 DOSAGES = pymatrix_to_Carray(dosageMatrixArr);
 /** note: I'm not sure if you *have* to do this conversion. This may work fine
     with DOSAGES = (double*) dosageMatrixArr->data. Right now this allocates
     extra data AND ALSO INCREFs the overlying PyObject. I think one of these
     can be omitted to save memory. **/
 Py_RETURN_TRUE;
}

static PyObject *grm_calcdistance(PyObject *self, PyObject *args) {
 /** Given column indeces (samples) of DOSAGES, and an array of row indeces (the variants missing for one or both),
     calculate the genetic distance **/

 int samI,samJ,nMissing, *missing;
 PyObject *missingObj;
 PyArrayObject *missingArr;

 if ( ! PyArg_ParseTuple(args,"iiOi",&samI,&samJ,&missingObj,&nMissing) ) return NULL;
 if ( NULL == missingObj ) return NULL;
 
 missingArr = pyvector(missingObj);
 missing = pyvector_to_Carray(missingArr);

 double replaced1[nMissing];
 double replaced2[nMissing];

 int missingVectorIdx;
 int missingVariantIdx;
 for ( missingVectorIdx = 0; missingVectorIdx < nMissing; missing++ ) {
  missingVariantIdx = missing[missingVectorIdx];
  replaced1[missingVariantIdx] = DOSAGES[missingVariantIdx][samI];
  replaced2[missingVariantIdx] = DOSAGES[missingVariantIdx][samJ];
  DOSAGES[missingVariantIdx][samI]=0.0;
  DOSAGES[missingVariantIdx][samJ]=0.0;
 }

 double distance = _calcDistance(samI,samJ);

 for ( missingVectorIdx = 0; missingVectorIdx < nMissing; missing++ ) {
  missingVariantIdx = missing[missingVectorIdx];
  DOSAGES[missingVariantIdx][samI] = replaced1[missingVariantIdx];
  DOSAGES[missingVariantIdx][samJ] = replaced2[missingVariantIdx];
 }

 return PyFloat_FromDouble(distance);
 
}

static PyObject *grm_calcdistance_nomissing(PyObject *self, PyObject *args) {
 /** Given only column and row indeces of DOSAGES, calculate the genetic distance **/
 int samI,samJ;
 if ( ! PyArg_ParseTuple(args,"ii",&samI,&samJ) ) return NULL;

 return PyFloat_FromDouble(_calcDistance(samI,samJ));
}

static double _calcDistance(int i, int j) {
 /** Calculate the distance between samples i and j. DOSAGE (global) is a SNP x SAMPLE matrix **/
 double dist = 0.0;
 int variant;
 for ( variant = 0; variant < N_VARIANTS; variant++ ) {
  dist = dist + DOSAGES[variant][i]*DOSAGES[variant][j];
 }
 dist = dist/N_VARIANTS;
 return dist;
}

PyArrayObject *pymatrix(PyObject *objin)  {
  return (PyArrayObject *) PyArray_GETCONTIGUOUS(objin);
}

PyArrayObject *pyvector(PyObject *objin)  {
  return (PyArrayObject *) PyArray_GETCONTIGUOUS(objin);
}

double **pymatrix_to_Carray(PyArrayObject *arrayin)  {
 double **c, *a;
 int i,n,m;
 n=arrayin->dimensions[0];
 m=arrayin->dimensions[1];
 c=ptrvector(n);
 a=(double *) arrayin->data;  /* pointer to arrayin data as double */
 for ( i=0; i<n; i++)  {
  c[i]=a+i*m;
 }
 return c;
}

double **ptrvector(long n)  {
 double **v;
 v=(double **)malloc((size_t) (n*sizeof(double)));
 if (!v)   {
  printf("In **ptrvector. Allocation of memory for double array failed.");
  exit(0);
 }
 return v;
}

double *pyvector_to_Carray(PyArrayObject *arrayin)  {
 int i,n;
 n=arrayin->dimensions[0];
 return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

void free_Carrayptrs(double **v)  {
 free((char*) v);
}
