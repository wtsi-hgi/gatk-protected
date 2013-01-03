#include "Python.h"
#include "arrayobject.h"
#include <math.h>
#include <string.h>
#include "logitPredict.h"
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

static PyMethodDef logitPredictMethods[] = {
 {"error_out",(PyCFunction)error_out,METH_NOARGS,NULL},
 {"predict",lg_predict,METH_VARARGS,NULL},
 {"predictorsTimesWeights",lg_calcXW,METH_VARARGS,NULL},
 {NULL,NULL}
};

#if PY_MAJOR_VERSION >= 3

static int logitPredict_traverse(PyObject *m,visitproc visit, void *arg) {
 Py_VISIT(GETSTATE(m)->error);
 return 0;
}

static int logitPredict_clear(PyObject *m) {
 Py_CLEAR(GETSTATE(m)->error);
 return 0;
}

static struct PyModuleDef logitPredictModule = {
 PyModuleDef_HEAD_INIT,
 "logitPredict",
 NULL,
 sizeof(struct module_state),
 logitPredictMethods,
 NULL,
 logitPredict_traverse,
 logitPredict_clear,
 NULL
};

#define INITERROR return NULL

PyObject* PyInit_logitPredict(void)

#else
#define INITERROR return

void init_logitPredict(void)
#endif
{
 #if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&logitPredictModule);
 #else
  PyObject *module = Py_InitModule("logitPredict",logitPredictMethods);
 #endif

 if ( module == NULL )
  INITERROR;
 struct module_state *st = GETSTATE(module);
 st->error = PyErr_NewException("logitPredict.Error",NULL,NULL);
 if ( st->error == NULL ) {
  Py_DECREF(module);
  INITERROR;
 }
 import_array();
 #if PY_MAJOR_VERSION >= 3
  return module;
 #endif
}

static PyObject *lg_predict(PyObject *self, PyObject *args) {
 PyObject *predictorMatrixObj, *betaVecObj;
 PyArrayObject *predictorMatrix, *betaVec;
 npy_double** predictors;
 npy_double* beta;
 npy_intp nRow,nCol;
 npy_int N;
 if ( ! PyArg_ParseTuple(args,"OOi",&predictorMatrixObj,&betaVecObj,&N) ) return NULL;
 if ( NULL == predictorMatrixObj ) return NULL;
 if ( NULL == betaVecObj ) return NULL;

 Py_INCREF(predictorMatrixObj);
 Py_INCREF(betaVecObj);
 
 /* from PyObject to Python Array */
 predictorMatrix = pymatrix(predictorMatrixObj);
 betaVec = pymatrix(betaVecObj);

 /* get the row and col sizes */
 nRow = predictorMatrix->dimensions[0];
 nCol = predictorMatrix->dimensions[1];
 if ( betaVec->dimensions[0] != nCol || betaVec->dimensions[1] != 1) {
  printf("Beta dimensions are incorrect.\n");
  return NULL;
 }
 /* now to contiguous c arrays */
 predictors = pymatrix_to_Carray(predictorMatrix);
 beta = pyvector_to_Carray(betaVec);

 /* generate the predictor matrix. Predict is N x K, beta is K x 1, calculate compwise.logit(Predict*beta) */
 PyArrayObject *predictionMat = PyArray_ZEROS(1,&nRow,NPY_DOUBLE,0);
 npy_double* prediction = (npy_double*) predictionMat->data;
 npy_double b;
 npy_intp c,r;
 for ( c = 0; c < nCol; c++ ) {
  b = beta[c];
  for (r = 0; r < nRow; r++ ) {
   prediction[r] = prediction[r] + b*predictors[r][c];
  }
 }
 /* convert product into logit into predition */
 npy_double t;
 double eps = 0.0001;
 double oneMinusEps = 0.9999;
 double ulim = 10.0;
 double llim = -10.0;
 for ( r = 0; r < nRow; r++ ) {
  if ( prediction[r] > ulim ) {
   prediction[r] = N*oneMinusEps;
  } else if ( prediction[r] < llim ) {
   prediction[r] = N*eps;
  } else {
   t = exp(prediction[r]);
   prediction[r] = (N*t)/(1+t);
  }
 }
 
 /* remember to DECREF the objects */
 Py_DECREF(predictorMatrixObj);
 Py_DECREF(betaVec);

 /* free objects */
 free((char*) predictors);

 return PyArray_Return(predictionMat);
}

static PyObject *lg_calcXW(PyObject *self, PyObject *args) {
 /**** Given the predictors (N x K matrix) and predictions (N x 1 vector), calculate the weights T(Predictors)*diag(predictions) ****/
 PyObject *predictorObj, *predictionObj;
 PyArrayObject *predictorMat, *predictionVec;
 double** predictMatrix;
 double* predictions;
 npy_intp nRow,nCol, N,dims[2];

 if ( ! PyArg_ParseTuple(args,"OOi",&predictorObj,&predictionObj,&N) ) return NULL;
 if ( NULL == predictorObj ) return NULL;
 if ( NULL == predictionObj ) return NULL;

 Py_INCREF(predictorObj);
 Py_INCREF(predictionObj);

 /* to python array */
 predictorMat = pymatrix(predictorObj);
 predictionVec = pyvector(predictionObj);

 /* row and col */
 nRow = predictorMat->dimensions[0];
 nCol = predictorMat->dimensions[1];
 dims[0] = nCol;
 dims[1] = nRow;

 if ( predictionVec->dimensions[0] != nRow ) { 
  printf("Prediction vector dimensions are wrong\n");
  return NULL;
 }

 /* to contiguous c arrays */
 predictMatrix = pymatrix_to_Carray(predictorMat);
 predictions = pyvector_to_Carray(predictionVec);

 /* now we bave a NxK matrix predictMatrix (X) and Nx1 vector (W). Want to calculate X.T*diag(W), which will be KxN */
 PyArrayObject *XW_py = PyArray_SimpleNew(2,&dims,NPY_DOUBLE);
 npy_double** XW = pymatrix_to_Carray(XW_py); 

 npy_intp r,c;
 npy_double W_npq;
 double N_doub,val;

 N_doub = (double)((int) N);

 for ( r = 0; r < nRow; r++ ) {
  W_npq = N_doub*predictions[r]*(1.0-predictions[r]);
  for ( c = 0; c < nCol; c++ ) {
   val = W_npq*predictMatrix[r][c];
   XW[c][r] = W_npq*predictMatrix[r][c];
  }
 }

 /* remember to DECREF */
 Py_DECREF(predictorObj);
 Py_DECREF(predictionObj);

 /* free objects */
 free((char*) predictMatrix);
 free((char*) XW);

 return PyArray_Return(XW_py); 
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
