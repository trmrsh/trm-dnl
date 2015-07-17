//
// Python/C interface file for dnl programs
//

#include <Python.h>
#include "numpy/arrayobject.h"

#include <vector>
#include <map>
#include <iostream>
#include "trm/subs.h"
#include "trm/formula.h"
#include "trm/buffer2d.h"
#include <stdexcept>

// A simple exception class
class Dnl_Error {
public:
    
    //! Default constructor
    Dnl_Error() {}
    
    //! Constructor which sets an error message with the Python interface
    Dnl_Error(PyObject* etype, const std::string& str) {
	PyErr_SetString(etype, str.c_str());
    }

};

/* Extracts a header item from a Dset. Dset is a sub-type of dictionaries
 * and stores header items in that form. If in Python a header item is defined 
 * as dset['Telescope']['Longitude'] = 123.45 then you would supply
 * header = "Telescope.Longitude" to this routine and it will recursively 
 * look for it, returning NULL if it can't find what you want.
 * 
 * dset    -- the Dset to search through
 * header  -- string defining the header item name; see above.
 *
 * Returns a PyObject* which it will be up to you to translate
 */

PyObject* get_header_value(PyObject* dset, const std::string& header){

    if(dset == NULL) return NULL;

    PyObject* value = NULL;
    std::string::size_type pdot = header.find('.');
    if(pdot == std::string::npos){
	value = PyDict_GetItemString(dset, header.c_str());
    }else{
	PyObject* next = PyDict_GetItemString(dset, header.substr(0,pdot).c_str());
	value = get_header_value(next, header.substr(pdot+1));
    }
    return value;
}

/* Sets a header item in a Dset. Dset is a sub-type of dictionaries
 * and stores header items in that form. If in Python a header item is defined 
 * as dset['Telescope']['Longitude'] = 123.45 then you would supply
 * header = "Telescope.Longitude" to this routine and it will recursively 
 * look for it, returning NULL if it can't find what you want. The header
 * item must exist.
 * 
 * dset    -- the Dset to search through
 * header  -- string defining the header item name; see above.
 * value   -- value to change the header item to.
 *
 * Returns a bool which is false if things have gone wrong
 */

bool set_header_value(PyObject* dset, const std::string& header, double value){

    if(dset == NULL) return false;

    bool rval;
    std::string::size_type pdot = header.find('.');
    if(pdot == std::string::npos){
	if(PyDict_SetItemString(dset, header.c_str(), PyFloat_FromDouble(value)) == 0)
	    rval = true;
	else
	    rval = false;
    }else{
	PyObject* next = PyDict_GetItemString(dset, header.substr(0,pdot).c_str());
	rval = set_header_value(next, header.substr(pdot+1), value);
    }
    return rval;
}

/* Func is a function object that stores everything needed to be able to
 * compute a model and its derivatives at a point x. The model is defined from
 * a file. The model is allowed to have parameters that are independent from
 * slot to slot so there is rather complex code to allow relevant values to be
 * stored and retrieved. Most of the hard work is done in the constructor from the file.
 */

class Func {

public:
  
    /* Constructor of the function object from an algebraic expression, maps of names & values, and
     * a vector of which ones are to be varied. This is not quite complete as the object must be
     * subsequently initialised with a particular slot (because in general there is a potential need for
     * header values) using 'slot_initialise'
     *
     *  plist --- a list of Dsets 
     *  fname --- file name
     */

    Func(PyObject* plist, const char* fname) : ready(false) {
    
	// Kick off with a few checks
	if(!PyList_Check(plist))
	    throw Dnl_Error(PyExc_TypeError, "Func::Func: second argument must be a list");

	// nlist used later, number of Dsets in the list
	Py_ssize_t nlist = PyList_Size(plist);
	if(nlist == 0)
	    throw Dnl_Error(PyExc_ValueError, "Func::Func: list pointed at by second argument has no elements");
	
	// Now read and interpret the fit file
	std::ifstream ifstr(fname);
	if(!ifstr)
	    throw Dnl_Error(PyExc_IOError, std::string("Func::Func: file = ") + fname + " could not be opened.");

	std::string line, express;
	std::string::size_type posn;
	bool more = false;
	int var_index = 0;

	// loop through the fit file, line by line
	while(ifstr){
	    getline(ifstr, line);
	    if(line[0] != '#'){ // skip comments
		if(more){
		    if(line == ""){
			more = false;
		    }else{
			std::string::size_type pos2 = line.find('\\'); // continuation
			if(pos2 == std::string::npos){
			    express += line;
			    more = false;
			}else{
			    express += line.substr(0,pos2);
			}
		    }
		    
		}else if((posn = line.find('=')) != std::string::npos){
		    
		    // Read a value setting line. These can be of the following forms:
		    //
		    // height = 10           # intialises a variable value
		    // height = 10 f         # sets a constant value (can end in 'f' or 'F')
		    // height = fit.height   # a separate variable from each slot initialised and stored from headers
		    // height = fit.height f # a constant grabbed from the headers of each slot.
		    
		    std::string::size_type posb = line.find(' ');
		    std::string name = line.substr(0,std::min(posb,posn));
		    
		    if(name == "Chisq"  || name == "Ndof"   || name == "PI"    || name == "TWOPI"  ||
		       name == "UNIT"   || name == "ZERO"   || name == "MUNIT" || name == "DAY"    ||
		       name == "VLIGHT" || name == "HALPHA" || name == "HBETA" || name == "HGAMMA" ||
		       name == "HDELTA")
			throw Dnl_Error(PyExc_IOError, "Func::Func: variable name = " + name + " has a special meaning and is not permitted.");
		    
		    std::istringstream istr(line.substr(posn+1));
		    double value;
		    istr >> value;
		    if(!istr){
			
			// Failed to read value as a double; proceed on the basis that it must be a header. 
			// Find start and end of header
			std::string::size_type start = line.find_first_not_of(" \t", posn+1);
			if(start == std::string::npos)
			    throw Dnl_Error(PyExc_IOError, "Func::Func: failed to translate line = " + line + " as a header item or a number.");

			std::string::size_type end = line.find_first_of(" \t\n", start);
			if(end == std::string::npos)
			    end = line.length();
			std::string hitem = line.substr(start, end-start);
		    
			// store the variable
			vals[name]        = 0.;
			this->var_to_head[name] = hitem;
		    
			// Get first Dset from the list
			PyObject* dset = NULL;
			
			if(end != line.length() && (line.find('s',end) != std::string::npos || line.find('S',end) != std::string::npos)){
			    
			    // get first item of list
			    dset = PyList_GetItem(plist, 0);
			    
			    // retrieve the value of the header item
			    PyObject* dptr = get_header_value(dset, hitem);
			    if(!dptr)
				throw Dnl_Error(PyExc_IOError, "Func::Func: could not find header item = " + hitem + " in the first Dset supplied");

			    double value = PyFloat_AsDouble(dptr);

			    if (PyErr_Occurred())
				throw Dnl_Error(PyExc_IOError, "Func::Func: could not interpret header item = " + hitem + " as a double.");
			
			    // store
			    vals[name] = value;
			    this->is_head.push_back(false);
			    
			    // is it variable or fixed?
			    if(line.find('f',end) == std::string::npos && line.find('F',end) == std::string::npos){
				this->vars.push_back(name);
				var_index++;
			    }
			    
			}else if(end != line.length() && (line.find('f',end) != std::string::npos || line.find('F',end) != std::string::npos)){
			    
			    // Values for model value = name taken from each slot independently, but kept fixed
			    this->hconst.push_back(name);
			    
			    // Now check that every spectrum has this header item.
			    for(Py_ssize_t i=0; i<nlist; i++){
				dset = PyList_GetItem(plist, i);
				if(!dset)
				    throw Dnl_Error(PyExc_TypeError, "Func::Func: second input is not a list");

				PyObject* dptr = get_header_value(dset, hitem);
				if(!dptr)
				    throw Dnl_Error(PyExc_ValueError, "Func::Func: could not find header item = " 
						    + hitem + " in the " + Subs::str(i+1) + "th Dset supplied");
			    } 
			
			}else{
			
			    for(Py_ssize_t i=0; i<nlist; i++){
				dset = PyList_GetItem(plist, i);
				if(!dset)
				    throw Dnl_Error(PyExc_TypeError, "Func::Func: second input is not a list");

				PyObject* dptr = get_header_value(dset, hitem);
				if(!dptr)
				    throw Dnl_Error(PyExc_ValueError, "Func::Func: could not find header item = " 
						    + hitem + " in the " + Subs::str(i+1) + "th Dset supplied");

				double value = PyFloat_AsDouble(dptr);

				if (PyErr_Occurred())
				    throw Dnl_Error(PyExc_ValueError, "Func::Func: could not interpret header item = " 
						    + hitem + " of the " + Subs::str(i+1) + "th Dset as a double.");

				this->hvals[name].push_back(value); // initialise value
				this->vars.push_back(name);
				this->is_head.push_back(true);
				this->var_ind_to_slot[var_index++] = i; 
			    } 
			}
			
		    }else{
			
			// Single value, initialised from file 
			vals[name] = value;
			this->is_head.push_back(false);
			
			// fixed or not?
			std::string::size_type posf;;
			if(!(
			       ((posf = line.find('f')) != std::string::npos && posf > posn) ||
			       ((posf = line.find('F')) != std::string::npos && posf > posn))){
			    this->vars.push_back(name);
			    var_index++;
			}
		    }
		    
		}else if((posn = line.find("Equation:")) != std::string::npos){
		    
		    // The equation reading line
		    std::string::size_type pos2 = line.find('\\');
		    if(pos2 == std::string::npos){
			express = line.substr(posn+9);
		    }else{
			express = line.substr(posn+9,pos2-posn-9);
			more = true;
		    }
		}
	    }
	}
	ifstr.close();
	
//	this->list_vars(std::cerr);

	// Now create the formula
	try{
	    this->form = Formula::Formula(express);
	    this->vals["x"] = 0.;
	    this->form.check(this->vals);
	}
	catch(const std::runtime_error& err){
	    throw Dnl_Error(PyExc_ValueError, std::string("Func::Func: ") + err.what());
	}   

	// Create derivative formulae
	std::string previous;
	for(size_t i=0; i<this->vars.size(); i++){
	    if(i == 0 || this->vars[i] != previous){
		this->derv[this->vars[i]] = this->form.deriv(this->vars[i]);
	    }
	    previous = this->vars[i];
	}
	this->can_vary.resize(this->vars.size());
    }
  
    // Returns number of variable parameters
    int nvar() const {return this->vars.size();}

    // Increments model values
    void add_vars(const Subs::Buffer1D<double>& dvals){
	for(size_t i=0; i<this->vars.size(); i++){
	    if(this->is_head[i]){
		this->hvals[this->vars[i]][this->var_ind_to_slot[i]] += dvals[i];
	    }else{
		vals[this->vars[i]] += dvals[i];
	    }
	}
    }
  
    // Decrements model values
    void sub_vars(const Subs::Buffer1D<double>& dvals){
	for(size_t i=0; i<this->vars.size(); i++){
	    if(this->is_head[i]){
		this->hvals[this->vars[i]][this->var_ind_to_slot[i]] += dvals[i];
	    }else{
		vals[this->vars[i]] += dvals[i];
	    }
	}
    }
  
    // Initialise model values from a given slot's Dset and its headers
    void initialise(PyObject* dlist, int slot){

	// get the dset
	PyObject* dset = PyList_GetItem(dlist, slot);
	if(!dset)
	    throw Dnl_Error(PyExc_ValueError, "initialise: failed to retrieve Dset " + Subs::str(slot));

	for(size_t i=0; i<this->can_vary.size(); i++){
	    if(!this->is_head[i] || this->var_ind_to_slot[i] == slot){
		this->can_vary[i] = true;
	    }else{
		this->can_vary[i] = false;
	    }
	}
	
	// first, the internal variables associated with headers that are being optimised
	for(std::map<std::string, std::vector<double> >::iterator cit=this->hvals.begin(); cit!=this->hvals.end(); cit++)
	    vals[cit->first] = cit->second[slot];
	
	// second, the values associated with headers that are constant
	for(size_t i=0; i<this->hconst.size(); i++){
	    PyObject* dptr = get_header_value(dset, this->var_to_head[this->hconst[i]]);
	    if(!dptr)
		throw Dnl_Error(PyExc_ValueError, "Func::initialise: could not find header item = " 
				+ this->var_to_head[this->hconst[i]] + " in the Dset supplied");

	    double value = PyFloat_AsDouble(dptr);
	    if (PyErr_Occurred())
		throw Dnl_Error(PyExc_ValueError, "Func::initialise: could not interpret header item = " 
				+ this->var_to_head[this->hconst[i]] + " as a double.");

	    vals[this->hconst[i]] = value;
	}
	    
	this->ready = true;
    }

    // Computes model at x 
    double operator()(double x) {
	if(!this->ready)
	    throw Dnl_Error(PyExc_ValueError, "Func::operator(): Func object not yet initialised");
	vals["x"] = x;
	return this->form.value(vals);
    }   
  
    // Computes derivatives at x
    void deriv(double x, Subs::Buffer1D<double>& d) {
	if(!this->ready)
	    throw Dnl_Error(PyExc_ValueError, "Func::operator(): Func object not yet initialised");
	vals["x"] = x;
	for(size_t i=0; i<this->vars.size(); i++){
	    if(this->can_vary[i])
		d[i] = this->derv[this->vars[i]].value(vals);
	    else
		d[i] = 0.;
	}
    }
    
    // Substitutes variables with numbers
    void simplify(){
	this->vals.erase("x");
	this->form.subst(this->vals);
	for(std::map<std::string, Formula::Formula>::iterator cit=this->derv.begin(); cit != this->derv.end(); cit++)
	    cit->second.subst(this->vals);
    }
    
    // Gets value of variable with index number index
    double get_value(int index) {
	if(this->is_head[index])
	    return this->hvals[this->vars[index]][this->var_ind_to_slot[index]];
	else
	    return vals[this->vars[index]];
    }

    // Gets name of variable with index number index
    std::string get_name(int index) const {
	return this->vars[index];
    }

    // Is the variable obtained from the headers
    bool is_from_header(int index) const {
	return this->is_head[index];
    }

    // Gets header name equivalent to a given index (if any)
    std::string get_header_name(int index) {
	if(this->is_head[index])
	    return var_to_head[vars[index]];
	else
	    return std::string();
    }

    // Output names, types and values of all variables
    void list_vars(std::ostream& ostr) {
	for(int i=0; i<this->nvar(); i++){
	    ostr << "Variable = " << vars[i];
	    if(is_head[i])
		ostr << ", from headers of slot " << this->var_ind_to_slot[i];
	    else
		ostr << ", from model file";
	    ostr << ", value = " << get_value(i) << std::endl;
	}	    
    }

private:

    std::vector<std::string> vars;                       // name of model parameter associated with each variable
    std::vector<bool> is_head;                           // true/false for header variable or not
    std::map<std::string, std::vector<double> > hvals;   // values equivalent to header variables for each slot 
    std::map<std::string, double> vals;                  // model value for model name
    std::map<std::string, std::string> var_to_head;      // gives header name corresponding to a given a variable 
    std::map<size_t, int> var_ind_to_slot;               // given variable index of header variables, returns equivalent slot 
    std::vector<std::string> hconst;                     // model constant names to be set from headers for each slot
    std::vector<bool> can_vary;                          // true for each variable that has a potentially non-zero derivative
    bool ready;                                          // indicates that we can go.
    Formula::Formula form;                               // main formula
    std::map<std::string, Formula::Formula> derv;        // derivatives for each variable of the equation
    
};

/* Structure to handle the many pointers required to access a Dset 
 */

struct Dset {

    /* Constructor which gets pointers to the data in a dset.
     * dset  -- the dset to unravel
     */
    
    Dset(PyObject* dset){
	
	this->xaxis = PyObject_GetAttrString(dset, "x");
	if(!this->xaxis)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not retrieve X-axis from Dset");

	this->xdat = PyObject_GetAttrString(this->xaxis, "dat");
	if(!this->xdat)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not retrieve dat from X-axis of Dset");

	this->xobj = PyArray_FROM_OTF(this->xdat, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
	if(!this->xobj)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not convert x data of Dset to a double array");

	this->x = (double*) PyArray_DATA(this->xobj);
	
	this->yaxis = PyObject_GetAttrString(dset, "y");
	if(!this->yaxis)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not retrieve Y-axis from Dset");
	
	this->ydat = PyObject_GetAttrString(this->yaxis, "dat");
	if(!this->ydat)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*) :could not retrieve dat from Y-axis of Dset");
	
	this->ydobj = PyArray_FROM_OTF(this->ydat, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
	if(!this->ydobj)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not convert y data of Dset to a double array");
	
	this->y = (double*) PyArray_DATA(this->ydobj);
	
	this->yerr = PyObject_GetAttrString(this->yaxis, "err");
	if(!this->yerr)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not retrieve err from Y-axis of Dset");
	
	this->yeobj = PyArray_FROM_OTF(this->yerr, NPY_FLOAT, NPY_IN_ARRAY | NPY_FORCECAST);
	if(!this->yeobj)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not convert y errors of Dset to a float array");
	
	this->ye = (float*) PyArray_DATA(this->yeobj);
	this->npix = PyArray_Size(this->xdat);
	
	// data masks
	this->mask = PyObject_GetAttrString(dset, "mask");
	if(!this->mask)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not retrieve mask from Dset");
	
	this->mobj = PyArray_FROM_OTF(this->mask, NPY_BOOL, NPY_IN_ARRAY | NPY_FORCECAST);
	if(!this->mobj)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not convert mask of Dset to a boolean array");
	
	this->m = (bool*) PyArray_DATA(this->mobj);
	
	this->good = PyObject_GetAttrString(dset, "good");
	if(!this->good)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not retrieve good from Dset");
	
	this->gobj = PyArray_FROM_OTF(this->good, NPY_BOOL, NPY_IN_ARRAY | NPY_FORCECAST);
	if(!this->gobj)
	    throw Dnl_Error(PyExc_ValueError, "Dset(PyObject*): could not convert good of Dset to a boolean array");
	
	this->g = (bool*) PyArray_DATA(this->gobj);
    }
    
    // Destructor XDECREFs all PyObjects
    ~Dset(){
	Py_XDECREF(xaxis);
	Py_XDECREF(xdat);
	Py_XDECREF(xobj);
	Py_XDECREF(yaxis);
	Py_XDECREF(ydat);
	Py_XDECREF(yerr);
	Py_XDECREF(ydobj);
	Py_XDECREF(yeobj);
	Py_XDECREF(mask);
	Py_XDECREF(mobj);
	Py_XDECREF(good);
	Py_XDECREF(gobj);
    }

    // PyObjects that need DECREF-ing 
    PyObject* xaxis;
    PyObject* xdat; 
    PyObject* xobj; 
    PyObject* yaxis;
    PyObject* ydat; 
    PyObject* yerr; 
    PyObject* ydobj;
    PyObject* yeobj;
    PyObject* mask;
    PyObject* mobj;
    PyObject* good;
    PyObject* gobj;

    // Number of pixels and data pointers
    npy_intp npix;
    double* x;
    double* y;
    float* ye;
    bool *m, *g;
};
    

/* This function evaluates the fit coresponding to one dset
 * 
 * dlist -- list of Dsets
 * slot  -- the particular one to evaluate the fit for
 * func  -- function object
 * fit   -- pointer where the fit will be stored. Must have had space allocated already.
 */

void fit_func_eval(PyObject* dlist, int slot, const Func& func, double* fit){

    // get the dset
    PyObject* dset = PyList_GetItem(dlist, slot);
    if(!dset)
	throw Dnl_Error(PyExc_ValueError, "fit_func_eval: failed to retrieve Dset " + Subs::str(slot));
    
    // Get the data pointers
    Dset ddat(dset);
	
    // Copy the function so that it can be simplified to speed up evaluation
    Func tfunc(func);
    
    // Set any header item values
    tfunc.initialise(dlist, slot);
    
    // Simplify
    tfunc.simplify();  
    
    // Calculate 
    for(npy_intp n=0; n<ddat.npix; n++)
	fit[n] = tfunc(ddat.x[n]);

}

/* This function carries out pixel rejection
 */

int fit_func_rej(PyObject* dlist, const Func& func, double thresh){

    Py_ssize_t nlist = PyList_Size(dlist);
    if (PyErr_Occurred())
	throw Dnl_Error(PyExc_ValueError, "fit_func_rej: could not find size of Dset list");
    
    if (nlist == 0)
	throw Dnl_Error(PyExc_ValueError, "fit_func_rej: empty list entered");

    int nrej = 0;
    for(Py_ssize_t ns=0; ns<nlist; ns++){

	// get the dset
	PyObject* dset = PyList_GetItem(dlist, ns);
	if(!dset)
	    throw Dnl_Error(PyExc_ValueError, "fit_func_rej: failed to retrieve Dset " + Subs::str(ns));
	
	// Get the data pointers
	Dset ddat(dset);

	// Copy the function so that it can be simplified to speed up evaluation
	Func tfunc(func);

	// Set any header item values
	tfunc.initialise(dlist, ns);
	
	// Simplify
	tfunc.simplify();  

	// finally, get on and compute
	for(npy_intp n=0; n<ddat.npix; n++){
	    
	    if(ddat.m[n] && ddat.g[n]){
		double dev = abs((ddat.y[n] - tfunc(ddat.x[n]))/ddat.ye[n]);
		if(dev > thresh){
		    ddat.m[n] = false;
		    nrej++;
		}
	    }
	}
    }
    return nrej;
}

/* Uses the Levenburg-Marquardt method to fit a general function 
 *
 * dlist   --- list of Dsets to fit
 * func    --- function object
 * chisq   --- Chi**2, returned
 * lambda  --- multiplier, < 0 to initialise
 */
    
void fit_func_cof(PyObject* dlist, const Func& func, Subs::Buffer2D<double>& alpha, Subs::Buffer1D<double>& beta, double& chisq, int& ndata);
    
void fit_func(PyObject* dlist, Func& func, double& chisq, double& lambda, int& ndata, Subs::Buffer2D<double>& covar){

    static Subs::Buffer1D<double> beta, da;
    static Subs::Buffer2D<double> oneda, alpha;
    static double ochisq;

    if(lambda < 0.){
	beta.resize(func.nvar());
	da.resize(func.nvar());
	oneda.resize(func.nvar(),1);
	alpha.resize(func.nvar(),func.nvar());
	covar.resize(func.nvar(),func.nvar());
	lambda=0.001;
	fit_func_cof(dlist, func, alpha, beta, chisq, ndata);
	ochisq = chisq;
    }

    // Alter linearised fitting matrix by augmenting diagonal elements
    for (int j=0; j<func.nvar(); j++){
	for (int k=0; k<func.nvar(); k++) covar[j][k]=alpha[j][k];
	covar[j][j]=alpha[j][j]*(1.0+lambda);
	oneda[j][0]=beta[j];
    }

    // Matrix solution to find a (hopefully) better solution
    try {
	Subs::gaussj(covar, oneda);
    }
    catch(...){
	throw Dnl_Error(PyExc_ValueError, "fit_func: singular matrix in gaussj, probably because of degeneracy");
    }

    for (int j=0; j<func.nvar(); j++) da[j] = oneda[j][0];

    // lambda == 0 indicates convergence
    if(lambda == 0.) return;

    Subs::Buffer1D<double> save(da);
    func.add_vars(da);
    fit_func_cof(dlist, func, covar, da, chisq, ndata);

    // If Chi**2 reduced ==> success & reduce lambda, else increase lambda
    if(chisq < ochisq){
	lambda *= 0.1;
	ochisq  = chisq;
	alpha   = covar;
	beta    = da;
    }else{
	lambda *= 10.0;
	func.sub_vars(save);
	chisq   = ochisq;
    }
}

/* Routine needed by the Levenburg-Marquardt fit_func routine
 *
 * dlist   --- list of Dsets to fit
 * func    --- function object
 * alpha   --- 2D matrix used in Lev.Marq. method
 * beta    --- 1D vector used in Lev.Marq. method
 * chisq   --- Chi**2, returned
 */

void fit_func_cof(PyObject* dlist, const Func& func, Subs::Buffer2D<double>& alpha, Subs::Buffer1D<double>& beta, double& chisq, int& ndata){

    // for derivatives
    Subs::Buffer1D<double> dyda(func.nvar());
  
    // Initialise alpha and beta
    alpha = 0.;
    beta  = 0.;

    // Now the real work.
    chisq = 0.;
    double wt, wgt;

    Py_ssize_t nlist = PyList_Size(dlist);
    if (PyErr_Occurred())
	throw Dnl_Error(PyExc_ValueError, "fit_func_rej: could not find size of Dset list");

    if (nlist == 0)
	throw Dnl_Error(PyExc_ValueError, "fit_func_rej: empty list entered");

    ndata = 0;
    for(Py_ssize_t ns=0; ns<nlist; ns++){

	// get the dset
	PyObject* dset = PyList_GetItem(dlist, ns);
	if(!dset)
	    throw Dnl_Error(PyExc_ValueError, "fit_func_cof: failed to retrieve Dset " + Subs::str(ns));

	// Get the data pointers
	Dset ddat(dset);

	// Copy the function so that it can be simplified to speed up evaluation
	Func tfunc(func);

	// Set any header item values
	tfunc.initialise(dlist, ns);
	
	// Simplify
	tfunc.simplify();  

	// finally, get on and compute
	for(npy_intp n=0; n<ddat.npix; n++){

	    if(ddat.m[n] && ddat.g[n]){
		double diff = ddat.y[n] - tfunc(ddat.x[n]);
		tfunc.deriv(ddat.x[n], dyda);
		wgt = 1./Subs::sqr(ddat.ye[n]);
	    
		for(int j=0, l=0; l<tfunc.nvar(); l++){
		    wt = wgt*dyda[l];
		    for(int k=0, m=0; m<=l; m++)
			alpha[j][k++] += wt*dyda[m];
		    beta[j++]  += wt*diff;
		}
		chisq += wgt*Subs::sqr(diff);
		ndata++;
	    }
	}
    }
  
    // symmetrise
    for(int j=1; j<func.nvar(); j++)
	for(int k=0; k<j; k++) alpha[k][j] = alpha[j][k];

}

// The function called by Python
// Arguments:
//  1) List of Dsets to fit
//  2) Name of file containing model
//  3) Rejection threshold, sigmas
//  4) Number of iterations to carry out per reject cycle

static PyObject* 
ffit(PyObject *self, PyObject *args)
{
    PyObject* dlist;
    char* fname;
    int nits = -1, output=0, ndelta=2;
    double thresh, delta = 1.e-3;
    if(!PyArg_ParseTuple(args, "Osd|iid", &dlist, &fname, &thresh, &output, &nits, &delta, &ndelta))
	return NULL;

    // A few checks

    if(delta <= 0.){
	PyErr_SetString(PyExc_ValueError, "dnl.ffit: delta <= 0.");
	return NULL;
    }
    if(thresh <= 1.){
	PyErr_SetString(PyExc_ValueError, "dnl.thresh: thresh <= 1.");
	return NULL;
    }
    if(ndelta <= 0){
	PyErr_SetString(PyExc_ValueError, "dnl.ffit: ndelta <= 0");
	return NULL;
    }

    PyObject* vtup  = NULL;
    PyObject* vec   = NULL;
    PyObject* vdict = NULL;
    PyObject* cdict = NULL;

    try { 

	// Create the function object
	Func func(dlist, fname);

	Subs::Buffer2D<double> covar;

	double chisq=0., cold=1000.;
	int ndata, nrej = 1, nrejc = 0, nrejtot = 0, nitt=0, nsame = 0;
	while(nrej > 0){
	    nrej = 0;

	    double lambda=-1.;
	    int ncount = 0;
	    while((nits == -1 || ncount < nits) && (chisq > cold || chisq < cold - delta || nsame < ndelta)){
		fit_func(dlist, func, chisq, lambda, ndata, covar);
		if(ncount == 0) cold = chisq + 2.;
		if(chisq >= cold - delta && chisq <= cold){
		    nsame++;
		}else{
		    nsame = 0;
		}
		cold = chisq;
		nitt++;
		ncount++;
		if(output)
		    std::cerr << "Cycle " << ncount << ", chi**2 = " << chisq << ", lambda = " << lambda << std::endl;
	    }

	    // Carry out rejection; should probably subtract number of variables from ndata in next line
	    nrej = fit_func_rej(dlist, func, thresh*sqrt(chisq/ndata));
	    nrejtot += nrej;
	    nrejc++;
	    if(output)
		std::cerr << "Reject cycle " << nrejc << " rejected " << nrej << " pixels" << std::endl;
	}

	// Compute covariances
	double lambda = 0.;
	fit_func(dlist, func, chisq, lambda, ndata, covar);

	if(output)
	    func.list_vars(std::cerr);

	// Create dictionaries for storage of variable names and values and covariances 
	// for return to the calling program
	Py_ssize_t nvar = func.nvar();
	if(nvar){

	    // names and values
	    vdict = PyDict_New();
	    if(!vdict)
		throw Dnl_Error(PyExc_MemoryError, "ffit: failed to create dictionary of variable values and names.");
	    for(Py_ssize_t i=0; i<nvar; i++)
		PyDict_SetItemString(vdict, func.get_name(i).c_str(), PyFloat_FromDouble(func.get_value(i)));

	    // covariance dictionaries
	    cdict = PyDict_New();
	    if(!cdict)
		throw Dnl_Error(PyExc_MemoryError, "dnl.ffit: failed to create dictionary of covariances.");
	    for(Py_ssize_t i=0; i<nvar; i++){
		PyObject* adict = PyDict_New();
		if(!adict)
		    throw Dnl_Error(PyExc_MemoryError, "dnl.ffit: failed to create sub-dictionary of covariances.");
		for(Py_ssize_t j=0; j<nvar; j++)
		    PyDict_SetItemString(adict, func.get_name(j).c_str(), PyFloat_FromDouble(covar[int(i)][int(j)]));
		PyDict_SetItemString(cdict, func.get_name(i).c_str(), adict);
	    }
	}

	// Create tuple for storage of output fits
	Py_ssize_t ndset = PyList_GET_SIZE(dlist);

	vtup = PyTuple_New(ndset);
	if(!vtup)
	    throw Dnl_Error(PyExc_MemoryError, "ffit: failed to create tuple of fits.");

	for(Py_ssize_t i=0; i<ndset; i++){
	    
	    PyObject* dset = PyList_GetItem(dlist, i);
	    if(!dset)
		throw Dnl_Error(PyExc_ValueError, "ffit: failed to retrieve Dset " + Subs::str(i));

	    // Get the data size
	    Dset ddat(dset);

	    npy_intp dims[1] = {ddat.npix};
	    vec = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	    if(!vec)
		throw Dnl_Error(PyExc_MemoryError, "ffit: failed to create fit vector for Dset " + Subs::str(i));

	    PyArrayObject* farr = (PyArrayObject*)vec;
	    double* fit = (double*)farr->data;
	    fit_func_eval(dlist, i, func, fit);
	    PyTuple_SetItem(vtup, i, vec); 

	    // Update the values of any variables taken from the headers
	    for(int j=0; j<nvar; j++){
		if(func.is_from_header(j)){
		    if(!set_header_value(dset, func.get_header_name(j), func.get_value(j)))
			throw Dnl_Error(PyExc_RuntimeError, "ffit: failed to set value of header item = " + 
					func.get_header_name(j) + " for Dset " + Subs::str(i));
		}
	    }
	}

	// Return point
	return Py_BuildValue("OOOdiii",vtup,vdict,cdict,chisq,ndata,nrejtot,nitt);
    }
    catch(const Dnl_Error& err){
	Py_XDECREF(vtup);
	Py_XDECREF(vec);
	Py_XDECREF(vdict);
	// some code needed here to deal with sub-cdict dictionaries.
	Py_XDECREF(cdict);
	return NULL;
    }
    catch(const Subs::Subs_Error& err){
	std::string message = "ffit: " + err;
	PyErr_SetString(PyExc_TypeError, message.c_str());
	Py_XDECREF(vtup);
	Py_XDECREF(vec);
	Py_XDECREF(cdict);
	return NULL;
    }
};

//----------------------------------------------------------------------------------------
// The methods

static PyMethodDef DanalMethods[] = {

    {"ffit", ffit, METH_VARARGS, 
     "(ftup,vdict,cdict,chisq,ndata,nrtot,nitt) = ffit(dlist, fname, thresh, output=0, nits=-1, delta=1.e-3, ndelta=2)\n\n"
     "Fits a list of Dsets\n\n"
     "Arguments:\n"
     "dlist   -- list of Dsets to fit\n"
     "fname   -- name of file defining the fit\n"
     "thresh  -- rejection threshold, sigmas\n"
     "output  -- set true for output\n"
     "nits    -- number of Levenberg-Marquardt iterations per reject cycle, -1 for continue until convergence\n"
     "delta   -- minimum decrease in chi**2 needed\n\n"
     "ndelta  -- number of times minimum decrease in chi**2 must fail to be achieved to stop iterations\n\n"
     "Returns:\n"
     "ftup    -- tuple of fits to each Dset, each fit is a numpy 1D array\n"
     "vdict   -- dictionary of variable names and values\n"
     "cdict   -- dictionary of covariances\n"
     "chisq   -- chi**2 of final fit\n"
     "ndata   -- number of data points fitted, excluding rejects\n"
     "nrtot   -- number of data points rejected\n"
     "nitt    -- total number of iterations\n\n"
     "Format of model input file:\n"
     "An example model file is as follows:\n\n"
     "Equation: a + b*x + h*exp(-sqr((x-x0)/sig)/2)\n\n"
     "a   = 1.0\n"
     "b   = 0.1\n"
     "h   = height\n"
     "x0  = 101 f\n"
     "sig = 10\n\n"
     "This one fits a straight line plus an gaussian to data. Initial values are provided for the variables\n"
     "and those which are to kept fixed are indicated with an 'f'. Initial values can be obtained from header\n"
     "header items as in this case the height of the gaussian will be. This allows a different initial value\n"
     "for each spectrum. Note the equation line can be continued using the standard '\' continuation character.\n"
    },
    {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
init_dnl(void)
{
    (void) Py_InitModule("_dnl", DanalMethods);
    import_array();
}
