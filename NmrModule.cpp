///////////////////////////////////////////////////////////////////////////////////////////////
//
// NmrModule.cpp
// Provides functionality for simulation of cross relaxation during
// ROESY spin-locks, see comments and references below.
//
///////////////////////////////////////////////////////////////////////////////////////////////


#include "Python.h"
#include <boost/numeric/odeint.hpp>

#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include "numpy/arrayobject.h"

using namespace std;
using namespace boost::numeric::odeint;

// type for magnetizations and derivatives
typedef boost::array< double , 6 > state_type;



///////////////////////////////////////////////////////////////////////////////////////////////

class Ode_sys
{
    
    /* Class implementing the ODEs as a function object/functor (see ref. 1) */
    
public:
    
    // Constructor, with member initializaiton list (ref. 2) setting all the constants
    Ode_sys( double r_w_a, double r_R1a, double r_R2a, double r_w_b, double r_R1b, double r_R2b, double r_w1, double r_s, double r_u):
    w_a(r_w_a), R1a(r_R1a), R2a(r_R2a), w_b(r_w_b), R1b(r_R1b), R2b(r_R2b), w1(r_w1), s(r_s), u(r_u) { }
    
    
    // operator overloading of ()  see (ref.1 )
    void operator()( const state_type &x , state_type &dmdt, const double /* t */ )
    {
        
        /* implementation of eq. 5 from Allard, P., Helgstand, M. and Hard, T.,
         Journal of Magnetic Resonance, 129: 19-29 (1997), notice conversion to rad/s
         for some parameters, spin-lock along y-axis. */
        
        // magnetizaitons
        double MAx = x[0];
        double MAy = x[1];
        double MAz = x[2];
        double MBx = x[3];
        double MBy = x[4];
        double MBz = x[5];
        
        // The acutal equations
        
        dmdt[0] = -MAx * R2a - MAy * w_a * 2 * M_PI + w1 * 2 * M_PI * MAz - u * MBx;     // dMAx/dt
        
        dmdt[1] = w_a * 2 * M_PI * MAx - R2a * MAy - u * MBy;                            // dMAy/dt
        
        dmdt[2] = corr_a - w1 * 2 * M_PI * MAx - R1a * MAz - s * MBz;                    // dMAz/dt
        
        dmdt[3] = -u * MAx - R2b * MBx - w_b * 2 * M_PI * MBy + w1 * 2 * M_PI * MBz;     // dMBx/dt
        
        dmdt[4] = -u * MAy + w_b  * 2 * M_PI * MBx - R2b * MBy;                          // dMBy/dt
        
        dmdt[5] = corr_b - s * MAz - w1 * 2 * M_PI * MBx - R1b * MBz;                    // dMBz/dt
        
    }

    
private:
    
    // constants
    double w_a;        // offset spin A relaive to spin lock
    double R1a;
    double R2a;
    
    double w_b;        // offset spin B relaive to spin lock
    double R1b;
    double R2b;
    
    double w1;         // spin lock strength
    double u;          // transverse cross relaxation rate constant
    double s;          // longitsuinal cross relaxation rate constant
    
    // equilibrium magnetizations
    const double MAz_zero = 1.0;
    const double MBz_zero = 1.0;
    
    // some necessary recalculaitons
    double corr_a = R1a * MAz_zero + s * MBz_zero;
    double corr_b = s * MAz_zero + R1b * MBz_zero;
    
};







///////////////////////////////////////////////////////////////////////////////////////////////

class Magnetization_saver
{
    /* A Magentization saver object is the same same as an "observer" in (ref. 1)
     it is called by the integration routine and saves the state of the magnetizaiton. */

    
public:
    
    // Constructor, recieve pointer to array allocated in allard97_wrapper
    Magnetization_saver(double *&r_array_p) : mag_n_time_holder(r_array_p), call_counter(0) {}
    
    // overload ()
    void operator()( const state_type &x, double t)
    {
    
        mag_n_time_holder[call_counter*7+0] = t;
        mag_n_time_holder[call_counter*7+1] = x[0];
        mag_n_time_holder[call_counter*7+2] = x[1];
        mag_n_time_holder[call_counter*7+3] = x[2];
        mag_n_time_holder[call_counter*7+4] = x[3];
        mag_n_time_holder[call_counter*7+5] = x[4];
        mag_n_time_holder[call_counter*7+1] = x[5];
    
        call_counter+=1;
    }


private:
    
    double *mag_n_time_holder;   // pointer to start of array
    int call_counter;
    
};







///////////////////////////////////////////////////////////////////////////////////////////////
// Start of Python module part (ref. 3-5)
///////////////////////////////////////////////////////////////////////////////////////////////



static PyObject* allard97_wrapper(PyObject *dummy, PyObject *args)
{
    
    /* C++ function that recieves python arguments convert them 
     C++ and performs integration and returns a PyObject in this case
     a numpy 2D array  */
    
    
    double r_w_a;        // offset A
    double r_R1a;
    double r_R2a;
    
    double r_w_b;        // offset B
    double r_R1b;
    double r_R2b;
    
    double r_mt;         // mixing time
    double r_w1;         // spin-lock strength
    double r_s;          // longitudinal cross-relaxation constant
    double r_u;          // transverse cross-relaxation constant
    
    double r_max;        // starting magnetizations
    double r_may;
    double r_maz;
    double r_mbx;
    double r_mby;
    double r_mbz;
    
   
    /* PyArg_ParseTuple will extract the values from python and put them in cpp variables above, see (ref. 6) */
    if (!PyArg_ParseTuple(args,"dddddddddddddddd", &r_mt, &r_w_a, &r_R1a, &r_R2a, &r_w_b, &r_R1b, &r_R2b, &r_w1, &r_s, &r_u,
                          &r_max, &r_may, &r_maz, &r_mbx, &r_mby, &r_mbz)) return NULL;
    
    
    /* Create Ode_sys function object with the wnated constants.
    (offset A, R1a, R2a, offset B, R1b, R2b, SL-strength,
     longitudinal cross relaxation, transverse cross relaxation)*/
    Ode_sys allard97(r_w_a, r_R1a, r_R2a, r_w_b, r_R1b, r_R2b, r_w1, r_s, r_u);
    
    
    /* initial value of magentizations at start of spin lock */
    state_type mag_column_vector = {{ r_max , r_may , r_maz, r_mbx , r_mby , r_mbz}};
    
    
    
    
    /* INTEGRATION PART:
    -----------------------
     
    See (ref. 7) which is one Boost example that uses the dopri5 stepper. */
    
    
    // declare stepper type
    typedef runge_kutta_dopri5< state_type > dopri5_type;
    
    // declare integrator type with chosen stepper
    typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
    
    // integrator object to use
    controlled_dopri5_type dopri5;
    
    // starting time point
    double start_time = 0.0;
    
    // step size, dopri5 will ignore and adjust this though
    double dt = 0.001;
    
    // time points for observer calls
    int N = r_mt/0.0001;
    vector<double> ob_calls(N);
    
    for(int i = 0; i <= N; i++)
    {
        ob_calls[i] = 0.0001 * i;
    }
    
    // number of steps integrator will take
    int steps;
    
    // create a 1D array with new
    double *arr = new double[N*7];
    
    // create observer/saver object, notice passing of pointer to
    // newly allocated array to the observer/saver
    Magnetization_saver mag_saver(arr);
    
    // do integration, passing integrator, function object, initial magnetization, start time, end time, step size and observer
    steps = integrate_times( dopri5 ,allard97 ,mag_column_vector , ob_calls.begin(), ob_calls.end(), dt, mag_saver);
    
    
    
    
    /* RETURN NUMPY ARRAY PART:
     -----------------------
     
     A tricky part, see ref 8-11, reference 11 = really good solved the problem I think */

    // dimensions for the numpy array, apparently long int is necesary
    long int dimensions[2];
    dimensions[0] = N;                // num rows
    dimensions[1] = 7;                // num cols
    
    
    // Create the numpy array
    PyObject *evolution_during_mixing = PyArray_SimpleNewFromData(2, dimensions, NPY_FLOAT64, arr);
    
    
    /* Ref. 11 solved my problems, solution = typecast to PyArrayObject and use PyArray_ENABLEFLAGS
     to make the array object own its data and python release the memory after.  */
    PyArray_ENABLEFLAGS((PyArrayObject*)evolution_during_mixing, NPY_ARRAY_OWNDATA);
    
    
    // return the created array
    return evolution_during_mixing;
    
}





///////////////////////////////////////////////////////////////////////////////////////////////

static PyMethodDef NmrMethods[] = {
    
    /* an array of structures each structure defines one function see (ref. 4)
    {"python_function_name", C function name (in this file),
     METH_VARARGS = flag says keyeords arguments not used,
     "short descritive string"}  */
    
    {"allard97_roesy_mixer" , allard97_wrapper,
        METH_VARARGS, "allard97_roesy_mixer describes magnetization evolution during ROESY spin-lock."},
    
    {NULL, NULL, 0, NULL} /* Sentinel = necessary by convention security */
    
};





///////////////////////////////////////////////////////////////////////////////////////////////

static struct PyModuleDef NmrModule = {
    
    // As describes in the docs (ref. 3-5)
    PyModuleDef_HEAD_INIT, "NMR_utils", NULL, -1, NmrMethods
    
};





///////////////////////////////////////////////////////////////////////////////////////////////

PyMODINIT_FUNC // Initiation function

PyInit_NMR_utils(void)
{
    import_array()   // necessary to work with numpy array
    return PyModule_Create(&NmrModule);
}




///////////////////////////////////////////////////////////////////////////////////////////////

/*
REFERECNES:
1. https://www.boost.org/doc/libs/1_60_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/short_example.html
2. https://en.cppreference.com/w/cpp/language/initializer_list
3. https://docs.python.org/3/extending/extending.html
4. https://docs.scipy.org/doc/numpy/reference/c-api.html
5. https://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html
6. https://docs.python.org/3/c-api/arg.html#c.PyArg_ParseTuple
7. https://github.com/headmyshoulder/odeint-v2/blob/master/examples/elliptic_functions.cpp
8. https://docs.scipy.org/doc/numpy-1.14.0/reference/c-api.array.html
9. https://stackoverflow.com/questions/30357115/pyarray-simplenewfromdata-example
10.https://en.wikipedia.org/wiki/Row-_and_column-major_order
11.http://acooke.org/cute/ExampleCod0.html   
*/


