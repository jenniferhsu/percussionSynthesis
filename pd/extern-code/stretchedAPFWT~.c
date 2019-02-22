#include "m_pd.h"
#include "math.h"
#include "stdlib.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ stretchedAPFWT~ ----------------------------- */
#define WAVETABLE_LENGTH 2048   // length of wavetable


/* tilde object to take absolute value. */

static t_class *stretchedAPFWT_class;

typedef struct _stretchedAPFWT
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f0;    	/* sounding frequency */
    t_float x_b0;       /* stretched APF coefficient/timbre control */
    t_float x_sr;       /* sample rate */
    t_float x_T;        /* sample period (1/sr) */
    t_float x_n;        /* sample counter */
} t_stretchedAPFWT;

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *stretchedAPFWT_perform(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *out = (t_float *)(w[2]);
    t_stretchedAPFWT * x = (t_stretchedAPFWT *)(w[3]);
    int bufferSize = (int)(w[4]);

    int i;
    float theta = 0.0f;
    float s;
    float c;
    float angleH = 0.0f;
    float realH = 0.0f;
    float samplePeriod = x->x_sr * (1.0f/x->x_f0);
    for(i = 0; i < bufferSize; i++)
    {

        theta = 2.0f * M_PI * x->x_f0 * x->x_n * x->x_T;
        angleH = theta - (2.0f * atan(x->x_b0 * sin(theta) / (1.0f + x->x_b0 * cos(theta))));
        realH = cos(angleH);

        // output sample (use the real value)
        *(out + i) = realH;

        // increment and wrap sample counter
        x->x_n += 1.0f;
        if(x->x_n >= samplePeriod)
        {
            //post("%f", x->x_n);
            x->x_n = 0.0f;
        }

    }

    return (w+5);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */
static void stretchedAPFWT_dsp(t_stretchedAPFWT *x, t_signal **sp)
{
    
    x->x_sr = sp[0]->s_sr;
    x->x_T = 1.0f / x->x_sr;

    dsp_add(stretchedAPFWT_perform, 4, sp[0]->s_vec, sp[1]->s_vec, x, sp[0]->s_n);
}


// right inlet callback to set the stretched APF coefficient
void stretchedAPFWT_setb0(t_stretchedAPFWT *x, t_floatarg f)
{
    x->x_b0 = f;
}


static void *stretchedAPFWT_new(void)
{
    t_stretchedAPFWT*x = (t_stretchedAPFWT *)pd_new(stretchedAPFWT_class);
    // create right inlet for stretched APF coefficient b0
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("b0"));
    // create outlet for outgoing audio signal
    outlet_new(&x->x_obj, gensym("signal"));

    x->x_f0 = 304.0f;
    x->x_b0 = 0.9f;

    // initialize internal variables
    x->x_n = 0.0f;

    // read in sine wavetable
    FILE *myFile;
    myFile = fopen("/Applications/Pd-0.48-0.app/Contents/Resources/extra/sinWavetable.txt", "r");
    if (myFile == NULL){
        post("Error Reading File\n");
        exit(0);
    }

    // read file into array
    float sinWavetable[WAVETABLE_LENGTH];
    int i;
    for (i = 0; i < WAVETABLE_LENGTH; i++) {
        fscanf(myFile, "%f,", &sinWavetable[i] );
    }
    // print out wavetable contents
    // post("wavetable is: \n");
    // for (i = 0; i < WAVETABLE_LENGTH; i++) {
    //     post("%f ", sinWavetable[i]);
    // }
    fclose(myFile);


    return (x);
}

    /* this routine, which must have exactly this name (with the "~" replaced
    by "_tilde) is called when the code is first loaded, and tells Pd how
    to build the "class". */
void stretchedAPFWT_tilde_setup(void)
{
    stretchedAPFWT_class = class_new(gensym("stretchedAPFWT~"), (t_newmethod)stretchedAPFWT_new, 0,
    	sizeof(t_stretchedAPFWT), 0, A_DEFFLOAT, 0);
	    /* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(stretchedAPFWT_class, t_stretchedAPFWT, x_f0);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(stretchedAPFWT_class, (t_method)stretchedAPFWT_dsp, gensym("dsp"), 0);
    
    // register callback method for right inlet float
    class_addmethod(stretchedAPFWT_class, (t_method)stretchedAPFWT_setb0, gensym("b0"), A_FLOAT, 0);

}
