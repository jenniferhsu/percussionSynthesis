#include "m_pd.h"
#include "math.h"
#include "stdlib.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ stretchedAPFMS~ ----------------------------- */
#define TWO_PI 8.0f * atan(1.0f)
#define WAVETABLE_LENGTH 2048   // length of wavetable
#define WAVETABLE_LENGTH_OVER_TWO_PI 2048 * (1.0f / (8.0f * atan(1.0f)))


/* tilde object to take absolute value. */

static t_class *stretchedAPFMS_class;

typedef struct _stretchedAPFMS
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f0;    	/* sounding frequency */
    t_float x_Theta;    /* oscillator phase */
    t_float x_b0;       /* stretched APF coefficient/timbre control */
    t_float x_sr;       /* sample rate */
    t_float x_T;        /* sample period (1/sr) */
    t_float x_sinWavetable[WAVETABLE_LENGTH];   /* the sine wavetable to read from */
} t_stretchedAPFMS;

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *stretchedAPFMS_perform(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *out = (t_float *)(w[2]);
    t_stretchedAPFMS * x = (t_stretchedAPFMS *)(w[3]);
    int bufferSize = (int)(w[4]);

    int i;
    float s;
    float c;
    float offset = TWO_PI * x->x_f0 * x->x_T;
    float angleH = 0.0f;
    float realH = 0.0f;
    float samplePeriod = x->x_sr * (1.0f/x->x_f0);
    for(i = 0; i < bufferSize; i++)
    {

        //theta = 2.0f * M_PI * x->x_f0 * x->x_n * x->x_T;
        x->x_Theta = x->x_Theta + offset;
        angleH = x->x_Theta - (2.0f * atan(x->x_b0 * sin(x->x_Theta) / (1.0f + x->x_b0 * cos(x->x_Theta))));
        realH = cos(angleH);

        // output sample (use the real value)
        *(out + i) = realH;

        // increment and wrap theta
        x->x_Theta += offset;
        while(x->x_Theta < 0)
            x->x_Theta += TWO_PI;
        while(x->x_Theta >= TWO_PI)
            x->x_Theta -= TWO_PI;

    }

    return (w+5);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */
static void stretchedAPFMS_dsp(t_stretchedAPFMS *x, t_signal **sp)
{
    
    x->x_sr = sp[0]->s_sr;
    x->x_T = 1.0f / x->x_sr;

    dsp_add(stretchedAPFMS_perform, 4, sp[0]->s_vec, sp[1]->s_vec, x, sp[0]->s_n);
}


// right inlet callback to set the stretched APF coefficient
void stretchedAPFMS_setb0(t_stretchedAPFMS *x, t_floatarg f)
{
    x->x_b0 = f;
}


static void *stretchedAPFMS_new(void)
{
    t_stretchedAPFMS*x = (t_stretchedAPFMS *)pd_new(stretchedAPFMS_class);
    // create right inlet for stretched APF coefficient b0
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("b0"));
    // create outlet for outgoing audio signal
    outlet_new(&x->x_obj, gensym("signal"));

    // initialize internal variables
    x->x_f0 = 304.0f;
    x->x_Theta = 0.0f;
    x->x_b0 = 0.9f;

    // read in sine wavetable
    FILE *myFile;
    myFile = fopen("/Applications/Pd-0.48-0.app/Contents/Resources/extra/sinWavetable.txt", "r");
    if (myFile == NULL){
        post("Error Reading File\n");
        exit(0);
    }

    // read file into array
    int i;
    for (i = 0; i < WAVETABLE_LENGTH; i++) {
        fscanf(myFile, "%f,", &x->x_sinWavetable[i] );
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
void stretchedAPFMS_tilde_setup(void)
{
    stretchedAPFMS_class = class_new(gensym("stretchedAPFMS~"), (t_newmethod)stretchedAPFMS_new, 0,
    	sizeof(t_stretchedAPFMS), 0, A_DEFFLOAT, 0);
	    /* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(stretchedAPFMS_class, t_stretchedAPFMS, x_f0);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(stretchedAPFMS_class, (t_method)stretchedAPFMS_dsp, gensym("dsp"), 0);
    
    // register callback method for right inlet float
    class_addmethod(stretchedAPFMS_class, (t_method)stretchedAPFMS_setb0, gensym("b0"), A_FLOAT, 0);

}
