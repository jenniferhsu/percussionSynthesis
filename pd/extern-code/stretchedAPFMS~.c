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
#define ENV_WAVETABLE_LENGTH 2048          // envelope wavetable length
#define ENV_WAVETABLE_T60N 1536             // envelope wavetable T60 in samples


/* tilde object to take absolute value. */

static t_class *stretchedAPFMS_class;

typedef struct _stretchedAPFMS
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f0;    	/* sounding frequency */
    t_float x_ThetaS;   /* oscillator phase for sin*/
    t_float x_ThetaC;   /* oscillator phase for cos*/
    t_float x_b0;       /* stretched APF coefficient/timbre control */

    t_float x_A;        /* list of initial amplitudes */
    t_float x_einc;     /* list of envelope wavetable increments */
    t_float x_eind;     /* list of envelope indices */

    t_float x_sr;       /* sample rate */
    t_float x_T;        /* sample period (1/sr) */

    t_float * x_sinWavetable;   /* the sine wavetable to read from */
    t_float * x_decayExpWavetable;   /* the decaying exp wavetable to read from */
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
    float offset = TWO_PI * x->x_f0 * x->x_T;
    float index, frac, i0, i1;
    float sinTheta, cosTheta;
    float angleH = 0.0f;
    float realH = 0.0f;
    //float samplePeriod = x->x_sr * (1.0f/x->x_f0);
    for(i = 0; i < bufferSize; i++)
    {

        /* sin(x->xTheta) */
        index = x->x_ThetaS * WAVETABLE_LENGTH_OVER_TWO_PI;
        i0 = trunc(index);
        frac = index - i0;
        i1 = i0 + 1.0f;
        if(i1 > WAVETABLE_LENGTH)
            i1 = 0;
        sinTheta = ((1.0f - frac) * x->x_sinWavetable[(int)i0]) + (frac * x->x_sinWavetable[(int)i1]);

        /* cos(x->xTheta) */
        index = x->x_ThetaC * WAVETABLE_LENGTH_OVER_TWO_PI;
        i0 = trunc(index);
        frac = index - i0;
        i1 = i0 + 1.0f;
        if(i1 > WAVETABLE_LENGTH)
            i1 = 0;
        cosTheta = ((1.0f - frac) * x->x_sinWavetable[(int)i0]) + (frac * x->x_sinWavetable[(int)i1]);

        angleH = x->x_ThetaS - (2.0f * atan(x->x_b0 * sinTheta / (1.0f + x->x_b0 * cosTheta)));
        realH = cos(angleH);

        // output sample (use the real value)
        *(out + i) = realH;

        // increment offset and wrap Theta variables
        x->x_ThetaS += offset;
        while(x->x_ThetaS < 0)
            x->x_ThetaS += TWO_PI;
        while(x->x_ThetaS >= TWO_PI)
            x->x_ThetaS -= TWO_PI;

        x->x_ThetaC = M_PI_2 - x->x_ThetaS;
        while(x->x_ThetaC < 0)
            x->x_ThetaC += TWO_PI;
        while(x->x_ThetaC >= TWO_PI)
            x->x_ThetaC -= TWO_PI;

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

void stretchedAPFMS_bang(t_stretchedAPFMS *x)
{

    // reset envelope index
    x->x_eind = 0.0f;

    // reset oscillator phase
    x->x_ThetaS = 0.0f;
    x->x_ThetaC = 0.0f;
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
    x->x_ThetaS = 0.0f;
    x->x_ThetaC = 0.0f;
    x->x_b0 = 0.9f;

    x->x_A = 1.0f;
    // initialize t60 and envelope variables for 1.0sec for now
    t_float T60 = 1.0f;
    t_float T60N = T60 * 44100.0f; // so we have something working by default
    x->x_einc = (float)ENV_WAVETABLE_T60N / T60N;
    x->x_eind = 0.0f;

    /* sine wavetable */
    // read in sine wavetable
    FILE *myFile;
    myFile = fopen("/Applications/Pd-0.48-0.app/Contents/Resources/extra/sinWavetable.txt", "r");
    if (myFile == NULL){
        post("Error Reading File\n");
        exit(0);
    }
    // read file into array
    x->x_sinWavetable = (t_float *)malloc(WAVETABLE_LENGTH * sizeof(t_float)); 
    int i;
    for (i = 0; i < WAVETABLE_LENGTH; i++) {
        fscanf(myFile, "%f", &x->x_sinWavetable[i] );
    }
    // print out wavetable contents
    // post("wavetable is: \n");
    // for (i = 0; i < WAVETABLE_LENGTH; i++) {
    //     post("%f ", x->x_sinWavetable[i]);
    // }
    fclose(myFile);

    /* decaying exponential wavetable */
    // read in decaying exponential wavetable
    FILE *decayExpWavetableFile;
    decayExpWavetableFile = fopen("/Applications/Pd-0.48-0.app/Contents/Resources/extra/decayExpWavetable.txt", "r");
    if (decayExpWavetableFile == NULL){
        post("Error Reading File\n");
        exit(0);
    }
    // read file into array
    x->x_decayExpWavetable = (t_float *)malloc(ENV_WAVETABLE_LENGTH * sizeof(t_float)); 
    for (i = 0; i < ENV_WAVETABLE_LENGTH; i++) {
        fscanf(decayExpWavetableFile, "%f", &x->x_decayExpWavetable[i] );
    }
    fclose(decayExpWavetableFile);


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
    class_addbang(stretchedAPFMS_class, stretchedAPFMS_bang);
    
    // register callback method for right inlet float
    class_addmethod(stretchedAPFMS_class, (t_method)stretchedAPFMS_setb0, gensym("b0"), A_FLOAT, 0);

}
