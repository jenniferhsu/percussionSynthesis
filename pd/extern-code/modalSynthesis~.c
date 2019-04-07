#include "m_pd.h"
#include "math.h"
#include "stdlib.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif


/* ------------------------ modalSynthesis~ ----------------------------- */
#define TWO_PI 8.0f * atan(1.0f)
#define INVERSE_TWO_PI 1.0f / (8.0f * atan(1.0f))
#define WAVETABLE_LENGTH 2048   // length of wavetable
#define WAVETABLE_LENGTH_OVER_TWO_PI 2048 * (1.0f / (8.0f * atan(1.0f)))
#define ENV_WAVETABLE_LENGTH 2048          // envelope wavetable length
#define ENV_WAVETABLE_T60N 1536             // envelope wavetable T60 in samples
#define MAXT60SEC 10

/* tilde object to take absolute value. */

static t_class *modalSynthesis_class;

typedef struct _modalSynthesis
{
    t_object x_obj;     /* obligatory header */
    t_float x_f0;
    //t_float * x_f;        /* list of modal/fundamental frequencies */
    //t_int x_Nf;         /* number of modal frequencies */
    t_float x_A;        /* initial amplitude */
    t_float x_T60;      /* desired T60 time for amplitude envelope */
    t_float x_ThetaS;   /* instantaneous phase for sine function */
    t_float x_ThetaC;   /* instantaneous phase for cos function */
    t_float x_sr;       /* sample rate */
    t_float x_T;        /* sample period (1/sr) */
    t_float x_einc;     /* envelope wavetable increment */
    t_float x_eind;       /* envelope index */
    t_float * x_sinWavetable;   /* the sine wavetable to read from */
    t_float * x_decayExpWavetable;   /* the decaying exp wavetable to read from */

} t_modalSynthesis;

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *modalSynthesis_perform(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *out = (t_float *)(w[2]);
    t_modalSynthesis * x = (t_modalSynthesis *)(w[3]);
    int bufferSize = (int)(w[4]);

    int i;
    float sample;
    float re = 0.0f;
    //float im = 0.0f;
    float offset = TWO_PI * x->x_f0 * x->x_T;
    float index, frac, i0, i1;
    float envVal;
    for(i = 0; i < bufferSize; i++)
    {
        
        /* simple sine wave test */

        /* find the indices for the wavetable from the inst. phase
            and linearly interpolate */
        // use this for sine:
        /* 
        index = x->x_ThetaS * WAVETABLE_LENGTH_OVER_TWO_PI;
        i0 = trunc(index);
        frac = index - i0;
        i1 = i0 + 1.0f;
        if(i1 > WAVETABLE_LENGTH)
            i1 = 0;
        im = ((1.0f - frac) * x->x_sinWavetable[(int)i0]) + (frac * x->x_sinWavetable[(int)i1]);
        */
        
        /* calculate cosine value using wavetable */
        // for cosine, we need to do cos = sin(pi/2 - Theta)
        index = x->x_ThetaC * WAVETABLE_LENGTH_OVER_TWO_PI;
        i0 = trunc(index);
        frac = index - i0;
        i1 = i0 + 1.0f;
        if(i1 > WAVETABLE_LENGTH)
            i1 = 0;
        re = ((1.0f - frac) * x->x_sinWavetable[(int)i0]) + (frac * x->x_sinWavetable[(int)i1]);

        /* calculate exponentially decaying envelope value using the wavetable */
        if(x->x_eind >= ENV_WAVETABLE_LENGTH) {
            envVal = 0.0f;
            x->x_eind = ENV_WAVETABLE_LENGTH-2;
        } else {
            i0 = trunc(x->x_eind);
            frac = x->x_eind - i0;
            i1 = i0 + 1.0f;
            if(i1 > ENV_WAVETABLE_LENGTH)
                i1 = ENV_WAVETABLE_LENGTH-1;
            envVal = ((1.0f - frac) * x->x_decayExpWavetable[(int)i0]) + (frac * x->x_decayExpWavetable[(int)i1]);
            x->x_eind += x->x_einc;
        }


        sample = x->x_A * envVal * re;
        //sample = x->x_A * re;
        //sample = re;
        *(out + i) = sample;

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
static void modalSynthesis_dsp(t_modalSynthesis *x, t_signal **sp)
{
    
    x->x_sr = sp[0]->s_sr;
    x->x_T = 1.0f / x->x_sr;

    dsp_add(modalSynthesis_perform, 4, sp[0]->s_vec, sp[1]->s_vec, x, sp[0]->s_n);
}

void modalSynthesis_bang(t_modalSynthesis *x)
{
    x->x_eind = 0.0f;
}


// right inlet callback to set the Feedback FM coefficient
void modalSynthesis_setA(t_modalSynthesis *x, t_floatarg f)
{
    x->x_A = f;
}

// right inlet callback to set the Feedback FM coefficient
void modalSynthesis_setT60(t_modalSynthesis *x, t_floatarg f)
{
    x->x_T60 = f;
    t_float T60N = x->x_T60 * x->x_sr;

    x->x_einc = (float)ENV_WAVETABLE_T60N / T60N;
}


static void *modalSynthesis_new(void)
{
    t_modalSynthesis *x = (t_modalSynthesis *)pd_new(modalSynthesis_class);
    // create right inlet for inital amplitude A
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("A"));
    // create right inlet for T60
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("T60"));
    // create outlet for outgoing audio signal
    outlet_new(&x->x_obj, gensym("signal"));

    // initialize internal variables
    x->x_f0 = 220.0f;
    x->x_A = 1.0f;
    x->x_T60 = 0.8f;
    x->x_ThetaS = 0.0f;
    x->x_ThetaC = 0.0f;

    /*
    x->x_Nf = 7;
    x->x_fi = malloc( sizeof(t_float) * x->x_Nf );
    int i;
    for(i = 0; i < x->x_Nf; i++) {
        if(i == 0) {
            x->x_f[i] = 440.0f;
        } else {
            x->x_f[i] = 440.0f * (2.0f * i + 3)^2 / 3.011^2;
        }
        post("%d: %f\n", i, x->x_f[i]);
    }
    */

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
        //fscanf(myFile, "%f,", &x->x_sinWavetable[i] );
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
void modalSynthesis_tilde_setup(void)
{
    modalSynthesis_class = class_new(gensym("modalSynthesis~"), (t_newmethod)modalSynthesis_new, 0,
        sizeof(t_modalSynthesis), 0, A_DEFFLOAT, 0);
        /* this is magic to declare that the leftmost, "main" inlet
        takes signals; other signal inlets are done differently... */
        // fundamental frequency is set using the left inlet
    CLASS_MAINSIGNALIN(modalSynthesis_class, t_modalSynthesis, x_f0);
        /* here we tell Pd about the "dsp" method, which is called back
    when DSP is turned on. */
    class_addmethod(modalSynthesis_class, (t_method)modalSynthesis_dsp, gensym("dsp"), 0);
    class_addbang(modalSynthesis_class, modalSynthesis_bang);
    
    // register callback methods for right inlets
    class_addmethod(modalSynthesis_class, (t_method)modalSynthesis_setA, gensym("A"), A_FLOAT, 0);
    class_addmethod(modalSynthesis_class, (t_method)modalSynthesis_setT60, gensym("T60"), A_FLOAT, 0);


}


void modalSynthesis_tilde_free(t_modalSynthesis *x)
{
    free(x->x_sinWavetable);
    free(x->x_decayExpWavetable);
}
