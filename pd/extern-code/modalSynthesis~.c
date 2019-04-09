#include "m_pd.h"
#include "math.h"
#include "stdlib.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif


// TODO:
// why is only the highest frequency playing right now?
// fix setA and setT60 functions to take lists and modify variables correctly

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
    t_float x_f0;       /* modal frequency */

    // arrays for multiple sinusoids
    t_float * x_f;      /* list of modal/fundamental frequencies (length = Nf) */
    t_float * x_ThetaS; /* list of instantaneous phases for sine function (length = Nf) */
    t_float * x_ThetaC; /* list of instantaneous phase for cos function (length = Nf) */
    t_float * x_A;      /* list of initial amplitudes */
    t_float * x_einc;   /* list of envelope wavetable increments */
    t_float * x_eind;   /* list of envelope indices */
    t_int x_Nf;         /* number of modal frequencies */
    t_int x_NA;         /* number of envelope initial amplitudes */
    t_int x_NE;         /* number of envelope increments */
    t_int x_N;          /* holds minimum of x_Nf, x_NA, x_NE */
    


    t_float x_sr;       /* sample rate */
    t_float x_T;        /* sample period (1/sr) */
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

    int i, f;
    float sample;
    float re = 0.0f;
    //float im = 0.0f;
    //float offset = TWO_PI * x->x_f0 * x->x_T;
    float index, frac, i0, i1;
    float envVal;
    float oneOverN = 1.0f / (float)x->x_N;

    /* calculate phase offsets */
    float offset[x->x_Nf];
    for(f = 0; f < x->x_Nf; f++) {
        offset[f] = TWO_PI * x->x_f[f] * x->x_T;
    }

    for(i = 0; i < bufferSize; i++)
    {
        
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
        

        /* for an array of frequencies */
        sample = 0.0f;
        for(f = 0; f < x->x_N; f++) {
            /* calculate cosine value using wavetable */
            // for cosine, we need to do cos = sin(pi/2 - Theta)
        
            index = x->x_ThetaC[f] * WAVETABLE_LENGTH_OVER_TWO_PI;
            i0 = trunc(index);
            frac = index - i0;
            i1 = i0 + 1.0f;
            if(i1 > WAVETABLE_LENGTH)
                i1 = 0;
            re = ((1.0f - frac) * x->x_sinWavetable[(int)i0]) + (frac * x->x_sinWavetable[(int)i1]);


            /* calculate exponentially decaying envelope value using the wavetable */
            if(x->x_eind[f] >= ENV_WAVETABLE_LENGTH) {
                envVal = 0.0f;
                x->x_eind[f] = ENV_WAVETABLE_LENGTH-2;
            } else {
                i0 = trunc(x->x_eind[f]);
                frac = x->x_eind[f] - i0;
                i1 = i0 + 1.0f;
                if(i1 > ENV_WAVETABLE_LENGTH)
                    i1 = ENV_WAVETABLE_LENGTH-1;
                envVal = ((1.0f - frac) * x->x_decayExpWavetable[(int)i0]) + (frac * x->x_decayExpWavetable[(int)i1]);
                x->x_eind[f] += x->x_einc[f];
            }

            //sample = sample + re;
            //sample = x->x_A[f] * envVal * re;
            sample = sample + x->x_A[f] * envVal * re;

            // increment offset and wrap Theta variables
            x->x_ThetaS[f] += offset[f];
            while(x->x_ThetaS[f] < 0)
                x->x_ThetaS[f] += TWO_PI;
            while(x->x_ThetaS[f] >= TWO_PI)
                x->x_ThetaS[f] -= TWO_PI;

            x->x_ThetaC[f] = M_PI_2 - x->x_ThetaS[f];
            while(x->x_ThetaC[f] < 0)
                x->x_ThetaC[f] += TWO_PI;
            while(x->x_ThetaC[f] >= TWO_PI)
                x->x_ThetaC[f] -= TWO_PI;

        }

        *(out + i) = oneOverN * sample;

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

    int f;
    for(f = 0; f < x->x_Nf; f++) {

        // reset envelope index
        x->x_eind[f] = 0.0f;

        // reset oscillator phases
        x->x_ThetaS[f] = 0.0f;
        x->x_ThetaC[f] = 0.0f;
    }
    
}

// set the modal frequencies using a message like "setF0 1.0 0.9 0.8"
void modalSynthesis_setF0(t_modalSynthesis *x,t_symbol *selector, int argcount, t_atom *argvec)
{

    post("setF0: selector %s", selector->s_name);

    if(argcount > x->x_Nf) 
    {
        // allocate space for x->x_f, x->x_ThetaS, x->x_ThetaC
        free(x->x_f);
        free(x->x_ThetaS);
        free(x->x_ThetaC);
        x->x_f = (t_float *)malloc(argcount * sizeof(t_float));
        x->x_ThetaS = (t_float *)malloc(argcount * sizeof(t_float));
        x->x_ThetaC = (t_float *)malloc(argcount * sizeof(t_float));
    }

     x->x_Nf = argcount;

    // copy values into x->x_Nf
    int i;
    for (i = 0; i < x->x_Nf; i++)
    {
        if (argvec[i].a_type == A_FLOAT) {
            x->x_f[i] = argvec[i].a_w.w_float;
            x->x_ThetaS[i] = 0.0f;
            x->x_ThetaC[i] = 0.0f;
        } else if (argvec[i].a_type == A_SYMBOL) {
            post("F0 values must be floats.");
        }
    }

    // adjust x->x_N 
    x->x_N = fmin(x->x_Nf, x->x_NA);
    x->x_N = fmin(x->x_N, x->x_NE);

}


// set the amplitude envelope initial amplitude using a message like "setA 1.0 0.9 0.8"
void modalSynthesis_setA(t_modalSynthesis *x,t_symbol *selector, int argcount, t_atom *argvec)
{

    post("setA: selector %s", selector->s_name);

    if(argcount > x->x_NA) 
    {
        // allocate space for x->x_A
        free(x->x_A);
        x->x_A = (t_float *)malloc(argcount * sizeof(t_float));
    }

     x->x_NA = argcount;
    
   
    int i;
    for (i = 0; i < x->x_NA; i++)
    {
        if (argvec[i].a_type == A_FLOAT) {
            x->x_A[i] = argvec[i].a_w.w_float;
        } else if (argvec[i].a_type == A_SYMBOL) {
            post("A values must be floats.");
        }
    }

    // adjust x->x_N 
    x->x_N = fmin(x->x_Nf, x->x_NA);
    x->x_N = fmin(x->x_N, x->x_NE);
}

// set the T60 times for our amplitude envelopes using a message "setT60 1.0 0.9 0.8"
void modalSynthesis_setT60(t_modalSynthesis *x,t_symbol *selector, int argcount, t_atom *argvec)
{

    post("setT60: selector %s", selector->s_name);

    if(argcount > x->x_NE) 
    {
        // allocate space for x->x_einc and x->x_eind
        free(x->x_einc);
        free(x->x_eind);
        x->x_einc = (t_float *)malloc(argcount * sizeof(t_float));
        x->x_eind = (t_float *)malloc(argcount * sizeof(t_float));
    }

     x->x_NE = argcount;
    
    // figure out all new values for x->x_einc
    int i;
    t_float T60N;
    for (i = 0; i < x->x_NE; i++)
    {
        if (argvec[i].a_type == A_FLOAT) {
            T60N = argvec[i].a_w.w_float * x->x_sr;
            x->x_einc[i] = (float)ENV_WAVETABLE_T60N / T60N;
            //x->x_eind[i] = 0.0f;
        } else if (argvec[i].a_type == A_SYMBOL) {
            post("T60 values must be floats.");
        }
    }

    // adjust x->x_N 
    x->x_N = fmin(x->x_Nf, x->x_NA);
    x->x_N = fmin(x->x_N, x->x_NE);
}


static void *modalSynthesis_new(void)
{
    t_modalSynthesis *x = (t_modalSynthesis *)pd_new(modalSynthesis_class);

    // create outlet for outgoing audio signal
    outlet_new(&x->x_obj, gensym("signal"));

    // initialize internal variables
    //x->x_f0 = 220.0f;
    //x->x_A = 1.0f;
    //x->x_T60 = 0.8f;
    // x->x_ThetaS = 0.0f;
    // x->x_ThetaC = 0.0f;


    x->x_Nf = 7;
    x->x_NA = 7;
    x->x_NE = 7;
    x->x_N = fmin(x->x_Nf, x->x_NA);
    x->x_N = fmin(x->x_N, x->x_NE);

    x->x_f = malloc( sizeof(t_float) * x->x_Nf );
    x->x_ThetaS = malloc( sizeof(t_float) * x->x_Nf );
    x->x_ThetaC = malloc( sizeof(t_float) * x->x_Nf );
    x->x_A = malloc( sizeof(t_float) * x->x_NA );
    x->x_einc = malloc( sizeof(t_float) * x->x_NE );
    x->x_eind = malloc( sizeof(t_float) * x->x_NE );

    int f;
    for(f = 0; f < x->x_N; f++) {
        if(f == 0) {
            x->x_f[f] = 440.0f;
        } else {
            x->x_f[f] = 440.0f * powf(2.0f * f + 3.0f, 2.0f) / powf(3.011f, 2.0f);
        }
        x->x_ThetaS[f] = 0.0f;
        x->x_ThetaC[f] = 0.0f;
        x->x_A[f] = 1.0f;

        // initialize t60 and envelope variables for 1.0sec for now
        t_float T60 = 1.0f;
        t_float T60N = T60 * 44100.0f; // so we have something working by default
        x->x_einc[f] = (float)ENV_WAVETABLE_T60N / T60N;
        //x->x_eind[f] = 0.0f;
        x->x_eind[f] = ENV_WAVETABLE_LENGTH;
    }
    

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
    class_addmethod(modalSynthesis_class, (t_method)modalSynthesis_setF0, gensym("setF0"), A_GIMME, 0);
    class_addmethod(modalSynthesis_class, (t_method)modalSynthesis_setA, gensym("setA"), A_GIMME, 0);
    class_addmethod(modalSynthesis_class, (t_method)modalSynthesis_setT60, gensym("setT60"), A_GIMME, 0);

}


void modalSynthesis_tilde_free(t_modalSynthesis *x)
{
    free(x->x_sinWavetable);
    free(x->x_decayExpWavetable);

    free(x->x_f);
    free(x->x_ThetaS);
    free(x->x_ThetaC);
    free(x->x_A);
    free(x->x_einc);
    free(x->x_eind);
}
