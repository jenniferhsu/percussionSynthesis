#include "m_pd.h"
#include "math.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ feedbackFM~ ----------------------------- */

/* tilde object to take absolute value. */

static t_class *feedbackFM_class;

typedef struct _feedbackFM
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f;    	/* place to hold inlet's value if it's set by message */
    t_float x_B;        /* feedback FM coefficient */
    t_float x_sr;       /* sample rate */
    t_float x_T;        /* sample period (1/sr) */
    t_int x_cnt;        /* sample counter */

    t_float x_reMinus1;  /* Real(h(n-1)) */
    t_float x_imMinus1;  /* Imag(h(n-1)) */

} t_feedbackFM;

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *feedbackFM_perform(t_int *w)
{
    t_float *in = (t_float *)(w[1]);
    t_float *out = (t_float *)(w[2]);
    t_feedbackFM * x = (t_feedbackFM *)(w[3]);
    int bufferSize = (int)(w[4]);

    int i;
    float sample;
    float re = 0.0f;
    float im = 0.0f;
    float c = 0.0f;
    float s = 0.0f;
    float theta = 0.0f;
    for(i = 0; i < bufferSize; i++)
    {
        // microphone test
        /*
        sample = *(in + i);
        *(out + i) = sample;
        */


        // simple sine wave test
        /*
        theta = 2.0f * M_PI * x->x_f * x->x_cnt * x->x_T;
        re = cos(theta);
        im = sin(theta);
        sample = re;
        *(out + i) = sample;

        x->x_cnt++;
        if(x->x_cnt > x->x_sr)
            x->x_cnt = 0;
        */

        // feedback FM
        theta = 2.0f * M_PI * x->x_f * x->x_T * (1.0f + (x->x_B * x->x_reMinus1));
        c = cos(theta);
        s = sin(theta);

        // rotation multiplication
        re = *(in + i) + (c * x->x_reMinus1) - (s * x->x_imMinus1);
        im = (c * x->x_imMinus1) + (s * x->x_reMinus1);

        // output sample (use the real value)
        sample = re;
        *(out + i) = sample;

        // save values for next iteration
        x->x_reMinus1 = re;
        x->x_imMinus1 = im;

    }

    return (w+5);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */
static void feedbackFM_dsp(t_feedbackFM *x, t_signal **sp)
{
    
    x->x_sr = sp[0]->s_sr;
    x->x_T = 1.0f / x->x_sr;

    dsp_add(feedbackFM_perform, 4, sp[0]->s_vec, sp[1]->s_vec, x, sp[0]->s_n);
}


// right inlet callback to set the Feedback FM coefficient
void feedbackFM_setB(t_feedbackFM *x, t_floatarg f)
{
    x->x_B = f;
}

static void *feedbackFM_new(void)
{
    t_feedbackFM *x = (t_feedbackFM *)pd_new(feedbackFM_class);
    // create right inlet for Feedback FM coefficient B
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("B"));
    // create outlet for outgoing audio signal
    outlet_new(&x->x_obj, gensym("signal"));

    //x->x_f = 0;
    // these will be inlets so the user can set them
    x->x_f = 698.0f;
    x->x_B = 0.9f;
    
    // next step: set input to be an impulse to set it off

    // initialize internal variables
    //x->x_reMinus1 = 1.0f;
    x->x_reMinus1 = 0.0f;
    x->x_imMinus1 = 0.0f;

    return (x);
}

    /* this routine, which must have exactly this name (with the "~" replaced
    by "_tilde) is called when the code is first loaded, and tells Pd how
    to build the "class". */
void feedbackFM_tilde_setup(void)
{
    feedbackFM_class = class_new(gensym("feedbackFM~"), (t_newmethod)feedbackFM_new, 0,
    	sizeof(t_feedbackFM), 0, A_DEFFLOAT, 0);
	    /* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(feedbackFM_class, t_feedbackFM, x_f);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(feedbackFM_class, (t_method)feedbackFM_dsp, gensym("dsp"), 0);
    
    // register callback method for right inlet float
    class_addmethod(feedbackFM_class, (t_method)feedbackFM_setB, gensym("B"), A_FLOAT, 0);
}
