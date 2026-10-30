"""The terms of the toy model, one file a term (the owner's message 119: "Are you ready to add further terms and details
to the model?"). The line is the Fourier transform of a function C(t). A term gives one of two things, or both:

    width(T_celsius, values, const)        a Lorentzian width in MHz, which adds to those of the other terms: together
                                           they are the factor exp(-pi Gamma |t|) of C(t)
    factor(t, T_celsius, values, const)    one more factor of C(t), on the times t in microseconds

`const.toml` names each term with the number of the detail at which it enters, under `details`, and `model.py` takes the
terms of a detail from there. So a new term is a new file here and a line of each of the two files: no program of the
model changes, and the kept fits of the details before it stay as they are. No number of the physics is written in a term."""
