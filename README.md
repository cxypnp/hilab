# hilab
Code for the HI-LAB project on the event plane, UrQMD, etc.

# UrQMD
To run UrQMD on latest gFortran complier, before "make":
Line 3733 and 3744 of make22.f should be:

    integer io,i,im,ip,ir,xmax
    real*8 e,x(minnuc:maxmes)    

The type of xmax is not correct for new gFortran.