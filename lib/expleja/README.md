#`expleja` - **The action of the matrix exponential function**

##About
`expleja` contains the MATLAB function 'expleja` for directly computing expm(hA)V, where A is an n-by-n matrix and V is a tall, skinny an n-by-m matrix.

The basic functionality is as follows:

* expleja(t,A,b) computes expm(t*A)b
* expleja(t,A,B) computes expm(t*A)B, for a skinny n-by-m matrix B

The function works for any matrix A and only use matrix-vector products with A and A^*.

Function `demo_expleja.m` runs some simple test code. In Octave you can also use `demo expleja`. 

The description of the complete interface can be found in `expleja.m` by typing `help expleja` in MATLAB or OCTAVE as well as with `Contents.m`.

Detail on the underlying algorithms can be found in:

M. Caliari, P. Kandolf, A. Ostermann, S. Rainer, "[The Leja method revisited: backward error analysis for the matrix exponential](http://epubs.siam.org/doi/abs/10.1137/15M1027620)" SIAM Journal on Scientific Computation, 38 (3), A1639-A1661(2016) 

auxiliary/normest1.m is needed for GNU Octave <= 4.0.x
